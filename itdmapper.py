#!/usr/bin/env python3
"""FLT3 ITD caller for indexed, coordinate-sorted BBmerged BAM files."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
import tempfile
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from itdmapper_settings import load_settings, validate_settings

try:
    import pysam
except ImportError:  # pragma: no cover - helper tests can run without pysam
    pysam = None


TRANSCRIPT = "NM_004119.3"
TARGET_C_START = 1705
TARGET_C_END = 1942
_PACKAGED_SETTINGS = load_settings()
DEFAULT_CLUSTER_TOLERANCE = int(_PACKAGED_SETTINGS["amplicon"]["endpoint_tolerance"])
DEFAULT_MIN_CLUSTER_FRACTION = float(_PACKAGED_SETTINGS["amplicon"]["min_cluster_fraction"])
DEFAULT_MIN_CLUSTER_SUPPORT = int(_PACKAGED_SETTINGS["amplicon"]["min_cluster_support"])
DEFAULT_END_ANCHOR_TOLERANCE = int(_PACKAGED_SETTINGS["amplicon"]["end_anchor_tolerance"])
DEFAULT_MIN_OVERLAP_FRACTION = float(_PACKAGED_SETTINGS["amplicon"]["min_overlap_fraction"])
DEFAULT_MIN_MAPQ = int(_PACKAGED_SETTINGS["reads"]["min_mapq"])
DEFAULT_MIN_MEAN_BASE_QUALITY = float(_PACKAGED_SETTINGS["reads"]["min_mean_base_quality"])
DEFAULT_EXCLUDE_DUPLICATES = bool(_PACKAGED_SETTINGS["reads"]["exclude_duplicates"])
DEFAULT_EXCLUDE_QCFAIL = bool(_PACKAGED_SETTINGS["reads"]["exclude_qcfail"])


@dataclass(frozen=True)
class Amplicon:
    start: int  # 0-based genomic, inclusive
    end: int    # 0-based genomic, exclusive
    support: int

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def label(self) -> str:
        return f"{self.start + 1}-{self.end}"


@dataclass
class Flt3Reference:
    chrom: str
    genomic_start: int  # 1-based inclusive
    genomic_end: int    # 1-based inclusive
    strand: str
    sequence: str       # stored in FLT3/transcript orientation
    c_by_index: list[Optional[int]]
    index_by_c: dict[int, int]
    region_by_index: list[str]

    @classmethod
    def load(cls, path: Path) -> "Flt3Reference":
        with path.open() as handle:
            header = handle.readline().strip()
            sequence = "".join(line.strip() for line in handle if line.strip())

        if not header.startswith(">"):
            raise ValueError(f"{path} is not a FASTA file")
        match = re.search(
            r"range=([^: ]+):(\d+)-(\d+).*\bstrand=([+-])",
            header,
        )
        if not match:
            raise ValueError(
                "FLT3 FASTA header must contain range=chrom:start-end and strand=+/-"
            )

        chrom, start, end, strand = match.groups()
        genomic_start, genomic_end = int(start), int(end)
        expected = genomic_end - genomic_start + 1
        if len(sequence) != expected:
            raise ValueError(
                f"FLT3 FASTA length {len(sequence)} does not match header span {expected}"
            )

        c_by_index: list[Optional[int]] = [None] * len(sequence)
        index_by_c: dict[int, int] = {}
        region_by_index: list[str] = ["other"] * len(sequence)
        c_pos = 0
        exon_number = 0
        in_coding_exon = False
        for idx, base in enumerate(sequence):
            if base.isupper():
                if not in_coding_exon:
                    exon_number += 1
                    in_coding_exon = True
                c_pos += 1
                c_by_index[idx] = c_pos
                index_by_c[c_pos] = idx
                region_by_index[idx] = f"exon{exon_number}"
            else:
                if in_coding_exon:
                    in_coding_exon = False
                if exon_number > 0:
                    region_by_index[idx] = f"intron{exon_number}"

        for required in (TARGET_C_START, TARGET_C_END):
            if required not in index_by_c:
                raise ValueError(f"CDS coordinate c.{required} not present in {path}")

        return cls(
            chrom=chrom,
            genomic_start=genomic_start,
            genomic_end=genomic_end,
            strand=strand,
            sequence=sequence,
            c_by_index=c_by_index,
            index_by_c=index_by_c,
            region_by_index=region_by_index,
        )

    def genomic_pos_for_index(self, idx: int) -> int:
        if self.strand == "-":
            return self.genomic_end - idx
        return self.genomic_start + idx

    def genomic_pos_for_c(self, c_pos: int) -> int:
        return self.genomic_pos_for_index(self.index_by_c[c_pos])

    def sequence_index_for_genomic_pos(self, genomic_pos: int) -> int:
        if not self.genomic_start <= genomic_pos <= self.genomic_end:
            raise ValueError(f"Genomic position {genomic_pos} lies outside FLT3.fa")
        if self.strand == "-":
            return self.genomic_end - genomic_pos
        return genomic_pos - self.genomic_start

    def target_fetch_interval(self) -> tuple[int, int]:
        """Return hg19 target interval as 0-based half-open coordinates."""
        g1 = self.genomic_pos_for_c(TARGET_C_START)
        g2 = self.genomic_pos_for_c(TARGET_C_END)
        return min(g1, g2) - 1, max(g1, g2)

    def genomic_interval_for_region(self, region: str) -> tuple[int, int]:
        positions = [
            self.genomic_pos_for_index(idx)
            for idx, value in enumerate(self.region_by_index)
            if value == region
        ]
        if not positions:
            raise ValueError(f"Region {region!r} is not present in FLT3.fa")
        return min(positions) - 1, max(positions)

    def extract_interval(self, start0: int, end0: int) -> str:
        """Extract genomic interval in FLT3/transcript orientation."""
        low1 = start0 + 1
        high1 = end0
        if self.strand == "-":
            first_idx = self.sequence_index_for_genomic_pos(high1)
            last_idx = self.sequence_index_for_genomic_pos(low1)
        else:
            first_idx = self.sequence_index_for_genomic_pos(low1)
            last_idx = self.sequence_index_for_genomic_pos(high1)
        return self.sequence[first_idx:last_idx + 1]

    def annotation_rows(self, start0: int, end0: int) -> list[dict[str, object]]:
        """Build mergeITD annotation rows for an inferred amplicon."""
        rows: list[dict[str, object]] = []
        positions: Iterable[int]
        if self.strand == "-":
            positions = range(end0, start0, -1)  # 1-based positions: end0 .. start0+1
        else:
            positions = range(start0 + 1, end0 + 1)

        for amplicon_bp, genomic_pos in enumerate(positions, start=1):
            idx = self.sequence_index_for_genomic_pos(genomic_pos)
            c_pos = self.c_by_index[idx]
            rows.append(
                {
                    "amplicon_bp": amplicon_bp,
                    "region": self.region_by_index[idx],
                    "chr13_bp": genomic_pos,
                    "transcript_bp": c_pos if c_pos is not None else ".",
                    "protein_as": ((c_pos - 1) // 3 + 1) if c_pos is not None else ".",
                }
            )
        return rows


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def resolve_bam_contig(bam, reference_chrom: str) -> str:
    candidates = [reference_chrom]
    if reference_chrom.startswith("chr"):
        candidates.append(reference_chrom[3:])
    else:
        candidates.append(f"chr{reference_chrom}")
    for candidate in candidates:
        if candidate in bam.references:
            return candidate
    raise ValueError(
        f"BAM does not contain {reference_chrom!r} (or its chr/no-chr equivalent)"
    )


def validate_bam(bam_path: Path):
    if pysam is None:
        raise RuntimeError("pysam is required: pip install pysam")
    if not bam_path.is_file():
        raise ValueError(f"BAM not found: {bam_path}")
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        sort_order = bam.header.to_dict().get("HD", {}).get("SO")
        if sort_order != "coordinate":
            raise ValueError(
                f"BAM must declare coordinate sorting (HD:SO=coordinate); found {sort_order!r}"
            )
        if not bam.has_index():
            raise ValueError("BAM must have a readable .bai or .csi index")
    except Exception:
        bam.close()
        raise
    return bam


def _mean_base_quality(read) -> Optional[float]:
    qualities = read.query_qualities
    if qualities is None or len(qualities) == 0:
        return None
    return sum(qualities) / len(qualities)


def primary_target_reads(
    bam,
    contig: str,
    start0: int,
    end0: int,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_mean_base_quality: float = DEFAULT_MIN_MEAN_BASE_QUALITY,
    exclude_duplicates: bool = DEFAULT_EXCLUDE_DUPLICATES,
    exclude_qcfail: bool = DEFAULT_EXCLUDE_QCFAIL,
):
    candidates = []
    for read in bam.fetch(contig, start0, end0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None or read.reference_end is None:
            continue
        candidates.append(read)

    if not candidates:
        raise ValueError("No primary reads overlap FLT3 exons 14-15")
    if any(read.is_paired for read in candidates):
        raise ValueError(
            "Input BAM contains paired-end records; ITDmapper expects already-BBmerged single reads"
        )

    reads = []
    for read in candidates:
        if min_mapq > 0 and read.mapping_quality < min_mapq:
            continue
        if exclude_duplicates and read.is_duplicate:
            continue
        if exclude_qcfail and read.is_qcfail:
            continue
        if min_mean_base_quality > 0:
            mean_quality = _mean_base_quality(read)
            if mean_quality is None or mean_quality < min_mean_base_quality:
                continue
        reads.append(read)

    if not reads:
        raise ValueError("No primary reads pass configured read-quality filters")
    return reads


def _cluster_endpoint_pairs(
    pair_counts: Counter[tuple[int, int]], tolerance: int
) -> list[dict[str, object]]:
    clusters: list[dict[str, object]] = []
    for pair, count in pair_counts.most_common():
        matched = None
        for cluster in clusters:
            rep_start, rep_end = cluster["representative"]
            if abs(pair[0] - rep_start) <= tolerance and abs(pair[1] - rep_end) <= tolerance:
                matched = cluster
                break
        if matched is None:
            clusters.append(
                {
                    "representative": pair,
                    "support": count,
                    "best_exact_count": count,
                }
            )
        else:
            matched["support"] += count
            if count > matched["best_exact_count"]:
                matched["representative"] = pair
                matched["best_exact_count"] = count
    return clusters


def infer_amplicons(
    reads,
    exon14_interval: tuple[int, int],
    exon15_interval: tuple[int, int],
    tolerance: int = DEFAULT_CLUSTER_TOLERANCE,
    min_cluster_fraction: float = DEFAULT_MIN_CLUSTER_FRACTION,
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT,
) -> list[Amplicon]:
    """Infer amplicons from dominant genomic alignment endpoint pairs.

    Primary reads overlapping both exon 14 and exon 15 are used as anchors so
    shorter exon14-to-exon15 amplicons remain eligible while ITD/soft-clipped
    outliers cannot expand the WT reference. Near-identical endpoint pairs are
    clustered to tolerate small alignment shifts.
    """
    exon14_start, exon14_end = exon14_interval
    exon15_start, exon15_end = exon15_interval

    def overlaps(read, interval_start, interval_end):
        return read.reference_start < interval_end and read.reference_end > interval_start

    anchors = [
        read
        for read in reads
        if overlaps(read, exon14_start, exon14_end)
        and overlaps(read, exon15_start, exon15_end)
    ]
    if not anchors:
        anchors = reads

    pair_counts = Counter((r.reference_start, r.reference_end) for r in anchors)
    clusters = _cluster_endpoint_pairs(pair_counts, tolerance)
    if not clusters:
        raise ValueError("Unable to infer a FLT3 amplicon from BAM alignment endpoints")

    max_support = max(int(cluster["support"]) for cluster in clusters)
    min_support = max(min_cluster_support, math.ceil(max_support * min_cluster_fraction))

    amplicons = []
    for cluster in clusters:
        if int(cluster["support"]) < min_support:
            continue
        start, end = cluster["representative"]
        if start >= end:
            continue
        amplicons.append(Amplicon(start=start, end=end, support=int(cluster["support"])))

    amplicons.sort(key=lambda amp: (-amp.support, amp.start, amp.end))
    if not amplicons:
        raise ValueError("No sufficiently supported FLT3 amplicon endpoint cluster found")
    return amplicons


def assign_reads_to_amplicons(
    reads,
    amplicons: list[Amplicon],
    end_anchor_tolerance: int = DEFAULT_END_ANCHOR_TOLERANCE,
    min_overlap_fraction: float = DEFAULT_MIN_OVERLAP_FRACTION,
):
    assigned: list[list] = [[] for _ in amplicons]
    for read in reads:
        scores = [
            abs(read.reference_start - amp.start) + abs(read.reference_end - amp.end)
            for amp in amplicons
        ]
        idx = min(range(len(scores)), key=scores.__getitem__)
        amp = amplicons[idx]
        overlap = max(0, min(read.reference_end, amp.end) - max(read.reference_start, amp.start))
        if overlap == 0:
            continue
        # At least one amplicon end should remain anchored, or the genomic span
        # should still substantially overlap the inferred WT amplicon.
        end_anchor = min(
            abs(read.reference_start - amp.start),
            abs(read.reference_end - amp.end),
        ) <= end_anchor_tolerance
        if end_anchor or overlap >= min_overlap_fraction * amp.length:
            assigned[idx].append(read)
    return assigned


def sequence_counts(reads, strand: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    for read in reads:
        # BAM SEQ is stored in reference alignment orientation. FLT3.fa is in
        # transcript orientation (minus strand for this hg19 reference).
        sequence = read.query_sequence
        if strand == "-":
            sequence = reverse_complement(sequence)
        counts[sequence.upper()] += 1
    return counts


def write_reference(path: Path, sequence: str) -> None:
    path.write_text(sequence.upper() + "\n")


def write_annotation(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = ["amplicon_bp", "region", "chr13_bp", "transcript_bp", "protein_as"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_sequence_counts(path: Path, counts: Counter[str], min_read_copies: int) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sequence", "Counts", "SeqLength"])
        for sequence, count in counts.most_common():
            if count >= min_read_copies:
                writer.writerow([sequence, count, len(sequence)])


def run_mergeitd(
    runner: Path,
    workdir: Path,
    sample: str,
    reads_tsv: Path,
    reference_txt: Path,
    annotation_tsv: Path,
    min_insert_seq_length: int,
    min_total_reads: int,
    min_vaf: float,
    match_score: float,
    mismatch_score: float,
    gap_open_score: float,
    gap_extend_score: float,
    min_aligned_block: int,
    min_reference_fraction: float,
    max_indel_fraction: float,
) -> Path:
    output_dir = workdir / "mergeitd"
    cmd = [
        sys.executable,
        str(runner),
        "--sample", sample,
        "--reads", str(reads_tsv),
        "--reference", str(reference_txt),
        "--annotation", str(annotation_tsv),
        "--output-dir", str(output_dir),
        "--min-insert-seq-length", str(min_insert_seq_length),
        "--min-total-reads", str(min_total_reads),
        "--min-vaf", str(min_vaf),
        "--match-score", str(match_score),
        "--mismatch-score", str(mismatch_score),
        "--gap-open-score", str(gap_open_score),
        "--gap-extend-score", str(gap_extend_score),
        "--min-aligned-block", str(min_aligned_block),
        "--min-reference-fraction", str(min_reference_fraction),
        "--max-indel-fraction", str(max_indel_fraction),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        details = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"mergeITD runner failed: {details}")
    return output_dir / "filtered_mut_vaf.csv"


def read_calls(path: Path, amplicon: Amplicon, contig: str) -> list[dict[str, str]]:
    calls = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            hgvs = row.get("HGVS", "")
            if not hgvs:
                continue
            if row.get("insertRegion") not in {"exon14", "intron14", "exon15"}:
                continue
            row["HGVS"] = f"{TRANSCRIPT}:{hgvs}"
            row["amplicon"] = f"{contig}:{amplicon.label}"
            calls.append(row)
    return calls


def _add_bool_override(
    parser: argparse.ArgumentParser,
    name: str,
    dest: str,
    enable_help: str,
    disable_help: str,
) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        f"--{name}", dest=dest, action="store_true", default=None, help=enable_help
    )
    group.add_argument(
        f"--no-{name}", dest=dest, action="store_false", help=disable_help
    )


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Call FLT3 insertions from an indexed hg19 BBmerged BAM using mergeITD alignment semantics."
    )
    parser.add_argument("bam", type=Path, help="coordinate-sorted, indexed single-read BAM")
    parser.add_argument(
        "--settings",
        type=Path,
        default=None,
        help="partial TOML settings file; CLI options override file values",
    )
    parser.add_argument(
        "--flt3-reference",
        type=Path,
        default=script_dir / "annotation" / "FLT3.fa",
        help="hg19 FLT3 genomic FASTA in transcript orientation (default: annotation/FLT3.fa)",
    )

    parser.add_argument("--min-mapq", type=int, default=None)
    parser.add_argument("--min-mean-base-quality", type=float, default=None)
    _add_bool_override(
        parser,
        "exclude-duplicates",
        "exclude_duplicates",
        "Exclude BAM duplicate-flagged reads",
        "Include BAM duplicate-flagged reads",
    )
    _add_bool_override(
        parser,
        "exclude-qcfail",
        "exclude_qcfail",
        "Exclude BAM QC-fail-flagged reads",
        "Include BAM QC-fail-flagged reads",
    )
    parser.add_argument("--min-read-copies", type=int, default=None)

    parser.add_argument("--endpoint-tolerance", type=int, default=None)
    parser.add_argument("--min-cluster-fraction", type=float, default=None)
    parser.add_argument("--min-cluster-support", type=int, default=None)
    parser.add_argument("--target-fetch-flank", type=int, default=None)
    parser.add_argument("--end-anchor-tolerance", type=int, default=None)
    parser.add_argument("--min-overlap-fraction", type=float, default=None)

    parser.add_argument("--match-score", type=float, default=None)
    parser.add_argument("--mismatch-score", type=float, default=None)
    parser.add_argument("--gap-open-score", type=float, default=None)
    parser.add_argument("--gap-extend-score", type=float, default=None)
    parser.add_argument("--min-aligned-block", type=int, default=None)
    parser.add_argument("--min-reference-fraction", type=float, default=None)
    parser.add_argument("--max-indel-fraction", type=float, default=None)

    parser.add_argument("--min-insert-seq-length", type=int, default=None)
    parser.add_argument("--min-total-reads", type=int, default=None)
    parser.add_argument("--min-vaf", type=float, default=None, help="VAF percent, matching mergeITD")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def resolve_args(args: argparse.Namespace) -> argparse.Namespace:
    settings = load_settings(args.settings)
    mappings = {
        "min_mapq": ("reads", "min_mapq"),
        "min_mean_base_quality": ("reads", "min_mean_base_quality"),
        "exclude_duplicates": ("reads", "exclude_duplicates"),
        "exclude_qcfail": ("reads", "exclude_qcfail"),
        "min_read_copies": ("reads", "min_read_copies"),
        "endpoint_tolerance": ("amplicon", "endpoint_tolerance"),
        "min_cluster_fraction": ("amplicon", "min_cluster_fraction"),
        "min_cluster_support": ("amplicon", "min_cluster_support"),
        "target_fetch_flank": ("amplicon", "target_fetch_flank"),
        "end_anchor_tolerance": ("amplicon", "end_anchor_tolerance"),
        "min_overlap_fraction": ("amplicon", "min_overlap_fraction"),
        "match_score": ("alignment", "match"),
        "mismatch_score": ("alignment", "mismatch"),
        "gap_open_score": ("alignment", "gap_open"),
        "gap_extend_score": ("alignment", "gap_extend"),
        "min_aligned_block": ("alignment", "min_aligned_block"),
        "min_reference_fraction": ("alignment", "min_reference_fraction"),
        "max_indel_fraction": ("alignment", "max_indel_fraction"),
        "min_insert_seq_length": ("calling", "min_insert_length"),
        "min_total_reads": ("calling", "min_supporting_reads"),
        "min_vaf": ("calling", "min_vaf_percent"),
    }
    effective = deepcopy(settings)
    for attr, (section, key) in mappings.items():
        if getattr(args, attr) is None:
            setattr(args, attr, settings[section][key])
        effective[section][key] = getattr(args, attr)
    validate_settings(effective, source="CLI/settings")
    return args


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    runner = script_dir / "mergeitd_runner.py"

    try:
        args = resolve_args(args)
        reference = Flt3Reference.load(args.flt3_reference)
        bam = validate_bam(args.bam)
        try:
            contig = resolve_bam_contig(bam, reference.chrom)
            target_start0, target_end0 = reference.target_fetch_interval()
            exon14_interval = reference.genomic_interval_for_region("exon14")
            exon15_interval = reference.genomic_interval_for_region("exon15")
            fetch_start0 = max(0, target_start0 - args.target_fetch_flank)
            fetch_end0 = target_end0 + args.target_fetch_flank
            reads = primary_target_reads(
                bam,
                contig,
                fetch_start0,
                fetch_end0,
                min_mapq=args.min_mapq,
                min_mean_base_quality=args.min_mean_base_quality,
                exclude_duplicates=args.exclude_duplicates,
                exclude_qcfail=args.exclude_qcfail,
            )
            amplicons = infer_amplicons(
                reads,
                exon14_interval,
                exon15_interval,
                tolerance=args.endpoint_tolerance,
                min_cluster_fraction=args.min_cluster_fraction,
                min_cluster_support=args.min_cluster_support,
            )
            read_groups = assign_reads_to_amplicons(
                reads,
                amplicons,
                end_anchor_tolerance=args.end_anchor_tolerance,
                min_overlap_fraction=args.min_overlap_fraction,
            )
        finally:
            bam.close()

        all_calls: list[dict[str, str]] = []
        with tempfile.TemporaryDirectory(prefix="itdmapper_") as tmp:
            temp_root = Path(tmp)
            for idx, (amplicon, amp_reads) in enumerate(zip(amplicons, read_groups), start=1):
                if not amp_reads:
                    continue
                workdir = temp_root / f"amplicon_{idx}"
                workdir.mkdir()

                ref_seq = reference.extract_interval(amplicon.start, amplicon.end)
                anno_rows = reference.annotation_rows(amplicon.start, amplicon.end)
                counts = sequence_counts(amp_reads, reference.strand)

                reference_txt = workdir / "reference.txt"
                annotation_tsv = workdir / "annotation.tsv"
                reads_tsv = workdir / "reads.tsv"
                write_reference(reference_txt, ref_seq)
                write_annotation(annotation_tsv, anno_rows)
                write_sequence_counts(reads_tsv, counts, args.min_read_copies)

                if args.verbose:
                    print(
                        f"Inferred {contig}:{amplicon.start + 1}-{amplicon.end} "
                        f"({amplicon.support} anchor reads; {len(amp_reads)} reads assigned)",
                        file=sys.stderr,
                    )

                filtered = run_mergeitd(
                    runner=runner,
                    workdir=workdir,
                    sample=args.bam.stem,
                    reads_tsv=reads_tsv,
                    reference_txt=reference_txt,
                    annotation_tsv=annotation_tsv,
                    min_insert_seq_length=args.min_insert_seq_length,
                    min_total_reads=args.min_total_reads,
                    min_vaf=args.min_vaf,
                    match_score=args.match_score,
                    mismatch_score=args.mismatch_score,
                    gap_open_score=args.gap_open_score,
                    gap_extend_score=args.gap_extend_score,
                    min_aligned_block=args.min_aligned_block,
                    min_reference_fraction=args.min_reference_fraction,
                    max_indel_fraction=args.max_indel_fraction,
                )
                all_calls.extend(read_calls(filtered, amplicon, contig))

        if not all_calls:
            print("NO_FLT3_ITD")
            return 0

        # Multiple inferred amplicons may independently observe the same call.
        # Report the strongest observation for each HGVS string.
        best: dict[str, dict[str, str]] = {}
        for call in all_calls:
            key = call["HGVS"]
            if key not in best or float(call["vaf_percent"]) > float(best[key]["vaf_percent"]):
                best[key] = call

        print("HGVS\tnetInsert\tcounts\tvaf_percent\tcoverage\tamplicon")
        for call in sorted(best.values(), key=lambda row: float(row["vaf_percent"]), reverse=True):
            print(
                "\t".join(
                    [
                        call["HGVS"],
                        call["netInsert"],
                        call["counts"],
                        call["vaf_percent"],
                        call["coverage"],
                        call["amplicon"],
                    ]
                )
            )
        return 0
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
