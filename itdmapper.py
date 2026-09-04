#!/usr/bin/env python3
"""FLT3 ITD caller for BBmerged BAM or FASTQ inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
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

from Bio import Align, SeqIO

from hgvs_normalize import normalize_mergeitd_name
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


@dataclass(frozen=True)
class AlignmentMetrics:
    sequence: str
    score: float
    query_start: int
    query_end: int
    reference_start: int
    reference_end: int
    query_aligned_bases: int
    indel_bases: int

    @property
    def query_span(self) -> int:
        return self.query_end - self.query_start

    @property
    def reference_span(self) -> int:
        return self.reference_end - self.reference_start


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

    def hgvs_coord_for_index(self, idx: int) -> str:
        """Return c.-coordinate text for a full-reference sequence index.

        Intronic positions are numbered from the nearest flanking coding exon,
        using the downstream exon on an exact midpoint tie to match mergeITD's
        existing annotation convention.
        """
        if idx < 0 or idx >= len(self.sequence):
            raise ValueError(f"Reference index {idx} lies outside FLT3.fa")
        c_pos = self.c_by_index[idx]
        if c_pos is not None:
            return str(c_pos)

        left = idx - 1
        while left >= 0 and self.c_by_index[left] is None:
            left -= 1
        right = idx + 1
        while right < len(self.sequence) and self.c_by_index[right] is None:
            right += 1
        if left < 0 or right >= len(self.sequence):
            raise ValueError(
                f"Reference index {idx} cannot be represented as an intronic c. coordinate"
            )

        left_c = self.c_by_index[left]
        right_c = self.c_by_index[right]
        assert left_c is not None and right_c is not None
        left_distance = idx - left
        right_distance = right - idx
        if left_distance < right_distance:
            return f"{left_c}+{left_distance}"
        return f"{right_c}-{right_distance}"

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

    def reference_span_to_genomic(
        self,
        interval_start0: int,
        interval_end0: int,
        reference_start: int,
        reference_end: int,
    ) -> tuple[int, int]:
        """Convert coordinates within extract_interval() back to genomic coordinates."""
        if self.strand == "-":
            return interval_end0 - reference_end, interval_end0 - reference_start
        return interval_start0 + reference_start, interval_start0 + reference_end

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


def detect_input_format(path: Path, override: Optional[str] = None) -> str:
    if override is not None:
        return override
    name = path.name.lower()
    if name.endswith(".bam"):
        return "bam"
    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"
    raise ValueError(
        f"Unable to infer input format from {path.name!r}; use --input-format bam|fastq"
    )


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


def load_fastq_sequence_counts(
    fastq_path: Path,
    min_mean_base_quality: float,
    min_read_copies: int,
) -> tuple[Counter[str], int, int]:
    """Quality-filter raw FASTQ reads, then collapse exact sequences.

    Returns (filtered sequence counts, raw read count, quality-passing read count).
    The minimum-copy filter is applied after collapsing and before discovery or
    amplicon assignment.
    """
    if not fastq_path.is_file():
        raise ValueError(f"FASTQ not found: {fastq_path}")
    opener = gzip.open if fastq_path.name.lower().endswith(".gz") else open
    counts: Counter[str] = Counter()
    raw_reads = 0
    quality_passing_reads = 0
    try:
        with opener(fastq_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                raw_reads += 1
                qualities = record.letter_annotations.get("phred_quality", [])
                if min_mean_base_quality > 0:
                    if not qualities or (sum(qualities) / len(qualities)) < min_mean_base_quality:
                        continue
                sequence = str(record.seq).upper()
                if not sequence:
                    continue
                quality_passing_reads += 1
                counts[sequence] += 1
    except (ValueError, OSError) as exc:
        raise ValueError(f"Unable to read FASTQ {fastq_path}: {exc}") from exc

    if raw_reads == 0:
        raise ValueError("FASTQ contains no reads")
    if quality_passing_reads == 0:
        raise ValueError("No FASTQ reads pass configured read-quality filters")

    counts = Counter({seq: n for seq, n in counts.items() if n >= min_read_copies})
    if not counts:
        raise ValueError(
            "No FASTQ sequences remain after applying the minimum identical-read-copy filter"
        )
    return counts, raw_reads, quality_passing_reads


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


def _clusters_to_amplicons(
    clusters: list[dict[str, object]],
    min_cluster_fraction: float,
    min_cluster_support: int,
    error_message: str,
) -> list[Amplicon]:
    if not clusters:
        raise ValueError(error_message)
    max_support = max(int(cluster["support"]) for cluster in clusters)
    min_support = max(min_cluster_support, math.ceil(max_support * min_cluster_fraction))
    amplicons = []
    for cluster in clusters:
        support = int(cluster["support"])
        if support < min_support:
            continue
        start, end = cluster["representative"]
        if start < end:
            amplicons.append(Amplicon(start=start, end=end, support=support))
    amplicons.sort(key=lambda amp: (-amp.support, amp.start, amp.end))
    if not amplicons:
        raise ValueError(error_message)
    return amplicons


def infer_amplicons(
    reads,
    exon14_interval: tuple[int, int],
    exon15_interval: tuple[int, int],
    tolerance: int = DEFAULT_CLUSTER_TOLERANCE,
    min_cluster_fraction: float = DEFAULT_MIN_CLUSTER_FRACTION,
    min_cluster_support: int = DEFAULT_MIN_CLUSTER_SUPPORT,
) -> list[Amplicon]:
    """Infer BAM amplicons from dominant genomic alignment endpoint pairs."""
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
    return _clusters_to_amplicons(
        clusters,
        min_cluster_fraction,
        min_cluster_support,
        "No sufficiently supported FLT3 amplicon endpoint cluster found",
    )


def _make_aligner(match: float, mismatch: float, gap_open: float, gap_extend: float):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    return aligner


def _alignment_metrics(
    sequence: str,
    reference: str,
    aligner,
    min_aligned_block: int,
) -> Optional[AlignmentMetrics]:
    alignments = aligner.align(sequence, reference)
    if not alignments:
        return None
    aln = alignments[-1]
    coords = aln.coordinates
    q_blocks: list[tuple[int, int]] = []
    r_blocks: list[tuple[int, int]] = []
    for j in range(len(coords[0]) - 1):
        q0, q1 = int(coords[0][j]), int(coords[0][j + 1])
        r0, r1 = int(coords[1][j]), int(coords[1][j + 1])
        if q1 - q0 >= min_aligned_block and r1 - r0 >= min_aligned_block:
            q_blocks.append((q0, q1))
            r_blocks.append((r0, r1))
    if not q_blocks:
        return None

    indel_bases = 0
    for i in range(len(q_blocks) - 1):
        q_gap = q_blocks[i + 1][0] - q_blocks[i][1]
        r_gap = r_blocks[i + 1][0] - r_blocks[i][1]
        indel_bases += max(0, q_gap) + max(0, r_gap)

    return AlignmentMetrics(
        sequence=sequence,
        score=float(aln.score),
        query_start=q_blocks[0][0],
        query_end=q_blocks[-1][1],
        reference_start=r_blocks[0][0],
        reference_end=r_blocks[-1][1],
        query_aligned_bases=sum(q1 - q0 for q0, q1 in q_blocks),
        indel_bases=indel_bases,
    )


def _best_orientation_metrics(
    sequence: str,
    reference: str,
    aligner,
    min_aligned_block: int,
) -> Optional[AlignmentMetrics]:
    options = []
    for oriented in (sequence.upper(), reverse_complement(sequence.upper())):
        metrics = _alignment_metrics(oriented, reference, aligner, min_aligned_block)
        if metrics is not None:
            options.append(metrics)
    if not options:
        return None
    return max(
        options,
        key=lambda m: (m.score, m.query_aligned_bases, m.reference_span),
    )


def discover_amplicons_from_sequences(
    sequence_counts: Counter[str],
    reference: Flt3Reference,
    target_interval: tuple[int, int],
    discovery_flank: int,
    top_unique: int,
    endpoint_tolerance: int,
    min_cluster_fraction: float,
    min_cluster_support: int,
    match_score: float,
    mismatch_score: float,
    gap_open_score: float,
    gap_extend_score: float,
    min_aligned_block: int,
    min_query_fraction: float,
    max_indel_fraction: float,
) -> list[Amplicon]:
    """Infer FASTQ amplicons from the most abundant unique sequences.

    Discovery is performed against a broad FLT3 window around the exon14-15 target.
    Candidate amplicons must overlap the target but may extend beyond it.
    Endpoint support is weighted by the collapsed read count.
    """
    target_start0, target_end0 = target_interval
    reference_start0 = reference.genomic_start - 1
    reference_end0 = reference.genomic_end
    discovery_start0 = max(reference_start0, target_start0 - discovery_flank)
    discovery_end0 = min(reference_end0, target_end0 + discovery_flank)
    discovery_ref = reference.extract_interval(discovery_start0, discovery_end0).upper()
    aligner = _make_aligner(match_score, mismatch_score, gap_open_score, gap_extend_score)

    pair_counts: Counter[tuple[int, int]] = Counter()
    for sequence, count in sequence_counts.most_common(top_unique):
        metrics = _best_orientation_metrics(
            sequence, discovery_ref, aligner, min_aligned_block
        )
        if metrics is None:
            continue
        query_fraction = metrics.query_aligned_bases / len(metrics.sequence)
        denominator = metrics.indel_bases + metrics.reference_span
        indel_fraction = metrics.indel_bases / denominator if denominator else 1.0
        if query_fraction < min_query_fraction or indel_fraction > max_indel_fraction:
            continue
        start0, end0 = reference.reference_span_to_genomic(
            discovery_start0,
            discovery_end0,
            metrics.reference_start,
            metrics.reference_end,
        )
        if start0 >= target_end0 or end0 <= target_start0:
            continue
        pair_counts[(start0, end0)] += count

    clusters = _cluster_endpoint_pairs(pair_counts, endpoint_tolerance)
    return _clusters_to_amplicons(
        clusters,
        min_cluster_fraction,
        min_cluster_support,
        "Unable to infer a sufficiently supported FLT3 amplicon from FASTQ sequences",
    )


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
        end_anchor = min(
            abs(read.reference_start - amp.start),
            abs(read.reference_end - amp.end),
        ) <= end_anchor_tolerance
        if end_anchor or overlap >= min_overlap_fraction * amp.length:
            assigned[idx].append(read)
    return assigned


def assign_sequence_counts_to_amplicons(
    sequence_counts: Counter[str],
    amplicons: list[Amplicon],
    reference: Flt3Reference,
    match_score: float,
    mismatch_score: float,
    gap_open_score: float,
    gap_extend_score: float,
    min_aligned_block: int,
    min_reference_fraction: float,
    min_query_fraction: float,
    max_indel_fraction: float,
) -> tuple[list[Counter[str]], int, int]:
    """Assign collapsed FASTQ sequences to their best acceptable amplicon.

    Both read orientations are tested. A sequence must satisfy the configured
    reference coverage, query coverage and indel-fraction criteria against at
    least one amplicon. Accepted sequences are normalized to FLT3 transcript
    orientation before downstream mergeITD calling.
    """
    aligner = _make_aligner(match_score, mismatch_score, gap_open_score, gap_extend_score)
    amp_refs = [reference.extract_interval(amp.start, amp.end).upper() for amp in amplicons]
    assigned: list[Counter[str]] = [Counter() for _ in amplicons]
    accepted_reads = 0
    rejected_reads = 0

    for sequence, count in sequence_counts.items():
        candidates = []
        for idx, (amp, amp_ref) in enumerate(zip(amplicons, amp_refs)):
            metrics = _best_orientation_metrics(sequence, amp_ref, aligner, min_aligned_block)
            if metrics is None:
                continue
            reference_fraction = metrics.reference_span / amp.length
            query_fraction = metrics.query_aligned_bases / len(metrics.sequence)
            denominator = metrics.indel_bases + metrics.reference_span
            indel_fraction = metrics.indel_bases / denominator if denominator else 1.0
            if reference_fraction < min_reference_fraction:
                continue
            if query_fraction < min_query_fraction:
                continue
            if indel_fraction > max_indel_fraction:
                continue
            candidates.append(
                (
                    reference_fraction,
                    query_fraction,
                    metrics.score,
                    amp.support,
                    -idx,
                    idx,
                    metrics.sequence,
                )
            )

        if not candidates:
            rejected_reads += count
            continue
        best = max(candidates)
        idx = best[5]
        oriented_sequence = best[6]
        assigned[idx][oriented_sequence] += count
        accepted_reads += count

    return assigned, accepted_reads, rejected_reads


def sequence_counts(reads, strand: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    for read in reads:
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


def write_sequence_counts(path: Path, counts: Counter[str], min_read_copies: int = 1) -> None:
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


def read_calls(
    path: Path,
    amplicon: Amplicon,
    contig: str,
    reference: Flt3Reference,
) -> list[dict[str, str]]:
    calls = []
    amplicon_reference = reference.extract_interval(amplicon.start, amplicon.end).upper()
    first_genomic_pos = amplicon.end if reference.strand == "-" else amplicon.start + 1
    full_reference_offset = reference.sequence_index_for_genomic_pos(first_genomic_pos)

    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if not row.get("HGVS", ""):
                continue
            mutation_name = row.get("name", "")
            if not mutation_name:
                raise ValueError("mergeITD call is missing its raw mutation name")

            normalized = normalize_mergeitd_name(
                mutation_name=mutation_name,
                amplicon_reference=amplicon_reference,
                full_reference=reference.sequence,
                full_reference_offset=full_reference_offset,
                coord_for_index=reference.hgvs_coord_for_index,
            )
            row["HGVS"] = f"{TRANSCRIPT}:{normalized.hgvs}"

            site_idx = normalized.site_index
            if site_idx == len(reference.sequence):
                site_idx -= 1
            if 0 <= site_idx < len(reference.region_by_index):
                row["insertRegion"] = reference.region_by_index[site_idx]
            if row.get("insertRegion") not in {"exon14", "intron14", "exon15"}:
                continue

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
        description="Call FLT3 insertions from already-BBmerged BAM or FASTQ fragments."
    )
    parser.add_argument("input", type=Path, help="input .bam, .fastq/.fq, or gzipped FASTQ")
    parser.add_argument(
        "--input-format",
        choices=("bam", "fastq"),
        default=None,
        help="override automatic input-format detection",
    )
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
    parser.add_argument("--fastq-discovery-top-unique", type=int, default=None)

    parser.add_argument("--match-score", type=float, default=None)
    parser.add_argument("--mismatch-score", type=float, default=None)
    parser.add_argument("--gap-open-score", type=float, default=None)
    parser.add_argument("--gap-extend-score", type=float, default=None)
    parser.add_argument("--min-aligned-block", type=int, default=None)
    parser.add_argument("--min-reference-fraction", type=float, default=None)
    parser.add_argument("--min-query-fraction", type=float, default=None)
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
        "fastq_discovery_top_unique": ("fastq", "discovery_top_unique"),
        "match_score": ("alignment", "match"),
        "mismatch_score": ("alignment", "mismatch"),
        "gap_open_score": ("alignment", "gap_open"),
        "gap_extend_score": ("alignment", "gap_extend"),
        "min_aligned_block": ("alignment", "min_aligned_block"),
        "min_reference_fraction": ("alignment", "min_reference_fraction"),
        "min_query_fraction": ("alignment", "min_query_fraction"),
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


def _sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    runner = script_dir / "mergeitd_runner.py"

    try:
        args = resolve_args(args)
        input_format = detect_input_format(args.input, args.input_format)
        reference = Flt3Reference.load(args.flt3_reference)
        target_start0, target_end0 = reference.target_fetch_interval()
        contig = reference.chrom

        if input_format == "bam":
            bam = validate_bam(args.input)
            try:
                contig = resolve_bam_contig(bam, reference.chrom)
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
                bam_groups = assign_reads_to_amplicons(
                    reads,
                    amplicons,
                    end_anchor_tolerance=args.end_anchor_tolerance,
                    min_overlap_fraction=args.min_overlap_fraction,
                )
                read_groups = [sequence_counts(group, reference.strand) for group in bam_groups]
                min_output_copies = args.min_read_copies
            finally:
                bam.close()
        else:
            sequence_count_map, raw_reads, quality_passing_reads = load_fastq_sequence_counts(
                args.input,
                min_mean_base_quality=args.min_mean_base_quality,
                min_read_copies=args.min_read_copies,
            )
            amplicons = discover_amplicons_from_sequences(
                sequence_count_map,
                reference,
                (target_start0, target_end0),
                discovery_flank=args.target_fetch_flank,
                top_unique=args.fastq_discovery_top_unique,
                endpoint_tolerance=args.endpoint_tolerance,
                min_cluster_fraction=args.min_cluster_fraction,
                min_cluster_support=args.min_cluster_support,
                match_score=args.match_score,
                mismatch_score=args.mismatch_score,
                gap_open_score=args.gap_open_score,
                gap_extend_score=args.gap_extend_score,
                min_aligned_block=args.min_aligned_block,
                min_query_fraction=args.min_query_fraction,
                max_indel_fraction=args.max_indel_fraction,
            )
            read_groups, accepted_reads, rejected_reads = assign_sequence_counts_to_amplicons(
                sequence_count_map,
                amplicons,
                reference,
                match_score=args.match_score,
                mismatch_score=args.mismatch_score,
                gap_open_score=args.gap_open_score,
                gap_extend_score=args.gap_extend_score,
                min_aligned_block=args.min_aligned_block,
                min_reference_fraction=args.min_reference_fraction,
                min_query_fraction=args.min_query_fraction,
                max_indel_fraction=args.max_indel_fraction,
            )
            min_output_copies = 1  # already filtered before FASTQ discovery/assignment
            if accepted_reads == 0:
                raise ValueError("No FASTQ sequences satisfy alignment criteria for any inferred amplicon")
            if args.verbose:
                retained_after_copy_filter = sum(sequence_count_map.values())
                print(
                    f"FASTQ reads: {raw_reads}; quality-passing: {quality_passing_reads}; "
                    f"retained after copy filter: {retained_after_copy_filter}; "
                    f"assigned: {accepted_reads}; rejected by amplicon alignment: {rejected_reads}",
                    file=sys.stderr,
                )

        all_calls: list[dict[str, str]] = []
        with tempfile.TemporaryDirectory(prefix="itdmapper_") as tmp:
            temp_root = Path(tmp)
            for idx, (amplicon, counts) in enumerate(zip(amplicons, read_groups), start=1):
                if not counts:
                    continue
                workdir = temp_root / f"amplicon_{idx}"
                workdir.mkdir()

                ref_seq = reference.extract_interval(amplicon.start, amplicon.end)
                anno_rows = reference.annotation_rows(amplicon.start, amplicon.end)

                reference_txt = workdir / "reference.txt"
                annotation_tsv = workdir / "annotation.tsv"
                reads_tsv = workdir / "reads.tsv"
                write_reference(reference_txt, ref_seq)
                write_annotation(annotation_tsv, anno_rows)
                write_sequence_counts(reads_tsv, counts, min_output_copies)

                if args.verbose:
                    print(
                        f"Inferred {contig}:{amplicon.start + 1}-{amplicon.end} "
                        f"({amplicon.support} discovery-support reads; "
                        f"{sum(counts.values())} reads assigned)",
                        file=sys.stderr,
                    )

                filtered = run_mergeitd(
                    runner=runner,
                    workdir=workdir,
                    sample=_sample_name(args.input),
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
                all_calls.extend(read_calls(filtered, amplicon, contig, reference))

        if not all_calls:
            print("NO_FLT3_ITD")
            return 0

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
