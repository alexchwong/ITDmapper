import csv
import random
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from allele_reassignment import build_references, quantify, reassign_sequences, recount_sequences
from hgvs_normalize import ReferenceEdit, apply_edit
from itdmapper import _make_aligner, parse_args, resolve_args


class TestAlleleReassignment(unittest.TestCase):
    def setUp(self):
        rng = random.Random(42)
        self.wt = "".join(rng.choice("ACGT") for _ in range(240))
        self.edit = ReferenceEdit(120, 120, self.wt[90:120])
        self.candidates = {"ITD": self.edit}
        self.mutant = apply_edit(self.wt, self.edit)
        self.references = build_references([("amp", self.wt, 0)], self.candidates)
        self.aligner = _make_aligner(5, -15, -36, -0.5)

    def run_assignment(self, counts, **kwargs):
        return reassign_sequences(counts, self.references, self.candidates, self.aligner, **kwargs)

    def test_counts_rescue_short_flank_and_ambiguity(self):
        counts = Counter({(self.wt, "+"): 90, (self.mutant, "-"): 8,
                          (self.mutant[:155], "+"): 2, (self.wt[:80], "+"): 7})
        rows, assigned = self.run_assignment(counts)
        support, depth = quantify(self.candidates, assigned, 3)
        self.assertEqual(support["ITD"], 10)
        self.assertEqual(depth["ITD"], 100)
        self.assertEqual(sum(row["count"] for row in rows), 107)
        self.assertEqual(next(row for row in rows if row["count"] == 7)["status"], "ambiguous")
        self.assertEqual(sum(count for match, strand, count in assigned), 100)

    def test_both_orientations(self):
        from Bio.Seq import Seq
        rows, assigned = self.run_assignment(
            Counter({(str(Seq(self.mutant).reverse_complement()), "?"): 5}), both_orientations=True)
        self.assertEqual(rows[0]["status"], "assigned")
        self.assertEqual(assigned[0][0].sequence, self.mutant)

    def test_unexplained_terminal_sequence_rejected(self):
        rows, assigned = self.run_assignment(Counter({(self.mutant + "N" * 150, "+"): 4}))
        self.assertEqual(rows[0]["reason"], "insufficient_query_coverage")
        self.assertFalse(assigned)

    def test_duplicate_amplicons_do_not_create_allele_ties(self):
        self.references = build_references([("amp", self.wt, 0), ("other", self.wt, 0)], self.candidates)
        rows, assigned = self.run_assignment(Counter({(self.mutant, "+"): 5}))
        self.assertEqual(rows[0]["status"], "assigned")
        self.assertEqual(sum(count for match, strand, count in assigned), 5)

    def test_competing_itds_and_candidate_specific_coverage(self):
        edit2 = ReferenceEdit(170, 170, self.wt[150:170])
        self.candidates["other"] = edit2
        self.references = build_references([("amp", self.wt, 0)], self.candidates)
        rows, assigned = self.run_assignment(Counter({(self.mutant[:155], "+"): 2,
            (apply_edit(self.wt, edit2), "+"): 3, (self.wt, "+"): 10}))
        support, coverage = quantify(self.candidates, assigned, 3)
        self.assertEqual(support, Counter({"ITD": 2, "other": 3}))
        self.assertEqual(coverage, Counter({"ITD": 15, "other": 13}))

    def test_margin_and_anchor_requirements(self):
        rows, assigned = self.run_assignment(Counter({(self.mutant, "+"): 2}), min_score_margin=1000)
        self.assertEqual(rows[0]["reason"], "ambiguous_score")
        rows, assigned = self.run_assignment(Counter({(self.mutant[:152], "+"): 2}))
        self.assertEqual(rows[0]["reason"], "insufficient_variant_anchors")

    def test_settings_opt_in_and_validation(self):
        self.assertFalse(resolve_args(parse_args(["sample.bam"])).recalculate_vaf)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "settings.toml"
            path.write_text("[reassignment]\nenabled = true\nmin_score_margin = 12\n")
            args = resolve_args(parse_args(["sample.bam", "--settings", str(path)]))
            self.assertTrue(args.recalculate_vaf)
            self.assertEqual(args.reassignment_min_score_margin, 12)
            self.assertFalse(resolve_args(parse_args(["sample.bam", "--settings", str(path),
                                                      "--no-recalculate-vaf"])).recalculate_vaf)
        for options in [["--reassignment-min-score-margin", "-1"],
                        ["--reassignment-min-query-fraction", "1.1"],
                        ["--reassignment-min-anchor", "0"], ["--reassignment-audit", "audit.tsv"],
                        ["--recalculate-vaf", "--reassignment-audit", "sample.bam"]]:
            with self.assertRaises(ValueError):
                resolve_args(parse_args(["sample.bam"] + options))

    def test_wt_only_and_error_filter(self):
        refs = build_references([("amp", self.wt, 0)], {})
        rows, assigned = reassign_sequences(Counter({(self.wt, "+"): 12}), refs, {}, self.aligner)
        self.assertEqual(rows[0]["allele"], "WT")
        self.assertEqual(rows[0]["status"], "assigned")
        damaged = list(self.mutant)
        damaged[20] = "A" if damaged[20] != "A" else "C"
        rows, assigned = self.run_assignment(Counter({("".join(damaged), "+"): 1}), max_error_fraction=0)
        self.assertEqual(rows[0]["reason"], "excess_alignment_errors")

    def test_bam_exports_agree_with_recalculated_counts(self):
        import pysam
        from itdmapper import Flt3Reference, reverse_complement
        from Bio import SeqIO
        ref = Flt3Reference.load(ROOT / "annotation" / "FLT3.fa")
        start, end = ref.target_fetch_interval()
        wt = ref.extract_interval(start, end).upper()
        mutant = wt[:70] + wt[52:70] + wt[70:]
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "reads.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": ref.chrom, "LN": ref.genomic_end + 1000}]}
            with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
                groups = [(wt, [(0, len(wt))], 20),
                          (mutant, [(0, len(wt) - 70), (1, 18), (0, 70)], 6)]
                for seq, cigar, count in groups:
                    for i in range(count):
                        read = pysam.AlignedSegment()
                        read.query_name = f"read{handle.tell()}_{i}"
                        read.query_sequence = reverse_complement(seq)
                        read.reference_id = 0
                        read.reference_start = start
                        read.mapping_quality = 60
                        read.flag = 16 if i % 2 else 0
                        read.cigartuples = cigar
                        read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                        handle.write(read)
            pysam.index(str(bam))
            command = [sys.executable, str(ROOT / "itdmapper.py"), str(bam),
                       "--recalculate-vaf", "--variant-fasta", "--variant-alignments"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            calls = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
            self.assertEqual(len(calls), 1, result.stdout)
            self.assertEqual(int(calls[0]["counts"]), 6)
            self.assertEqual(int(calls[0]["coverage"]), 26)
            records = list(SeqIO.parse(Path(tmp) / "reads_variant001.fasta", "fasta"))
            self.assertEqual(len(records), 2)
            self.assertTrue(all("count3_" in record.id for record in records))
            text = (Path(tmp) / "reads_variant001_alignments.txt").read_text()
            self.assertEqual(text.count("assigned_count=3"), 2)
            self.assertIn("recount_mode=balanced", text)
            for mode in ("lenient", "balanced", "strict"):
                result = subprocess.run([sys.executable, str(ROOT / "itdmapper.py"), str(bam),
                    f"--variant-alignments={mode}"], capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stderr)
                calls = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
                self.assertEqual(calls[0]["counts"], "6")
                self.assertEqual(calls[0]["coverage"], "26")
                text = (Path(tmp) / "reads_variant001_alignments.txt").read_text()
                self.assertIn(f"recount_mode={mode}", text)
                with (Path(tmp) / "reads_reassignment_calls.tsv").open() as handle:
                    summary = next(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(int(summary["coverage"]), int(summary["baseline_coverage"])
                    - int(summary["lost_coverage"]) + int(summary["rescued_coverage"]))

    def test_fastq_wt_only_writes_empty_call_audit(self):
        from itdmapper import Flt3Reference
        ref = Flt3Reference.load(ROOT / "annotation" / "FLT3.fa")
        start, end = ref.target_fetch_interval()
        wt = ref.extract_interval(start, end)
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "wt.fastq"
            fastq.write_text("".join(f"@r{i}\n{wt}\n+\n{'I' * len(wt)}\n" for i in range(5)))
            result = subprocess.run([sys.executable, str(ROOT / "itdmapper.py"), str(fastq),
                                     "--recalculate-vaf"], capture_output=True, text=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), "NO_FLT3_ITD")
            with (Path(tmp) / "wt_reassignment_calls.tsv").open() as handle:
                self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), [])
            with (Path(tmp) / "wt_reassignment.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["allele"], "WT")
            self.assertEqual(rows[0]["count"], "5")

    def test_fastq_end_to_end_and_final_filtering_after_discovery(self):
        from itdmapper import Flt3Reference
        ref = Flt3Reference.load(ROOT / "annotation" / "FLT3.fa")
        start, end = ref.target_fetch_interval()
        wt = ref.extract_interval(start - 40, end + 40).upper()
        mutant = wt[:120] + wt[90:120] + wt[120:]
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "reads.fastq"
            fastq.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n"
                for i, seq in enumerate([wt] * 20 + [mutant] * 5 + [mutant[:165]] * 2)))
            result = subprocess.run([sys.executable, str(ROOT / "itdmapper.py"), str(fastq),
                "--recalculate-vaf", "--min-total-reads", "6", "--min-read-copies", "1"],
                capture_output=True, text=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            calls = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
            self.assertEqual(len(calls), 1, result.stdout)
            self.assertEqual(int(calls[0]["counts"]), 7)
            self.assertEqual(int(calls[0]["coverage"]), 27)
            audit = Path(tmp) / "reads_reassignment.tsv"
            with audit.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(sum(int(row["count"]) for row in rows), 27)

    def test_mode_cli_defaults_aliases_and_conflicts(self):
        self.assertEqual(resolve_args(parse_args(["sample.bam", "--recalculate-vaf"])).recount_mode, "balanced")
        self.assertFalse(resolve_args(parse_args(["--variant-alignments", "sample.bam"])).recalculate_vaf)
        for mode in ("lenient", "balanced", "strict"):
            args = resolve_args(parse_args(["sample.bam", f"--variant-alignments={mode}"]))
            self.assertTrue(args.recalculate_vaf)
            self.assertEqual(args.recount_mode, mode)
        for options in [["--variant-alignments=strict", "--no-recalculate-vaf"],
                        ["--variant-alignments=strict", "--recount-mode", "lenient"],
                        ["--recount-mode", "balanced"]]:
            with self.assertRaises(ValueError):
                resolve_args(parse_args(["sample.bam"] + options))
        with tempfile.TemporaryDirectory() as tmp:
            settings = Path(tmp) / "custom.toml"
            settings.write_text('[reassignment]\nenabled = true\nmode = "strict"\nmin_score_margin = 7\n')
            args = resolve_args(parse_args(["sample.bam", "--settings", str(settings)]))
            self.assertEqual(args.reassignment_min_score_margin, 7)
            self.assertEqual(args.reassignment_min_query_fraction, 0.9)
            args = resolve_args(parse_args(["sample.bam", "--settings", str(settings),
                                           "--reassignment-min-query-fraction", "0.6"]))
            self.assertEqual(args.reassignment_min_query_fraction, 0.6)
            settings.write_text('[reassignment]\nmode = "invented"\n')
            with self.assertRaisesRegex(ValueError, "mode"):
                resolve_args(parse_args(["sample.bam", "--settings", str(settings)]))

    def test_ambiguous_reads_remain_coverage_in_balanced(self):
        buckets = Counter({(self.wt, "+", "amp"): 9, (self.mutant, "+", "amp"): 1})
        baseline = {("amp", self.wt): (0, 240), ("amp", self.mutant): (0, 240)}
        # An artificially huge margin also admits the gapped ITD alignment for
        # WT reads. Balanced conservatively requires *all* plausible paths to
        # establish site coverage; lenient may retain the accepted baseline.
        for mode, expected_depth in [("lenient", 10), ("balanced", 1), ("strict", 0)]:
            audit, coverage, assignments = recount_sequences(buckets, self.references,
                self.candidates, self.aligner, baseline, mode=mode, min_score_margin=1000)
            self.assertFalse(assignments)
            self.assertEqual(sum(r["count"] for r in coverage if r["coverage_eligible"]), expected_depth)
            self.assertEqual(sum(r["count"] for r in audit), 10)
            self.assertTrue(all(r["reason"] == "ambiguous_score" for r in audit))

    def test_poor_fits_excluded_in_every_mode(self):
        sequence = "N" * 240
        for mode in ("lenient", "balanced", "strict"):
            audit, coverage, assignments = recount_sequences(Counter({(sequence, "+", "amp"): 4}),
                self.references, self.candidates, self.aligner, {("amp", sequence): (0, 240)}, mode=mode)
            self.assertFalse(assignments)
            self.assertFalse(coverage[0]["coverage_eligible"])
            self.assertEqual(coverage[0]["coverage_reason"], "poor_alignment_all_alternatives")

    def test_flank_mismatch_can_cover_without_exact_anchor(self):
        sequence = self.wt[:119] + ("A" if self.wt[119] != "A" else "C") + self.wt[120:]
        for mode, expected in [("balanced", True), ("strict", False)]:
            audit, coverage, assignments = recount_sequences(Counter({(sequence, "+", "amp"): 1}),
                self.references, self.candidates, self.aligner, {}, mode=mode)
            self.assertEqual(coverage[0]["coverage_eligible"], expected)

    def test_lenient_baseline_preservation_and_amplicon_constraint(self):
        short = self.wt[:119]
        references = self.references + build_references([("other", self.wt, 0)], self.candidates)
        for mode, expected in [("lenient", True), ("balanced", False)]:
            audit, coverage, assignments = recount_sequences(Counter({(short, "+", "amp"): 2}),
                references, self.candidates, self.aligner, {("amp", short): (0, 240)}, mode=mode)
            self.assertEqual(coverage[0]["coverage_eligible"], expected)
            self.assertEqual(audit[0]["amplicon"], "amp")
            self.assertEqual(len(coverage), 1)

    def test_unrelated_candidate_does_not_erase_coverage(self):
        # A read covering the main site ties WT with a candidate outside the read.
        short = self.wt[:155]
        candidates = dict(self.candidates, other=ReferenceEdit(190, 190, "GATTACA"))
        refs = build_references([("amp", self.wt, 0)], candidates)
        audit, coverage, assignments = recount_sequences(Counter({(short, "+", "amp"): 10}),
            refs, candidates, self.aligner, {}, mode="balanced")
        target = next(row for row in coverage if row["HGVS"] == "ITD")
        self.assertEqual(audit[0]["status"], "ambiguous")
        self.assertTrue(target["coverage_eligible"])
        self.assertFalse(target["unresolved"])

    def test_support_and_coverage_reconcile_for_rescued_fragment(self):
        audit, coverage, assignments = recount_sequences(Counter({(self.mutant, "+", ""): 5}),
            self.references, self.candidates, self.aligner, {}, mode="balanced")
        self.assertEqual(coverage[0]["coverage_reason"], "rescued")
        self.assertTrue(coverage[0]["supporting"])
        self.assertTrue(coverage[0]["coverage_eligible"])
        self.assertEqual(sum(n for match, strand, n in assignments), 5)

    def test_boundary_rescue_recovers_truncated_itd_in_all_modes(self):
        refs = build_references([("short", self.wt[:120], 0), ("long", self.wt, 0)], self.candidates)
        for mode in ("lenient", "balanced", "strict"):
            audit, coverage, assignments = recount_sequences(
                Counter({(self.mutant, "+", "short"): 7}), refs, self.candidates,
                self.aligner, {("short", self.mutant): (0, 120)}, mode=mode)
            self.assertEqual(audit[0]["original_amplicon"], "short")
            self.assertEqual(audit[0]["amplicon"], "long")
            self.assertTrue(audit[0]["amplicon_boundary_rescue"])
            self.assertEqual(audit[0]["status"], "assigned")
            self.assertEqual(sum(n for match, strand, n in assignments), 7)
            self.assertEqual(sum(r["count"] for r in coverage if r["supporting"]), 7)
            self.assertEqual(sum(r["count"] for r in coverage if r["coverage_eligible"]), 7)
            self.assertEqual(coverage[0]["coverage_reason"], "amplicon_boundary_rescue")

    def test_boundary_rescue_requires_confidence_and_never_moves_wt(self):
        refs = build_references([("short", self.wt[:120], 0), ("long", self.wt, 0)], self.candidates)
        for sequence, overrides in [(self.wt[:120], {}), ("N" * 270, {}),
                                    (self.mutant, {"min_score_margin": 10000})]:
            audit, coverage, assignments = recount_sequences(
                Counter({(sequence, "+", "short"): 2}), refs, self.candidates,
                self.aligner, {}, mode="balanced", **overrides)
            self.assertEqual(audit[0]["amplicon"], "short")
            self.assertFalse(audit[0]["amplicon_boundary_rescue"])

    def test_boundary_rescue_does_not_move_internal_variants(self):
        refs = build_references([("short", self.wt[:180], 0), ("long", self.wt, 0)], self.candidates)
        audit, coverage, assignments = recount_sequences(
            Counter({(self.mutant, "+", "short"): 1}), refs, self.candidates,
            self.aligner, {}, mode="lenient")
        self.assertEqual(audit[0]["amplicon"], "short")
        self.assertFalse(audit[0]["amplicon_boundary_rescue"])

    def test_boundary_rescue_requires_containment_and_supports_left_boundary(self):
        refs = build_references([("short", self.wt[120:], 120),
                                 ("overlap", self.wt[:180], 0)], self.candidates)
        audit, coverage, assignments = recount_sequences(
            Counter({(self.mutant, "+", "short"): 1}), refs, self.candidates,
            self.aligner, {}, mode="lenient")
        self.assertFalse(audit[0]["amplicon_boundary_rescue"])
        refs += build_references([("long", self.wt, 0)], self.candidates)
        audit, coverage, assignments = recount_sequences(
            Counter({(self.mutant, "+", "short"): 1}), refs, self.candidates,
            self.aligner, {}, mode="strict")
        self.assertEqual(audit[0]["amplicon"], "long")
        self.assertTrue(audit[0]["amplicon_boundary_rescue"])

    def test_boundary_rescue_records_source_coverage_loss_without_double_counting(self):
        candidates = dict(self.candidates, other=ReferenceEdit(60, 60, "GATTACA"))
        refs = build_references([("short", self.wt[:120], 0), ("long", self.wt, 0)], candidates)
        audit, coverage, assignments = recount_sequences(
            Counter({(self.mutant, "+", "short"): 4}), refs, candidates,
            self.aligner, {("short", self.mutant): (0, 120)}, mode="balanced")
        source = next(r for r in coverage if r["HGVS"] == "other" and r["amplicon"] == "short")
        self.assertTrue(source["baseline_covered"])
        self.assertFalse(source["coverage_eligible"])
        self.assertEqual(source["coverage_reason"], "transferred_to_containing_amplicon")
        self.assertEqual(sum(r["count"] for r in coverage if r["HGVS"] == "other" and
                             r["coverage_eligible"]), 4)

    def test_boundary_rescue_does_not_borrow_baseline_from_identical_sequence(self):
        refs = build_references([("short", self.wt[:120], 0), ("long", self.wt, 0)], self.candidates)
        audit, coverage, assignments = recount_sequences(
            Counter({(self.mutant, "+", "short"): 4, (self.mutant, "+", "long"): 3}),
            refs, self.candidates, self.aligner, {("long", self.mutant): (0, 240)}, mode="balanced")
        self.assertEqual(sum(r["count"] for r in coverage if r["baseline_covered"]), 3)
        self.assertEqual(sum(r["count"] for r in coverage if r["coverage_eligible"]), 7)
        self.assertEqual(sum(r["count"] for r in coverage if r["coverage_eligible"] and
                             not r["baseline_covered"]), 4)


if __name__ == "__main__":
    unittest.main()