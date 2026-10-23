import csv
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import Mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itdmapper import (
    _make_aligner,
    collect_variant_sequences,
    parse_args,
    read_variant_references,
    resolve_args,
    variant_sequence_records,
    write_variant_alignments,
    write_variant_fastas,
)


class TestVariantAlignments(unittest.TestCase):
    def test_flag_and_fastq_validation(self):
        self.assertFalse(parse_args(["sample.bam"]).variant_alignments)
        self.assertTrue(parse_args(["sample.bam", "--variant-alignments"]).variant_alignments)
        result = subprocess.run(
            [sys.executable, str(ROOT / "itdmapper.py"), "missing.fastq", "--variant-alignments"],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--variant-alignments requires BAM input", result.stderr)

    def test_provenance_collection_across_amplicons(self):
        variants, provenance = {}, {}
        counts = Counter({("ACGT", "+"): 3, ("ACGT", "-"): 2})
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "alignClasses.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["Sequence", "Aligned", "HGVS"])
                writer.writeheader()
                writer.writerows([
                    {"Sequence": "ACGT", "Aligned": True, "HGVS": "mutation;other"},
                    {"Sequence": "AAAA", "Aligned": False, "HGVS": "mutation"},
                ])
            calls = [{"name": "mutation", "HGVS": "variant"},
                     {"name": "other", "HGVS": "variant"}]
            for source in [("chr13:1-4", "ACGT"), ("chr13:2-5", "CGTA")]:
                collect_variant_sequences(path, calls, counts, variants, provenance,
                                          {"mutation": source, "other": source})
            self.assertEqual(variants["variant"], Counter({("ACGT", "+"): 6, ("ACGT", "-"): 4}))
            self.assertEqual(provenance["variant"][("ACGT", "+")], Counter({
                ("chr13:1-4", "ACGT"): 3, ("chr13:2-5", "CGTA"): 3,
            }))
            with self.assertRaisesRegex(ValueError, "Mutant references"):
                collect_variant_sequences(path, calls, counts, variants, provenance)

    def test_emitted_mutant_references_follow_raw_call_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            filtered = Path(tmp) / "filtered_mut_vaf.csv"
            filtered.write_text("name,HGVS\n2_3dup,variant1\n4_5ins[2]TT,variant2\n")
            fasta = Path(tmp) / "sample_ampliconome.fa"
            fasta.write_text(">WildType\nACGTAC\n>mutation_0\nACGTGTAC\n>mutation_1\nACGTTTAC\n")
            self.assertEqual(read_variant_references(filtered, fasta, "chr13:1-6"), {
                "2_3dup": ("chr13:1-6", "ACGTGTAC"),
                "4_5ins[2]TT": ("chr13:1-6", "ACGTTTAC"),
            })
            fasta.write_text(">WildType\nACGTAC\n")
            with self.assertRaisesRegex(ValueError, "Missing mutant reference"):
                read_variant_references(filtered, fasta, "chr13:1-6")

    def test_assignment_selects_matching_mutant_not_wt_or_another_variant(self):
        with tempfile.TemporaryDirectory() as tmp:
            alignment = Path(tmp) / "alignClasses.csv"
            alignment.write_text("Sequence,Aligned,HGVS\nACGTGTAC,True,2_3dup\n")
            variants, provenance = {}, {}
            collect_variant_sequences(
                alignment,
                [{"name": "2_3dup", "HGVS": "variant1"},
                 {"name": "4_5ins[2]TT", "HGVS": "variant2"}],
                Counter({("ACGTGTAC", "+"): 3}), variants, provenance,
                {"2_3dup": ("chr13:1-6", "ACGTGTAC"),
                 "4_5ins[2]TT": ("chr13:1-6", "ACGTTTAC")},
            )
            aligner = Mock(wraps=_make_aligner(5, -15, -36, -0.5))
            paths = write_variant_alignments(Path(tmp) / "sample.bam", [{"HGVS": "variant1"}],
                                             variants, provenance, aligner)
            self.assertEqual(str(aligner.align.call_args.args[1].seq), "ACGTGTAC")
            self.assertIn("reference_type=FLT3-ITD", paths[0].read_text())
            self.assertNotIn("variant2", variants)

    def test_native_rendering_and_runtime_scoring(self):
        cases = [
            ("ACGCAATGTG", "ACGCAATGTG"),  # identical
            ("ACGCAATGTG", "ACGCAGTGTG"),  # mismatch
            ("ACGCAATGTG", "ACGCTGTG"),    # insertion
            ("ACGCTGTG", "ACGCAATGTG"),    # deletion
            ("TTACGCAATGTG", "ACGCAATGTG"),  # terminal overhang
            ("ACGCAATGTG" * 15, "ACGCAATGTG" * 15),  # wrapped
            ("AAAAAA", "AAAAA"),  # equivalent optimal alignments
        ]
        for scores in [[], ["--match-score", "3", "--mismatch-score", "-2",
                            "--gap-open-score", "-4", "--gap-extend-score", "-1"]]:
            args = resolve_args(parse_args(["sample.bam"] + scores))
            aligner = _make_aligner(args.match_score, args.mismatch_score,
                                    args.gap_open_score, args.gap_extend_score)
            for sequence, reference in cases:
                with self.subTest(scores=scores, sequence=sequence):
                    with tempfile.TemporaryDirectory() as tmp:
                        calls = [{"HGVS": "variant"}]
                        variants = {"variant": Counter({(sequence, "+"): 2})}
                        provenance = {"variant": {(sequence, "+"): Counter({
                            ("chr13:1-10", reference): 2,
                        })}}
                        paths = write_variant_alignments(Path(tmp) / "sample.bam", calls,
                                                         variants, provenance, aligner)
                        expected = str(aligner.align(SeqRecord(Seq(sequence), id="read"),
                                                    SeqRecord(Seq(reference), id="reference"))[-1])
                        self.assertIn(expected, paths[0].read_text())
                        self.assertEqual(paths[0].name, "sample_variant001_alignments.txt")

    def test_order_counts_multiple_references_and_cache(self):
        variants = {"high": Counter({("ACGT", "-"): 20, ("ACGT", "+"): 2,
                                     ("AAAA", "+"): 5})}
        variants["low"] = variants["high"].copy()
        sources = {
            ("AAAA", "+"): Counter({("chr13:1-4", "AAAA"): 5}),
            ("ACGT", "+"): Counter({("chr13:1-4", "ACGT"): 1,
                                      ("chr13:2-5", "CGTA"): 1}),
            ("ACGT", "-"): Counter({("chr13:1-4", "ACGT"): 20}),
        }
        calls = [{"HGVS": "high"}, {"HGVS": "low"}]
        aligner = Mock(wraps=_make_aligner(5, -15, -36, -0.5))
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "sample.bam"
            fastas = write_variant_fastas(bam, calls, variants)
            paths = write_variant_alignments(bam, calls, variants,
                                             {"high": sources, "low": sources}, aligner)
            self.assertEqual(aligner.align.call_count, 3)
            for fasta, path, call in zip(fastas, paths, calls):
                text = path.read_text()
                ids = [line for line in text.splitlines() if line.startswith("sequence")]
                records = variant_sequence_records(call["HGVS"], variants[call["HGVS"]])
                self.assertEqual(ids, [record.id for _, record in records])
                self.assertEqual(fasta.read_text(), "".join(
                    f">{record.id} {record.description}\n{record.seq}\n" for _, record in records
                ))
                self.assertIn("amplicon=chr13:2-5 assigned_count=1", text)
                self.assertIn("amplicon=chr13:1-4 assigned_count=20", text)
            self.assertEqual(write_variant_alignments(bam, [], {}, {}, aligner), [])
            with self.assertRaisesRegex(ValueError, "Incomplete mutant reference provenance"):
                write_variant_alignments(bam, calls, variants, {}, aligner)


if __name__ == "__main__":
    unittest.main()