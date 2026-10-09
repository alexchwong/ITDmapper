import csv
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from Bio import SeqIO
from itdmapper import (
    Flt3Reference,
    collect_variant_sequences,
    parse_args,
    pysam,
    reverse_complement,
    strand_sequence_counts,
    write_variant_fastas,
)


class TestVariantFasta(unittest.TestCase):
    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args(["sample.bam"]).variant_fasta)
        self.assertTrue(parse_args(["sample.bam", "--variant-fasta"]).variant_fasta)

    def test_counts_preserve_original_strand_in_common_orientation(self):
        reads = [SimpleNamespace(query_sequence="AACCGT", is_reverse=False)] * 3
        reads += [SimpleNamespace(query_sequence="AACCGT", is_reverse=True)] * 2
        self.assertEqual(
            strand_sequence_counts(reads, "-"),
            Counter({("ACGGTT", "+"): 3, ("ACGGTT", "-"): 2}),
        )

    def test_assignment_tokens_filters_and_merging(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "alignClasses.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["Sequence", "Aligned", "HGVS"])
                writer.writeheader()
                writer.writerows([
                    {"Sequence": "ACGT", "Aligned": True, "HGVS": "mutation;other"},
                    {"Sequence": "AAAA", "Aligned": True, "HGVS": "mutation_suffix"},
                    {"Sequence": "CCCC", "Aligned": False, "HGVS": "mutation"},
                ])
            calls = [{"name": "mutation", "HGVS": "variant1"},
                     {"name": "other", "HGVS": "variant2"},
                     {"name": "other", "HGVS": "variant1"}]
            counts = Counter({("ACGT", "+"): 4, ("ACGT", "-"): 2,
                              ("AAAA", "+"): 9, ("CCCC", "+"): 8})
            variants = {}
            collect_variant_sequences(path, calls, counts, variants)
            self.assertEqual(variants["variant1"], Counter({("ACGT", "+"): 4, ("ACGT", "-"): 2}))
            self.assertEqual(variants["variant2"], variants["variant1"])
            collect_variant_sequences(path, calls, counts, variants)
            self.assertEqual(variants["variant1"][("ACGT", "+")], 8)

    def test_writer_names_order_and_no_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "sample.bam"
            calls = [{"HGVS": "high"}, {"HGVS": "low"}]
            variants = {"high": Counter({("ACGT", "-"): 20, ("ACGT", "+"): 2,
                                         ("AAAA", "+"): 5}),
                        "low": Counter({("CCCC", "-"): 1})}
            paths = write_variant_fastas(bam, calls, variants)
            self.assertEqual([p.name for p in paths], ["sample_variant001.fasta", "sample_variant002.fasta"])
            with paths[0].open() as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            self.assertEqual([r.id for r in records], [
                "sequence001_count5_orientation+", "sequence002_count2_orientation+",
                "sequence003_count20_orientation-",
            ])
            self.assertEqual([str(r.seq) for r in records], ["AAAA", "ACGT", "ACGT"])
            self.assertIn("HGVS=high", records[0].description)
            self.assertEqual(write_variant_fastas(bam, [], {}), [])

    def test_fastq_export_rejected(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "itdmapper.py"), "missing.fastq", "--variant-fasta"],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires BAM input", result.stderr)

    @unittest.skipIf(pysam is None, "pysam is not installed")
    def test_synthetic_bam_end_to_end(self):
        reference = Flt3Reference.load(ROOT / "annotation" / "FLT3.fa")
        start, end = reference.target_fetch_interval()
        wt = reference.extract_interval(start, end).upper()
        mutants = [(70, 18, 12), (140, 24, 8)]
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "synthetic.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": reference.chrom, "LN": reference.genomic_end + 1000}]}
            expected = {}
            with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
                groups = [(wt, [(0, len(wt))], 80)]
                for pos, size, count in mutants:
                    seq = wt[:pos] + wt[pos - size:pos] + wt[pos:]
                    genomic_pos = len(wt) - pos
                    groups.append((seq, [(0, genomic_pos), (1, size), (0, pos)], count))
                    expected[seq] = count
                for group, (seq, cigar, count) in enumerate(groups):
                    for i in range(count):
                        read = pysam.AlignedSegment()
                        read.query_name = f"group{group}_read{i}"
                        read.query_sequence = reverse_complement(seq)
                        read.flag = 16 if i % 2 else 0
                        read.reference_id = 0
                        read.reference_start = start
                        read.mapping_quality = 60
                        read.cigartuples = cigar
                        read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                        handle.write(read)
            pysam.index(str(bam))
            command = [sys.executable, str(ROOT / "itdmapper.py"), str(bam)]
            baseline = subprocess.run(command, capture_output=True, text=True, timeout=120)
            self.assertEqual(baseline.returncode, 0, baseline.stderr)
            self.assertEqual(list(Path(tmp).glob("*.fasta")), [])
            result = subprocess.run(command + ["--variant-fasta"], capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, baseline.stdout)
            calls = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
            self.assertEqual(len(calls), 2, result.stdout)
            self.assertGreater(float(calls[0]["vaf_percent"]), float(calls[1]["vaf_percent"]))
            for number, call in enumerate(calls, start=1):
                path = Path(tmp) / f"synthetic_variant{number:03d}.fasta"
                with path.open() as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                self.assertEqual(len(records), 2)
                seq = str(records[0].seq)
                self.assertEqual(str(records[1].seq), seq)
                self.assertIn(f"count{expected[seq] // 2}_orientation+", records[0].id)
                self.assertIn(f"count{expected[seq] // 2}_orientation-", records[1].id)
                self.assertIn(f"HGVS={call['HGVS']}", records[0].description)


if __name__ == "__main__":
    unittest.main()