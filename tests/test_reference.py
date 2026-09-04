import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from itdmapper import Flt3Reference


class TestReference(unittest.TestCase):
    def test_repo_flt3_reference_coordinates(self):
        fasta = ROOT / "annotation" / "FLT3.fa"
        if not fasta.exists():
            self.skipTest("Run from the ITDmapper repository with annotation/FLT3.fa present")
        ref = Flt3Reference.load(fasta)
        self.assertEqual(ref.chrom, "chr13")
        self.assertEqual(ref.strand, "-")
        self.assertEqual(ref.genomic_pos_for_c(1705), 28608351)
        self.assertEqual(ref.genomic_pos_for_c(1837), 28608219)
        self.assertEqual(ref.genomic_pos_for_c(1838), 28608128)
        self.assertEqual(ref.genomic_pos_for_c(1942), 28608024)
        self.assertEqual(ref.region_by_index[ref.index_by_c[1705]], "exon14")
        self.assertEqual(ref.region_by_index[ref.index_by_c[1838]], "exon15")

    def test_repo_tsv_agrees_with_generated_annotation(self):
        fasta = ROOT / "annotation" / "FLT3.fa"
        tsv = ROOT / "annotation" / "amplicon_kayser.tsv"
        if not fasta.exists() or not tsv.exists():
            self.skipTest("Run from the ITDmapper repository with annotations present")

        ref = Flt3Reference.load(fasta)
        rows = ref.annotation_rows(28608021, 28608352)  # TSV 1-based 28608022..28608352
        by_genomic = {int(row["chr13_bp"]): row for row in rows}

        import csv
        with tsv.open() as handle:
            for old in csv.DictReader(handle, delimiter="\t"):
                g = int(old["chr13_bp"])
                generated = by_genomic[g]
                self.assertEqual(str(generated["transcript_bp"]), old["transcript_bp"])
                if "exon14" in old["region"]:
                    self.assertEqual(generated["region"], "exon14")
                elif "intron14" in old["region"]:
                    self.assertEqual(generated["region"], "intron14")
                elif "exon15" in old["region"]:
                    self.assertEqual(generated["region"], "exon15")


if __name__ == "__main__":
    unittest.main()
