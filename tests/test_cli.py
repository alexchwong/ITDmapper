import importlib.util
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class TestCLI(unittest.TestCase):
    def test_fixture_bam_runs_end_to_end(self):
        if importlib.util.find_spec("pysam") is None:
            self.skipTest("pysam is not installed")
        bam = ROOT / "fixtures" / "flt3_subset_downsampled.bam"
        bai = ROOT / "fixtures" / "flt3_subset_downsampled.bam.bai"
        fasta = ROOT / "annotation" / "FLT3.fa"
        if not bam.exists() or not bai.exists() or not fasta.exists():
            self.skipTest("Repository fixture/reference files are not present")

        result = subprocess.run(
            [sys.executable, str(ROOT / "itdmapper.py"), str(bam)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr + result.stdout)
        self.assertTrue(
            result.stdout.startswith("NO_FLT3_ITD") or result.stdout.startswith("HGVS\t"),
            msg=result.stdout,
        )


if __name__ == "__main__":
    unittest.main()
