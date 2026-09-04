import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class TestRunner(unittest.TestCase):
    def test_detects_synthetic_tandem_duplication(self):
        random.seed(4)
        ref = "".join(random.choice("ACGT") for _ in range(120))
        duplicated = ref[55:73]
        mutant = ref[:73] + duplicated + ref[73:]

        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            (work / "reference.txt").write_text(ref + "\n")
            with (work / "annotation.tsv").open("w") as handle:
                handle.write("amplicon_bp\tregion\tchr13_bp\ttranscript_bp\tprotein_as\n")
                for i in range(len(ref)):
                    handle.write(f"{i+1}\texon1\t{i+1}\t{i+1}\t{(i // 3) + 1}\n")
            with (work / "reads.tsv").open("w") as handle:
                handle.write("Sequence\tCounts\tSeqLength\n")
                handle.write(f"{ref}\t90\t{len(ref)}\n")
                handle.write(f"{mutant}\t10\t{len(mutant)}\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "mergeitd_runner.py"),
                    "--sample", "synthetic",
                    "--reads", str(work / "reads.tsv"),
                    "--reference", str(work / "reference.txt"),
                    "--annotation", str(work / "annotation.tsv"),
                    "--output-dir", str(work / "out"),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr + result.stdout)
            calls = (work / "out" / "filtered_mut_vaf.csv").read_text()
            self.assertIn("55_72dup", calls)
            self.assertIn("c.56_73dup", calls)


if __name__ == "__main__":
    unittest.main()
