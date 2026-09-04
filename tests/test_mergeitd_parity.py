import inspect
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import mergeitd_core


def normalized_source(obj):
    lines = inspect.getsource(obj).splitlines()
    return "\n".join(line.rstrip() for line in lines).strip()


class TestMergeITDParity(unittest.TestCase):
    def test_extracted_functions_match_legacy_source(self):
        if not (ROOT / "mergeitd.py").exists():
            self.skipTest("Place changed files in the ITDmapper repository to run parity test")
        import mergeitd

        names = [
            "annotateCoords",
            "save_stats",
            "Mutation",
            "getHGVS",
            "generateMutSeq",
            "getRefLoc",
            "findSNP",
            "alignITD",
            "save_coverage",
            "plot_coverage",
        ]
        for name in names:
            with self.subTest(name=name):
                self.assertEqual(
                    normalized_source(getattr(mergeitd_core, name)),
                    normalized_source(getattr(mergeitd, name)),
                )


if __name__ == "__main__":
    unittest.main()
