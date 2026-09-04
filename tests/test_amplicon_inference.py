import unittest
from dataclasses import dataclass
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from itdmapper import Amplicon, infer_amplicons, assign_reads_to_amplicons, sequence_counts


@dataclass
class FakeRead:
    reference_start: int
    reference_end: int
    query_sequence: str = "AACCGGTT"
    is_unmapped: bool = False
    is_secondary: bool = False
    is_supplementary: bool = False
    is_paired: bool = False


class TestAmpliconInference(unittest.TestCase):
    def test_modal_endpoint_cluster_ignores_outlier(self):
        reads = [FakeRead(90, 210) for _ in range(10)]
        reads += [FakeRead(91, 210) for _ in range(2)]
        reads += [FakeRead(50, 300)]
        amps = infer_amplicons(reads, (100, 120), (180, 200))
        self.assertEqual(amps[0], Amplicon(90, 210, 12))
        self.assertEqual(len(amps), 1)

    def test_assigns_soft_clipped_read_by_anchored_end(self):
        amps = [Amplicon(90, 210, 10)]
        reads = [FakeRead(90, 210), FakeRead(90, 150)]
        groups = assign_reads_to_amplicons(reads, amps)
        self.assertEqual(len(groups[0]), 2)

    def test_minus_strand_sequence_orientation(self):
        reads = [FakeRead(90, 210, "AACCGT")]
        counts = sequence_counts(reads, "-")
        self.assertEqual(counts["ACGGTT"], 1)


if __name__ == "__main__":
    unittest.main()
