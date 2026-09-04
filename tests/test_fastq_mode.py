import gzip
import random
import tempfile
import unittest
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from itdmapper import (
    Amplicon,
    Flt3Reference,
    assign_sequence_counts_to_amplicons,
    detect_input_format,
    discover_amplicons_from_sequences,
    load_fastq_sequence_counts,
    reverse_complement,
)


def dna(length, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def fake_reference(sequence, strand="+"):
    start1 = 1001
    end1 = start1 + len(sequence) - 1
    return Flt3Reference(
        chrom="chr13",
        genomic_start=start1,
        genomic_end=end1,
        strand=strand,
        sequence=sequence,
        c_by_index=[None] * len(sequence),
        index_by_c={},
        region_by_index=["other"] * len(sequence),
    )


def fastq_record(name, sequence, quality_char="I"):
    return f"@{name}\n{sequence}\n+\n{quality_char * len(sequence)}\n"


class TestFastqMode(unittest.TestCase):
    def test_input_detection_and_override(self):
        self.assertEqual(detect_input_format(Path("x.bam")), "bam")
        self.assertEqual(detect_input_format(Path("x.fastq")), "fastq")
        self.assertEqual(detect_input_format(Path("x.fq.gz")), "fastq")
        self.assertEqual(detect_input_format(Path("odd.name"), "fastq"), "fastq")
        with self.assertRaisesRegex(ValueError, "Unable to infer input format"):
            detect_input_format(Path("x.txt"))

    def test_fastq_quality_filter_then_copy_filter(self):
        seq_a = "ACGTACGTACGT"
        seq_b = "TTTTCCCCAAAA"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq.gz"
            text = (
                fastq_record("a1", seq_a)
                + fastq_record("a2", seq_a)
                + fastq_record("b_low", seq_b, quality_char="!")
                + fastq_record("singleton", "GGGGAAAACCCC")
            )
            with gzip.open(path, "wt") as handle:
                handle.write(text)

            counts, raw_reads, quality_passing = load_fastq_sequence_counts(
                path, min_mean_base_quality=20, min_read_copies=2
            )
            self.assertEqual(raw_reads, 4)
            self.assertEqual(quality_passing, 3)
            self.assertEqual(counts, Counter({seq_a: 2}))

    def test_discovers_multiple_amplicons_from_top_unique_and_reverse_orientation(self):
        ref_seq = dna(700, seed=7)
        ref = fake_reference(ref_seq)
        # genomic 0-based reference coordinates are 1000..1700
        amp1 = (1100, 1350)
        amp2 = (1125, 1380)
        seq1 = ref.extract_interval(*amp1).upper()
        seq2 = ref.extract_interval(*amp2).upper()
        counts = Counter(
            {
                seq1: 100,
                reverse_complement(seq1): 40,
                seq2: 35,
                reverse_complement(seq2): 10,
            }
        )
        amplicons = discover_amplicons_from_sequences(
            counts,
            ref,
            target_interval=(1200, 1300),
            discovery_flank=250,
            top_unique=10,
            endpoint_tolerance=5,
            min_cluster_fraction=0.10,
            min_cluster_support=2,
            match_score=5,
            mismatch_score=-15,
            gap_open_score=-36,
            gap_extend_score=-0.5,
            min_aligned_block=6,
            min_query_fraction=0.5,
            max_indel_fraction=0.7,
        )
        observed = {(a.start, a.end): a.support for a in amplicons}
        self.assertEqual(observed[amp1], 140)
        self.assertEqual(observed[amp2], 45)

    def test_discovery_coordinate_conversion_on_minus_strand(self):
        ref_seq = dna(700, seed=17)
        ref = fake_reference(ref_seq, strand="-")
        expected = (1120, 1370)
        seq = ref.extract_interval(*expected).upper()
        amplicons = discover_amplicons_from_sequences(
            Counter({reverse_complement(seq): 50}),
            ref,
            target_interval=(1200, 1300),
            discovery_flank=250,
            top_unique=10,
            endpoint_tolerance=5,
            min_cluster_fraction=0.10,
            min_cluster_support=2,
            match_score=5,
            mismatch_score=-15,
            gap_open_score=-36,
            gap_extend_score=-0.5,
            min_aligned_block=6,
            min_query_fraction=0.3,
            max_indel_fraction=0.7,
        )
        self.assertEqual((amplicons[0].start, amplicons[0].end), expected)

    def test_assignment_requires_reference_and_query_coverage(self):
        ref_seq = dna(600, seed=21)
        ref = fake_reference(ref_seq)
        amp = Amplicon(start=1100, end=1300, support=100)
        full = ref.extract_interval(amp.start, amp.end).upper()
        partial = full[:100]  # 50% of amplicon, 100% of query
        tailed = full[:100] + dna(300, seed=31)  # query fraction below 0.3
        counts = Counter({full: 10, reverse_complement(partial): 5, tailed: 3})

        groups, accepted, rejected = assign_sequence_counts_to_amplicons(
            counts,
            [amp],
            ref,
            match_score=5,
            mismatch_score=-15,
            gap_open_score=-36,
            gap_extend_score=-0.5,
            min_aligned_block=6,
            min_reference_fraction=0.4,
            min_query_fraction=0.3,
            max_indel_fraction=0.7,
        )
        self.assertEqual(accepted, 15)
        self.assertEqual(rejected, 3)
        self.assertEqual(groups[0][full], 10)
        self.assertEqual(groups[0][partial], 5)

    def test_discovery_fails_when_candidates_do_not_overlap_target(self):
        ref_seq = dna(700, seed=41)
        ref = fake_reference(ref_seq)
        off_target = ref.extract_interval(1010, 1100).upper()
        with self.assertRaisesRegex(ValueError, "Unable to infer"):
            discover_amplicons_from_sequences(
                Counter({off_target: 20}),
                ref,
                target_interval=(1300, 1400),
                discovery_flank=400,
                top_unique=10,
                endpoint_tolerance=5,
                min_cluster_fraction=0.1,
                min_cluster_support=2,
                match_score=5,
                mismatch_score=-15,
                gap_open_score=-36,
                gap_extend_score=-0.5,
                min_aligned_block=6,
                min_query_fraction=0.5,
                max_indel_fraction=0.7,
            )


if __name__ == "__main__":
    unittest.main()
