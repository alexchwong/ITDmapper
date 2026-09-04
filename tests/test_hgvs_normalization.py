import csv
import tempfile
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from hgvs_normalize import (
    ReferenceEdit,
    apply_edit,
    classify_edit,
    normalize_3prime,
    normalize_mergeitd_name,
)
from itdmapper import Amplicon, Flt3Reference, read_calls


def simple_coord(idx: int) -> str:
    return str(idx + 1)


def synthetic_flt3() -> Flt3Reference:
    # c.1835-c.1837, 90-bp intron 14, c.1838-c.1840.
    intron = "C" * 90
    sequence = "GTA" + intron.lower() + "CCT"
    c_by_index = [1835, 1836, 1837] + [None] * len(intron) + [1838, 1839, 1840]
    index_by_c = {c: i for i, c in enumerate(c_by_index) if c is not None}
    region_by_index = ["exon14"] * 3 + ["intron14"] * len(intron) + ["exon15"] * 3
    return Flt3Reference(
        chrom="chr13",
        genomic_start=1001,
        genomic_end=1000 + len(sequence),
        strand="+",
        sequence=sequence,
        c_by_index=c_by_index,
        index_by_c=index_by_c,
        region_by_index=region_by_index,
    )


class TestHgvsNormalization(unittest.TestCase):
    def test_homopolymer_insertion_moves_to_most_3prime_and_becomes_dup(self):
        ref = "CCAAAAAGG"
        raw = ReferenceEdit(3, 3, "A")
        normalized = normalize_3prime(ref, raw)
        self.assertEqual(normalized, ReferenceEdit(7, 7, "A"))
        self.assertEqual(classify_edit(ref, normalized), "dup")
        self.assertEqual(apply_edit(ref, raw), apply_edit(ref, normalized))

    def test_tandem_repeat_shift_can_rotate_inserted_sequence(self):
        ref = "CATATATG"
        raw = ReferenceEdit(2, 2, "TA")
        normalized = normalize_3prime(ref, raw)
        self.assertGreater(normalized.start, raw.start)
        self.assertEqual(apply_edit(ref, raw), apply_edit(ref, normalized))

    def test_normalization_is_idempotent(self):
        ref = "CCAAAAAGG"
        once = normalize_3prime(ref, ReferenceEdit(3, 3, "A"))
        twice = normalize_3prime(ref, once)
        self.assertEqual(once, twice)

    def test_normalization_uses_full_reference_beyond_amplicon_end(self):
        full_ref = "CCAAAAAGG"
        amplicon_ref = full_ref[:5]
        result = normalize_mergeitd_name(
            mutation_name="3_4ins[1]A",
            amplicon_reference=amplicon_ref,
            full_reference=full_ref,
            full_reference_offset=0,
            coord_for_index=simple_coord,
        )
        self.assertEqual(result.edit.start, 7)
        self.assertEqual(result.hgvs, "c.7dup")

    def test_intronic_insertion_uses_adjacent_flanking_bases_without_plus_one_shift(self):
        ref = synthetic_flt3()
        # Intron starts at index 3: index 23 is c.1837+21 and index 24 is +22.
        self.assertEqual(ref.hgvs_coord_for_index(23), "1837+21")
        self.assertEqual(ref.hgvs_coord_for_index(24), "1837+22")
        result = normalize_mergeitd_name(
            mutation_name="24_25ins[2]TT",
            amplicon_reference=ref.sequence,
            full_reference=ref.sequence,
            full_reference_offset=0,
            coord_for_index=ref.hgvs_coord_for_index,
        )
        self.assertEqual(result.hgvs, "c.1837+21_1837+22ins[2]TT")
        self.assertNotIn("1837+22_1837+23", result.hgvs)

    def test_exonic_insertion_keeps_correct_interbase_boundary(self):
        ref = "AACGGA"
        result = normalize_mergeitd_name(
            mutation_name="2_3ins[2]TT",
            amplicon_reference=ref,
            full_reference=ref,
            full_reference_offset=0,
            coord_for_index=lambda idx: str(1829 + idx),
        )
        self.assertEqual(result.hgvs, "c.1830_1831ins[2]TT")
        self.assertNotIn("1831_1832", result.hgvs)

    def test_read_calls_normalizes_before_final_hgvs_output(self):
        ref = synthetic_flt3()
        amplicon = Amplicon(start=1000, end=1000 + len(ref.sequence), support=10)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "filtered_mut_vaf.csv"
            fieldnames = [
                "netInsert", "counts", "vaf_percent", "insertPos", "insertRegion",
                "coverage", "name", "HGVS", "co_mutations",
            ]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow({
                    "netInsert": "2",
                    "counts": "10",
                    "vaf_percent": "5.0",
                    "insertPos": "24",
                    "insertRegion": "intron14",
                    "coverage": "200",
                    "name": "24_25ins[2]TT",
                    "HGVS": "c.1837+22_1837+23ins[2]TT",
                    "co_mutations": "{}",
                })

            calls = read_calls(path, amplicon, "chr13", ref)
        self.assertEqual(
            calls[0]["HGVS"],
            "NM_004119.3:c.1837+21_1837+22ins[2]TT",
        )

    def test_read_calls_maps_minus_strand_amplicon_coordinates_to_full_reference(self):
        sequence = "ACGTCAGTCCGATGCTAACG"
        c_by_index = [1705 + i for i in range(len(sequence))]
        ref = Flt3Reference(
            chrom="chr13",
            genomic_start=1001,
            genomic_end=1020,
            strand="-",
            sequence=sequence,
            c_by_index=c_by_index,
            index_by_c={c: i for i, c in enumerate(c_by_index)},
            region_by_index=["exon14"] * len(sequence),
        )
        # Genomic interval 1007..1016 corresponds to full transcript indices 4..13.
        amplicon = Amplicon(start=1006, end=1016, support=10)
        self.assertEqual(ref.extract_interval(amplicon.start, amplicon.end), sequence[4:14])

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "filtered_mut_vaf.csv"
            fieldnames = [
                "netInsert", "counts", "vaf_percent", "insertPos", "insertRegion",
                "coverage", "name", "HGVS", "co_mutations",
            ]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow({
                    "netInsert": "2", "counts": "10", "vaf_percent": "5.0",
                    "insertPos": "3", "insertRegion": "exon14", "coverage": "200",
                    "name": "3_4ins[2]AA", "HGVS": "c.raw", "co_mutations": "{}",
                })
            calls = read_calls(path, amplicon, "chr13", ref)

        # Local boundary 3 + full offset 4 = full boundary 7, flanked by c.1711/c.1712.
        self.assertEqual(calls[0]["HGVS"], "NM_004119.3:c.1711_1712ins[2]AA")

    def test_existing_duplication_coordinate_semantics_are_preserved(self):
        duplicated = "ACGTTGCAACCTGATCGGTAC"
        ref = duplicated + "TGGCA"
        result = normalize_mergeitd_name(
            mutation_name="0_20dup",
            amplicon_reference=ref,
            full_reference=ref,
            full_reference_offset=0,
            coord_for_index=lambda idx: str(1784 + idx),
        )
        self.assertEqual(result.hgvs, "c.1784_1804dup")
        self.assertEqual(apply_edit(ref, result.edit), duplicated + duplicated + ref[21:])


if __name__ == "__main__":
    unittest.main()
