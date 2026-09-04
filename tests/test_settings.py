import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from itdmapper import parse_args, primary_target_reads, resolve_args
from itdmapper_settings import load_settings


@dataclass
class FakeRead:
    mapping_quality: int = 60
    query_qualities: tuple[int, ...] | None = (35, 35, 35, 35)
    is_duplicate: bool = False
    is_qcfail: bool = False
    is_unmapped: bool = False
    is_secondary: bool = False
    is_supplementary: bool = False
    is_paired: bool = False
    query_sequence: str | None = "ACGT"
    reference_end: int | None = 120


class FakeBam:
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, contig, start0, end0):
        return iter(self.reads)


class TestSettings(unittest.TestCase):
    def test_packaged_defaults_preserve_existing_quality_behaviour(self):
        settings = load_settings()
        self.assertEqual(settings["reads"]["min_mapq"], 0)
        self.assertEqual(settings["reads"]["min_mean_base_quality"], 0.0)
        self.assertFalse(settings["reads"]["exclude_duplicates"])
        self.assertFalse(settings["reads"]["exclude_qcfail"])
        self.assertEqual(settings["reads"]["min_read_copies"], 2)
        self.assertEqual(settings["calling"]["min_vaf_percent"], 0.006)

    def test_custom_file_is_partial_and_cli_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "custom.toml"
            path.write_text(
                "[reads]\nmin_mapq = 25\nmin_read_copies = 3\n\n"
                "[calling]\nmin_vaf_percent = 0.05\n"
            )
            args = parse_args([
                "sample.bam",
                "--settings", str(path),
                "--min-vaf", "0.2",
            ])
            args = resolve_args(args)
            self.assertEqual(args.min_mapq, 25)
            self.assertEqual(args.min_read_copies, 3)
            self.assertEqual(args.min_vaf, 0.2)
            self.assertEqual(args.min_insert_seq_length, 6)

    def test_unknown_setting_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.toml"
            path.write_text("[reads]\nmin_mappq = 20\n")
            with self.assertRaisesRegex(ValueError, "Unknown setting"):
                load_settings(path)

    def test_invalid_cli_override_fails_validation(self):
        args = parse_args(["sample.bam", "--min-mapq", "-1"])
        with self.assertRaisesRegex(ValueError, "min_mapq"):
            resolve_args(args)

    def test_quality_filters(self):
        reads = [
            FakeRead(mapping_quality=60, query_qualities=(35, 35, 35, 35)),
            FakeRead(mapping_quality=10, query_qualities=(35, 35, 35, 35)),
            FakeRead(mapping_quality=60, query_qualities=(10, 10, 10, 10)),
            FakeRead(mapping_quality=60, query_qualities=(35, 35, 35, 35), is_duplicate=True),
            FakeRead(mapping_quality=60, query_qualities=(35, 35, 35, 35), is_qcfail=True),
        ]
        kept = primary_target_reads(
            FakeBam(reads),
            "chr13",
            0,
            1000,
            min_mapq=20,
            min_mean_base_quality=25,
            exclude_duplicates=True,
            exclude_qcfail=True,
        )
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0].mapping_quality, 60)

    def test_quality_filters_disabled_preserve_duplicate_and_qcfail_reads(self):
        reads = [
            FakeRead(mapping_quality=0, query_qualities=None, is_duplicate=True),
            FakeRead(mapping_quality=0, query_qualities=None, is_qcfail=True),
        ]
        kept = primary_target_reads(FakeBam(reads), "chr13", 0, 1000)
        self.assertEqual(len(kept), 2)


if __name__ == "__main__":
    unittest.main()
