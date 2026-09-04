#!/usr/bin/env python3
"""External adapter that invokes the verbatim mergeITD core."""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import Align

from itdmapper_settings import load_settings, validate_settings
from mergeitd_core import alignITD, annotateCoords


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--settings", type=Path, default=None)
    parser.add_argument("--match-score", type=float, default=None)
    parser.add_argument("--mismatch-score", type=float, default=None)
    parser.add_argument("--gap-open-score", type=float, default=None)
    parser.add_argument("--gap-extend-score", type=float, default=None)
    parser.add_argument("--min-aligned-block", type=int, default=None)
    parser.add_argument("--min-reference-fraction", type=float, default=None)
    parser.add_argument("--max-indel-fraction", type=float, default=None)
    parser.add_argument("--min-insert-seq-length", type=int, default=None)
    parser.add_argument("--min-total-reads", type=int, default=None)
    parser.add_argument("--min-vaf", type=float, default=None)
    return parser.parse_args(argv)


def resolve_args(args: argparse.Namespace) -> argparse.Namespace:
    settings = load_settings(args.settings)
    mappings = {
        "match_score": ("alignment", "match"),
        "mismatch_score": ("alignment", "mismatch"),
        "gap_open_score": ("alignment", "gap_open"),
        "gap_extend_score": ("alignment", "gap_extend"),
        "min_aligned_block": ("alignment", "min_aligned_block"),
        "min_reference_fraction": ("alignment", "min_reference_fraction"),
        "max_indel_fraction": ("alignment", "max_indel_fraction"),
        "min_insert_seq_length": ("calling", "min_insert_length"),
        "min_total_reads": ("calling", "min_supporting_reads"),
        "min_vaf": ("calling", "min_vaf_percent"),
    }
    effective = deepcopy(settings)
    for attr, (section, key) in mappings.items():
        if getattr(args, attr) is None:
            setattr(args, attr, settings[section][key])
        effective[section][key] = getattr(args, attr)
    validate_settings(effective, source="CLI/settings")
    return args


def main() -> int:
    args = resolve_args(parse_args())
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prealigns = pd.read_csv(args.reads, sep="\t")
    if prealigns.empty:
        # Keep a valid shape so alignITD can still emit empty result files.
        prealigns = pd.DataFrame(columns=["Sequence", "Counts", "SeqLength"])

    ref = args.reference.read_text().strip().upper()
    anno = pd.read_csv(args.annotation, sep="\t")
    anno = annotateCoords(anno)

    config = {
        "SAMPLE": args.sample,
        "REF": ref,
        "ANNO": anno,
        "PROGRESSBAR": False,
        "PLOT": False,
        "COST_MATCH": args.match_score,
        "COST_MISMATCH": args.mismatch_score,
        "COST_GAPOPEN": args.gap_open_score,
        "COST_GAPEXTEND": args.gap_extend_score,
        "MIN_ALIGN_LEN": args.min_aligned_block,
        "MIN_REF_ALN_FRACTION": args.min_reference_fraction,
        "MAX_FRAC_INDEL": args.max_indel_fraction,
        "MIN_INSERT_SEQ_LENGTH": args.min_insert_seq_length,
        "MIN_TOTAL_READS": args.min_total_reads,
        "MIN_VAF": args.min_vaf,
        "ALIGN_FILE": str(args.output_dir / "alignClasses.csv"),
        "MUTATION_FILE": str(args.output_dir / "mutation_vaf.csv"),
        "MUTATION_FILE_FILTERED": str(args.output_dir / "filtered_mut_vaf.csv"),
        "NETINSERT_FILE": str(args.output_dir / "netInserts_vaf.csv"),
        "OME_FILE": str(args.output_dir / f"{args.sample}_ampliconome.fa"),
        "OUT_COV_FILE": str(args.output_dir / "coverage.txt"),
        "OUT_COV_PLOT": str(args.output_dir / "coverage.png"),
        "STATS_FILE": str(args.output_dir / "stats.txt"),
    }

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = config["COST_MATCH"]
    aligner.mismatch_score = config["COST_MISMATCH"]
    aligner.open_gap_score = config["COST_GAPOPEN"]
    aligner.extend_gap_score = config["COST_GAPEXTEND"]
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    config["ALIGNER"] = aligner

    alignITD(prealigns, config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
