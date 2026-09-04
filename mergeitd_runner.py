#!/usr/bin/env python3
"""External adapter that invokes the verbatim mergeITD core."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import Align

from mergeitd_core import alignITD, annotateCoords


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-insert-seq-length", type=int, default=6)
    parser.add_argument("--min-total-reads", type=int, default=1)
    parser.add_argument("--min-vaf", type=float, default=0.006)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
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
        "COST_MATCH": 5,
        "COST_MISMATCH": -15,
        "COST_GAPOPEN": -36,
        "COST_GAPEXTEND": -0.5,
        "MIN_ALIGN_LEN": 6,
        "MIN_REF_ALN_FRACTION": 0.4,
        "MAX_FRAC_INDEL": 0.7,
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
