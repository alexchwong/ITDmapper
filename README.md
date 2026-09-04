# ITDmapper

Call FLT3 internal tandem duplications from a coordinate-sorted, indexed **single-read BAM containing already-BBmerged fragments**.

## Setup

```bash
python3 -m venv .env
source .env/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The `.env` virtual environment is ignored by Git. Activate it in each new shell before running ITDmapper.

## Usage

```bash
source .env/bin/activate
python itdmapper.py sample.bam
```

The default reference is `annotation/FLT3.fa` (hg19 FLT3 genomic sequence stored in FLT3/transcript orientation). The BAM must be aligned to hg19 and have a `.bai` or `.csi` index.

`itdmapper.py` finds primary reads overlapping FLT3 exons 14-15, infers dominant amplicon endpoint cluster(s) from the BAM, extracts each inferred WT amplicon from `FLT3.fa`, and constructs transcript/HGVS annotation from the mixed-case FASTA. Uppercase coding runs are treated as exons; because the CDS starts in exon 1, successive uppercase runs provide exon numbering.

BAM `SEQ` is converted from genomic-reference orientation to FLT3 transcript orientation before analysis. Reads are collapsed to unique sequences/counts and passed to `mergeitd_runner.py` in a temporary directory. The runner invokes `mergeitd_core.py`, which contains the relevant mergeITD alignment/calling functions extracted from `mergeitd.py` without algorithmic changes.

Output is HGVS on `NM_004119.3`, restricted to net insertions in exon 14, intron 14, or exon 15:

```text
HGVS    netInsert    counts    vaf_percent    coverage    amplicon
NM_004119.3:c....
```

or:

```text
NO_FLT3_ITD
```

### mergeITD defaults retained

- match: `5`
- mismatch: `-15`
- gap open: `-36`
- gap extend: `-0.5`
- minimum aligned block: `6 bp`
- minimum reference aligned fraction: `0.4`
- maximum indel fraction: `0.7`
- minimum read copies: `2`
- minimum net insertion: `6 bp`
- minimum supporting reads: `1`
- minimum VAF: `0.006%`

The last four thresholds can be adjusted from the CLI where applicable; alignment scoring is intentionally fixed to the legacy mergeITD values.

## Tests

The test suite uses Python's built-in `unittest` framework. From the repository root, run:

```bash
.env/bin/python -m unittest discover -s tests -v
```

`test_mergeitd_parity.py` checks that the extracted core functions remain source-identical (apart from trailing whitespace) to `mergeitd.py`. `test_cli.py` uses the BAM under `fixtures/` when `pysam` is installed.
