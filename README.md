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

## Settings and quality filters

All runtime thresholds are defined in the packaged `config/default.toml`. A custom TOML file can contain only the values that need to change:

```bash
python itdmapper.py sample.bam --settings config/default.toml
```

Explicit CLI options override both the custom file and packaged defaults. For example:

```bash
python itdmapper.py sample.bam \
  --settings config/default.toml \
  --min-mapq 20 \
  --min-vaf 0.1
```

Unknown TOML sections or keys are rejected so misspelled settings do not silently fall back to defaults.

### Packaged defaults

Read filters:

- minimum MAPQ: `0` (disabled)
- minimum mean base quality: `0` (disabled)
- exclude duplicate-flagged reads: `false`
- exclude QC-fail-flagged reads: `false`
- minimum identical read copies: `2`

Amplicon inference/assignment:

- endpoint cluster tolerance: `5 bp`
- minimum cluster fraction relative to the dominant cluster: `0.10`
- minimum endpoint-cluster support: `2 reads`
- target fetch flank: `1000 bp`
- end-anchor tolerance: `10 bp`
- minimum overlap fraction when an end is not anchored: `0.50`

mergeITD alignment:

- match: `5`
- mismatch: `-15`
- gap open: `-36`
- gap extend: `-0.5`
- minimum aligned block: `6 bp`
- minimum reference aligned fraction: `0.4`
- maximum indel fraction: `0.7`

Call filtering:

- minimum net insertion: `6 bp`
- minimum supporting reads: `1`
- minimum VAF: `0.006%`

The new MAPQ, mean-base-quality, duplicate, and QC-fail filters default to disabled, preserving the prior BAM-input behaviour. Unmapped, secondary and supplementary alignments are always excluded, and paired-end records are still rejected because the expected input is an already-BBmerged single-read BAM.

### CLI overrides

The principal settings can also be changed directly:

```text
--min-mapq
--min-mean-base-quality
--exclude-duplicates / --no-exclude-duplicates
--exclude-qcfail / --no-exclude-qcfail
--min-read-copies
--endpoint-tolerance
--min-cluster-fraction
--min-cluster-support
--target-fetch-flank
--end-anchor-tolerance
--min-overlap-fraction
--match-score
--mismatch-score
--gap-open-score
--gap-extend-score
--min-aligned-block
--min-reference-fraction
--max-indel-fraction
--min-insert-seq-length
--min-total-reads
--min-vaf
```

## Tests

The test suite uses Python's built-in `unittest` framework. From the repository root, run:

```bash
.env/bin/python -m unittest discover -s tests -v
```

`test_mergeitd_parity.py` checks that the extracted core functions remain source-identical (apart from trailing whitespace) to `mergeitd.py`. `test_cli.py` uses the BAM under `fixtures/` when `pysam` is installed. `test_settings.py` checks settings precedence, validation, and read-quality filtering.
