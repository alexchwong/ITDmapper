# ITDmapper

Call FLT3 internal tandem duplications from already-BBmerged **single-fragment BAM or FASTQ** input.

## Setup

```bash
python3 -m venv .env
source .env/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The default reference is `annotation/FLT3.fa`, the hg19 FLT3 genomic sequence stored in FLT3/transcript orientation.

## BAM mode

The BAM must be coordinate-sorted, indexed, aligned to hg19, and contain already-BBmerged single-fragment records.

```bash
source .env/bin/activate
python itdmapper.py sample.bam
```

BAM mode finds primary reads overlapping the FLT3 exon 14-15 target, infers dominant amplicon endpoint cluster(s) from BAM genomic alignment endpoints, extracts each inferred WT amplicon from `FLT3.fa`, then realigns collapsed read sequences with mergeITD semantics for ITD/HGVS calling.

BAM `SEQ` is normalized to FLT3 transcript orientation before downstream analysis. The BAM CIGAR is not used as the final ITD alignment.

Example with explicit settings:

```bash
python itdmapper.py sample.bam \
  --settings config/default.toml \
  --min-mapq 20 \
  --min-vaf 0.1
```

## FASTQ mode

FASTQ mode expects already-BBmerged single-fragment reads. Plain and gzip-compressed FASTQ are accepted:

```bash
python itdmapper.py sample.fastq
python itdmapper.py sample.fastq.gz
```

FASTQ mode has no genomic alignment information. It therefore:

1. applies per-raw-read mean base-quality filtering;
2. collapses exact sequences and applies `min_read_copies` **before** discovery/assignment;
3. aligns the **10 most abundant unique sequences** by default against a broad FLT3 window around the exon 14-15 target;
4. tests both read orientations and clusters sequence-derived FLT3 alignment endpoints, weighted by read count, to infer one or more amplicons;
5. requires inferred amplicons to overlap the exon 14-15 target, but does **not** require their boundaries to lie inside exon 14-15;
6. aligns every retained unique sequence against every inferred amplicon;
7. rejects a sequence unless it passes the configured alignment criteria for at least one amplicon, then assigns it to the best passing amplicon and normalizes it to FLT3 transcript orientation;
8. sends the assigned sequence counts through the existing mergeITD/HGVS calling path.

Example with FASTQ-specific overrides:

```bash
python itdmapper.py sample.fastq.gz \
  --fastq-discovery-top-unique 10 \
  --min-reference-fraction 0.4 \
  --min-query-fraction 0.3 \
  --min-mean-base-quality 20
```

If no sufficiently supported FLT3 amplicon can be inferred, FASTQ mode exits with an error rather than silently falling back to a predefined amplicon.

## Input format detection

Input format is detected from `.bam`, `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`. It can be overridden explicitly:

```bash
python itdmapper.py unusual_filename --input-format fastq
python itdmapper.py unusual_filename --input-format bam
```

## Output

Output is HGVS on `NM_004119.3`, restricted to net insertions in exon 14, intron 14, or exon 15:

```text
HGVS    netInsert    counts    vaf_percent    coverage    amplicon
NM_004119.3:c....
```

or:

```text
NO_FLT3_ITD
```

Use `--verbose` to report inferred amplicons and assignment counts. FASTQ mode additionally reports raw reads, quality-passing reads, reads retained after the copy-number filter, and reads accepted/rejected by amplicon alignment.

## Settings and quality filters

All runtime thresholds are defined in packaged `config/default.toml`. A custom TOML file can contain only values that need to change. Explicit CLI options override both the custom file and packaged defaults.

### Packaged defaults

Read filters:

- minimum MAPQ: `0` (disabled; BAM only)
- minimum mean base quality: `0` (disabled; FASTQ filtering is per raw read)
- exclude duplicate-flagged reads: `false` (BAM only)
- exclude QC-fail-flagged reads: `false` (BAM only)
- minimum identical read copies: `2`

Amplicon inference/assignment:

- endpoint cluster tolerance: `5 bp`
- minimum cluster fraction relative to the dominant cluster: `0.10`
- minimum endpoint-cluster support: `2 reads`
- target/discovery flank: `1000 bp`
- BAM end-anchor tolerance: `10 bp`
- BAM minimum overlap fraction when an end is not anchored: `0.50`
- FASTQ discovery unique sequences: `10`

Alignment:

- match: `5`
- mismatch: `-15`
- gap open: `-36`
- gap extend: `-0.5`
- minimum aligned block: `6 bp`
- minimum inferred-amplicon reference fraction: `0.4`
- minimum query fraction: `0.3`
- maximum indel fraction: `0.7`

`min_query_fraction = 0.3` is chosen to remain compatible with the existing `max_indel_fraction = 0.7`: a higher query threshold such as 0.5 would become the limiting filter for very large insertions. It also rejects reads where less than 30% of the query participates in accepted alignment blocks. `min_reference_fraction = 0.4` preserves the existing mergeITD reference-coverage threshold.

Call filtering:

- minimum net insertion: `6 bp`
- minimum supporting reads: `1`
- minimum VAF: `0.006%`

### CLI overrides

```text
--input-format bam|fastq
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
--fastq-discovery-top-unique
--match-score
--mismatch-score
--gap-open-score
--gap-extend-score
--min-aligned-block
--min-reference-fraction
--min-query-fraction
--max-indel-fraction
--min-insert-seq-length
--min-total-reads
--min-vaf
```

## Tests

From the repository root:

```bash
.env/bin/python -m unittest discover -s tests -v
```

`test_fastq_mode.py` covers automatic format detection, per-read FASTQ quality filtering, pre-assignment copy filtering, multiple amplicon discovery, reverse-orientation reads, alignment acceptance/rejection, and failure when discovery candidates do not overlap the FLT3 target.

`test_mergeitd_parity.py` continues to check that the extracted mergeITD core remains source-identical to `mergeitd.py`; FASTQ support does not modify the mergeITD core.
