# Alignment and optional recounting

This describes the implemented behavior, not an assay-validated quantification
standard. Recounting is opt-in; the legacy mergeITD core is unchanged.

## Command-line interface

```bash
python itdmapper.py sample.bam --variant-alignments
python itdmapper.py sample.bam --variant-alignments=report
python itdmapper.py sample.bam --variant-alignments=lenient
python itdmapper.py sample.bam --variant-alignments=balanced
python itdmapper.py sample.bam --variant-alignments=strict
python itdmapper.py sample.fastq --recalculate-vaf --recount-mode balanced
```

The bare flag and `=report` request BAM-only TXT exports and do **not** enable
recounting. If recounting was independently enabled, they export that recount's
support. An explicit counting value enables recounting and BAM TXT exports.
`--recalculate-vaf` without a mode uses the TOML mode (`balanced` by packaged
default). `--recount-mode`
alone is an error unless recounting is enabled in TOML. Conflicting explicit
modes, or a counting-valued report flag with `--no-recalculate-vaf`, are errors.
The `=value` syntax is recommended; a bare flag before an input path remains
supported. FASTQ cannot request original-BAM-strand TXT or FASTA exports.

### Why are there two ways to select stringency?

They select **the same recount policy**, not two independent stringency levels.
The counting-valued `--variant-alignments` option is a convenience alias that
combines analysis with report generation. The separate recount options allow
analysis without TXT reports, including FASTQ analysis, which cannot export
original-BAM-strand reports.

For each `MODE` in `lenient`, `balanced`, or `strict`, these are equivalent:

```bash
python itdmapper.py sample.bam --variant-alignments=MODE
python itdmapper.py sample.bam --recalculate-vaf --recount-mode MODE --variant-alignments
```

`MODE` above denotes one of those three literal mode names, not an accepted value
itself. The following table separates analysis changes from output selection:

| Options | Enables recount? | Mode | Writes BAM TXT reports? |
|---|---|---|---|
| none | No, unless enabled in TOML | Not used without recount | No |
| `--variant-alignments` or `--variant-alignments=report` | Does not enable it | Uses independently selected recount, if enabled | Yes |
| `--recalculate-vaf` | Yes | TOML mode; packaged default balanced | No |
| `--recalculate-vaf --recount-mode strict` | Yes | strict | No |
| `--variant-alignments=strict` | Yes | strict | Yes |
| `--recalculate-vaf --recount-mode strict --variant-alignments` | Yes | strict | Yes |

In particular, `report` is an **output-only value**, not an alignment strictness
preset. It does not disable a recount enabled elsewhere. Conversely, a counting
value on the report flag changes counts/VAF as well as requesting reports.

For clarity, prefer separate analysis and output options in reproducible scripts:

```bash
# Recount only; also valid for FASTQ input.
python itdmapper.py sample.bam --recalculate-vaf --recount-mode balanced

# The same recount, with supporting-read alignment reports.
python itdmapper.py sample.bam --recalculate-vaf --recount-mode balanced --variant-alignments

# Legacy counts plus alignment reports, even if TOML enables recounting.
python itdmapper.py sample.bam --no-recalculate-vaf --variant-alignments
```

`--recalculate-vaf` is currently a **boolean flag**. It does not accept a value:
`--recalculate-vaf=balanced` is **not supported**. The suggested interface
simplification that would move mode values onto this flag has not been
implemented; no current option is deprecated.

Specifying the same mode through both interfaces is allowed. Different explicit
modes are rejected, for example `--variant-alignments=strict --recount-mode lenient`.
`--variant-alignments=strict --no-recalculate-vaf` is also rejected. CLI mode
selection overrides the TOML mode; explicit threshold overrides still apply to
either spelling, so the alias does not create a second set of thresholds.

## Input, discovery, and candidate references

Inputs are already-BBmerged single fragments. BAM input must be sorted/indexed.
Unmapped, secondary and supplementary BAM records are excluded. Configured MAPQ,
quality, duplicate, QC-fail and copy filters still apply. BAM records are fetched
around the FLT3 target; reads outside that fetch are not recoverable.

Identical sequences are collapsed for computation but weighted by original
fragment counts. BAM strand counts remain separate in exports/audits; the same
sequence/amplicon alignment computation is reused across strands. BAM SEQ,
including soft clips, is converted to transcript orientation. Hard-clipped bases
are absent and cannot be recovered. FASTQ tests both orientations. Reverse-
complement FASTQ records are not merged before the raw-sequence copy filter.

Discovery uses the existing WT-amplicon alignment and mergeITD mutation calling.
Candidates are normalized against the full FLT3 reference and keyed by HGVS/edit,
**not one per insertion length**. Different same-length alleles can compete.
All discovered target-region net insertions meeting `min_insert_length` enter
competition before final count/VAF filtering, including single-read candidates.
There is currently no top-N or candidate-support admission filter. Such artifacts
can increase ambiguity; mode selection does not solve candidate-list completeness.

WT and each compatible candidate are reconstructed within identical amplicon
boundaries. An edit must be internal to the amplicon to be reconstructed. An
external allele list and multi-ITD haplotypes are not supported. Sequence-equivalent
reconstructed alternatives are rejected rather than silently assigned duplicate
identities.

## Amplicon provenance and baseline coverage

BAM discovery assigns fragments by genomic endpoint proximity, requiring an
anchored end or sufficient overlap. FASTQ assigns sequences using WT-amplicon
alignments. Recount normally preserves those assignments for retained fragments,
even if another amplicon would score better. A narrowly defined boundary rescue
is the exception (below). Fragments rejected by that initial assignment
can compete across all inferred amplicons as a rescue; one amplicon is selected.
Tied rescue amplicons use deterministic label order, not multiple counting.

### Candidate-truncating amplicon boundaries

An inferred endpoint may be an alignment truncation rather than a true molecule
boundary. Locking every fragment to that interval can prevent testing the correct
ITD altogether. Recount therefore also tests a **strictly longer containing WT
interval** when a candidate edit intersects or lies within the selected mode's
anchor distance of the source interval and the source truncates that edit or its
required flanks. Merely overlapping or unrelated amplicons are not eligible.

All alternatives, including WT, in eligible containing amplicons participate in
score competition. Transfer occurs only when the top-scoring allele:

1. Is an ITD in a containing interval, with edit/flanks truncated by the source.
2. Passes the mode's query-participation and internal-error QC.
3. Has the required exact WT anchors on both sides of its edit.
4. Beats the best original-amplicon score by the configured margin and a positive
   difference, and beats the next distinct allele by the same rule.

If any requirement fails, the original amplicon's ordinary competition is used;
WT reads, ties and poor alignments are not transferred just to increase counts.
An internal candidate is not moved simply because a longer reference scores well.
No new reference interval or candidate allele is invented during this rescue.

Each transferred fragment has one destination assignment. Its sequence audit
retains `original_amplicon` and sets `amplicon_boundary_rescue=True`. Eligible
destination coverage has reason `amplicon_boundary_rescue`. Source candidate
coverage rows record `transferred_to_containing_amplicon` and are excluded,
preserving source baseline-loss accounting without duplicate coverage. Destination
baseline membership is not fabricated from the source span: new eligible coverage
is counted as rescued. Calls also report `boundary_rescued_support`.

The recount retains the caller's accepted WT alignment spans from `alignClasses`.
These are **not** simply BAM outer start/end coordinates. A baseline-covered
fragment has an accepted span containing the candidate's normalized insertion
boundary. Coordinates are full-reference, transcript-oriented, zero-based and
half-open. Amplicon labels are one-based inclusive genomic coordinates.

`discovery_coverage` is the original caller's coverage at its raw mutation position.
`baseline_coverage` is the retained fragment pool evaluated at the normalized
boundary in the selected recount amplicon. These may differ after normalization,
copy filtering or selection of another reporting amplicon. Both are exported to
avoid implying exact equivalence.

## Biopython alignment details

The core uses `Bio.Align.PairwiseAligner`. Effective match/mismatch and affine-gap
scores apply in all modes (packaged values: +5, -15, gap open -36, extension -0.5).
Global alignment has free terminal gaps on both sequences: it behaves as an
overlap alignment and may leave query or reference ends unexplained. A best score
alone is therefore insufficient for acceptance.

Recount places the **allele reference as target** and the read as query. It uses
score-only comparisons first, materializing near-optimal alignments for QC and
coverage. For each pair the first optimal Biopython path is inspected; the legacy
discovery/report renderer uses the last optimal path. Equivalent optimal paths
can place repeat gaps differently. Recount does not enumerate all optimal paths.

- Query participation: diagonal aligned query bases / full read length.
- Internal error fraction: (diagonal mismatches + internal gap bases) /
  (diagonal bases + internal gap bases). Terminal overhangs are not internal gaps.
- `N` is not an exact anchor match and contributes a mismatch in diagonal QC.
- Novel inserted allele bases do not project to WT coverage.
- Coverage in balanced mode requires aligned WT bases immediately before the
  edited interval and immediately after it. Mismatches at these bases are allowed.
  This is deliberately weaker than exact support anchoring.
- Confident ITD support requires exact consecutive WT matches of the configured
  anchor length on **both** sides of the normalized edit.

If the highest-scoring alignment fails QC, weaker alternatives are examined until
acceptable ones are found. They can establish coverage, but the read is not
confidently assigned to a weaker allele when the top-scoring allele failed QC.
If no alternative passes, the read is rejected in **every** mode, even if it
previously contributed coverage. Genomic placement is evidence of locus origin,
not a guarantee that the full sequence fits a closed allele list.

## Presets and overrides

| Requirement | lenient | balanced | strict |
|---|---:|---:|---:|
| Minimum query participation | 0.50 | 0.80 | 0.90 |
| Maximum internal error fraction | 0.25 | 0.15 | 0.08 |
| Exact flank bases for ITD support | 2 | 3 | 5 |
| Minimum distinct-allele score margin | 5 | 10 | 20 |
| Ambiguous site-covering reads in coverage | Yes | Yes | No |
| Preserve accepted baseline span | If any alternative passes QC | Must pass realigned site coverage | Insufficient alone |
| Coverage requires exact flank anchors | No | No | Yes |

These are provisional engineering presets, **not calibrated assay thresholds**.
All modes use the same scoring and candidate list; strictness changes QC and
counting, not a different alignment algorithm. Override individual parameters with
the existing `--reassignment-min-query-fraction`, `--reassignment-max-error-fraction`,
`--reassignment-min-anchor` and `--reassignment-min-score-margin` options.

```toml
[reassignment]
enabled = true
mode = "balanced"
# Optional explicit overrides; otherwise the mode's preset applies.
min_score_margin = 12.0
```

Precedence is CLI threshold > custom TOML threshold > selected mode preset.
Explicit mode CLI overrides custom TOML mode. Packaged defaults do not mask mode
presets. Exact score ties are ambiguous even when minimum score margin is zero.

## Distinguish assignment from coverage

An adequately aligned read can still be ambiguous. Plausible alternatives are
QC-passing alignments within the configured margin of the best QC-passing score,
including exact ties. Comparing different amplicons of the same allele does not
create distinct-allele uncertainty.

Lenient retains accepted baseline site coverage if any compatible alternative
passes QC. For rescued coverage, every plausible alignment must establish site
coverage. Balanced requires that site coverage in every plausible alignment even
for baseline fragments. Strict requires unique confident assignment and exact
flank anchors for the candidate being measured. Confident candidate support is
always included in its denominator.

Examples:

| Read | lenient/balanced coverage | Confident A support |
|---|---|---|
| Clearly A, with required anchors | Yes | Yes |
| Clearly WT and covers A site | Yes | No |
| Tied WT/B, A excluded, all plausible paths cover A | Yes | No |
| Tied WT/A, both paths cover A | Yes, unresolved A-compatible | No |
| Shared flank only, does not cover site | No (except lenient accepted baseline) | No |
| Every alternative fits poorly | No | No |

A strict denominator can shrink more than support, increasing VAF. Conversely,
including unresolved fragments in lenient/balanced coverage makes confident-support
VAF a lower assignment fraction; neither mode is automatically more accurate.

## Counting and reconciliation

Each original fragment contributes at most once per candidate denominator:

```text
coverage = baseline_coverage - lost_coverage + rescued_coverage
vaf_percent = 100 * confident_variant_support / coverage
assignment_upper_vaf_percent = 100 * (support + unresolved_support) / coverage
```

Unresolved support comprises eligible, non-confident fragments for which that
candidate remains a plausible QC-passing alternative. The upper value is an
**assignment bound, not a statistical confidence interval**. Unresolved reads
are never silently assigned WT or distributed by abundance. No confidence is
claimed about unknown alleles missing from the candidate list.

The highest-VAF compatible amplicon is still reported, not pooled. Final support
and VAF thresholds apply after recount. BAM FASTA/TXT exports contain confident
support from that selected amplicon only and agree with the tabular support count.

## Audits

With `--reassignment-audit /absolute/path/audit.tsv` (parent must exist):

- `audit.tsv`: sequence/strand/amplicon counts, mode, best/runner-up scores,
  QC metrics, assignment status/reason, original amplicon and boundary-rescue flag.
- `audit_coverage.tsv`: one row per sequence/strand/candidate/amplicon with
  baseline membership, coverage eligibility/reason, support and unresolved flags.
- `audit_calls.tsv`: original counts/VAF/coverage and recount counts, coverage,
  mode, baseline/lost/rescued coverage, unresolved support, boundary-rescued support
  and assignment bound.

Without an explicit audit path these are beside the input, using
`<sample>_reassignment` as the stem. Existing audit paths are overwritten. The
coverage audit is candidate-expanded: **do not sum it across candidates** to
count raw fragments. Each row has one primary coverage reason, not all possible
failure flags. Standard reasons include poor alignment to all alternatives,
insufficient site coverage, ambiguity (strict), and insufficient exact anchors
(strict). Accepted rows say retained, rescued, or amplicon_boundary_rescue;
source transfer rows say transferred_to_containing_amplicon.

Alignment reports render confident variant support, not all excluded or unresolved
reads. Use the coverage audit to explain denominator changes. Rendering has the
same scoring but legacy optimal-path selection and can differ in gap placement
from the first-path recount metrics.

## Validation requirements

Tests cover CLI compatibility/conflicts, mode/TOML overrides, poor alignments in
all modes, ambiguity retained in balanced coverage, unrelated-candidate ties,
flank mismatches, baseline preservation, rescue, amplicon constraints, count
conservation, and FASTQ/BAM export integration. Legacy source parity remains tested.
Real-BAM mode comparisons are recorded separately during development; their VAF
changes establish behavioral differences, not assay accuracy.

### Development comparison before boundary rescue: October 23, 2026

Both supplied BAMs were run with `config/torrent.toml`, without changing their
existing reports. Outputs and audits were written under
`/tmp/itdmapper-mode-validation` (temporary development artifacts).

| Sample / main variant | Mode | Support | Coverage | VAF (%) |
|---|---|---:|---:|---:|
| 2626702553_R_ / c.1780_1800dup | legacy | 174 | 4461 | 3.9005 |
| same | lenient | 165 | 4464 | 3.6962 |
| same | balanced | 162 | 4350 | 3.7241 |
| same | strict | 153 | 3380 | 4.5266 |
| 2627906181_R / c.1809_1810ins[72]… | legacy | 476 | 2780 | 17.1223 |
| same | lenient | 483 | 2779 | 17.3804 |
| same | balanced | 480 | 2709 | 17.7187 |
| same | strict | 470 | 2494 | 18.8452 |

For the first sample, lenient reconciles as `4461 - 2 + 5 = 4464`.
Balanced reconciles as `4461 - 114 + 3 = 4350`: 23 baseline fragments failed
alignment QC for every alternative, and 91 failed consistent realigned site
coverage. Strict loses 1081: 896 ambiguous assignments, 70 QC failures, and 115
insufficient exact-anchor fragments. Thus stricter exclusion can increase VAF.

The second sample reports five alternatives in lenient mode, two in balanced,
and one in strict (including the main insertion). The additional alternatives
started with very low discovery support and may be artifacts; this comparison
does not validate them as biological variants. Candidate admission filtering
remains a separate outstanding limitation. All sample/mode runs exited zero;
coverage reconciliation and support-plus-unresolved <= coverage were checked
for every summary candidate.

### Boundary-rescue correction: 2627906181_R

The pre-rescue results above locked 260 fragments previously supporting the main
72 bp insertion to `chr13:28608247-28608382` (136 bp WT reference). The normalized
insertion is at local boundary 136, exactly its end, so the main mutant reference
could not be constructed there. The containing interval
`chr13:28608158-28608382` (225 bp) can represent both the insertion and its flanks.
Rejecting these reads without considering that interval was not evidence that
they fit the correct allele poorly.

With the boundary-rescue rules above, fresh runs using `config/torrent.toml`
produce the following main-variant counts:

| Mode | Support before rescue | Support after rescue | Boundary-rescued support | Coverage after rescue | VAF (%) |
|---|---:|---:|---:|---:|---:|
| lenient | 483 | 746 | 263 | 3042 | 24.5233 |
| balanced | 480 | 740 | 260 | 2969 | 24.9242 |
| strict | 470 | 717 | 247 | 2741 | 26.1583 |

Artifacts are under `/tmp/itdmapper-boundary-rescue`. Existing sample reports were
not overwritten. Strict still deliberately excludes more reads through its score
margin, five-base exact anchors and 90% query-participation requirements; the
boundary fix does not relax these thresholds. Thus strict support need not equal
balanced support, nor is its higher VAF proof of greater accuracy. Candidate
admission filtering remains outstanding.

Joining sequence/strand audits against the earlier 740-support run confirms that
lenient and balanced retain all 740 previously supporting fragments. Strict
retains 717: the remaining outcomes are 4 score-ambiguous fragments, 9 with
insufficient exact anchors, and 10 with no accepted alignment in the final
competition. A failed boundary-rescue attempt falls back to the original
amplicon, so that last reason must not be read as proof that every containing
reference was biologically inappropriate. The mode's rescue criteria were not met.

Final coverage reconciles as `2780 - 1 + 263 = 3042` (lenient),
`2780 - 71 + 260 = 2969` (balanced), and `2780 - 286 + 247 = 2741` (strict).
The source interval remains explicit in the audit; identical sequences in two
original amplicons do not borrow each other's baseline membership.