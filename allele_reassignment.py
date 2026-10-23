"""Closed-list, count-weighted allele competition using Biopython alignments.

Discovery is deliberately independent of this module. Coordinates are always
0-based, half-open and in transcript orientation; target is the allele reference.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter

from Bio.Seq import Seq

from hgvs_normalize import ReferenceEdit, apply_edit


# Query participation, internal error fraction, support anchor length, score margin.
# These are engineering presets, not assay-calibrated acceptance thresholds.
RECOUNT_PRESETS = {
    "lenient": (0.5, 0.25, 2, 5.0),
    "balanced": (0.8, 0.15, 3, 10.0),
    "strict": (0.9, 0.08, 5, 20.0),
}


@dataclass(frozen=True)
class AlleleReference:
    name: str
    amplicon: str
    sequence: str
    offset: int
    edit: ReferenceEdit | None = None

    def wt_position(self, position: int) -> int | None:
        """Map an allele base back to the full WT reference, excluding novel bases."""
        if self.edit is None:
            return self.offset + position
        start = self.edit.start - self.offset
        if position < start:
            return self.offset + position
        if position < start + self.edit.inserted_length:
            return None
        return self.offset + position - self.edit.net_insert


def build_references(amplicons, candidates):
    """Build WT and every compatible candidate for each (label, WT, offset)."""
    result = []
    for label, wt, offset in amplicons:
        result.append(AlleleReference("WT", label, wt, offset))
        seen = {wt}
        for name, edit in sorted(candidates.items()):
            if not offset < edit.start <= edit.end < offset + len(wt):
                continue
            local = ReferenceEdit(edit.start - offset, edit.end - offset, edit.alt)
            sequence = apply_edit(wt, local)
            if sequence in seen:
                raise ValueError(f"Sequence-equivalent candidate alleles in {label}: {name}")
            seen.add(sequence)
            result.append(AlleleReference(name, label, sequence, offset, edit))
    return result


@dataclass
class Match:
    reference: AlleleReference
    sequence: str
    score: float
    query_fraction: float
    error_fraction: float
    matched_positions: set[int]
    aligned_positions: set[int]

    def covers(self, edit: ReferenceEdit, anchor: int) -> bool:
        positions = {self.reference.wt_position(p) for p in self.matched_positions}
        return all(p in positions for p in range(edit.start - anchor, edit.start)) and all(
            p in positions for p in range(edit.end, edit.end + anchor)
        )

    def spans(self, edit: ReferenceEdit) -> bool:
        """Coverage permits mismatches, but requires aligned WT bases on both sides."""
        positions = {self.reference.wt_position(p) for p in self.aligned_positions}
        return edit.start - 1 in positions and edit.end in positions


def alignment_match(reference, sequence, aligner):
    """Measure query participation and internal alignment errors, not free ends."""
    alignment = next(iter(aligner.align(reference.sequence, sequence)), None)
    if alignment is None:
        return None
    coordinates = alignment.coordinates
    diagonal = []
    matched = set()
    aligned = set()
    paired = mismatches = 0
    for i in range(coordinates.shape[1] - 1):
        r0, r1 = map(int, coordinates[0, i:i + 2])
        q0, q1 = map(int, coordinates[1, i:i + 2])
        if r1 > r0 and q1 > q0:
            aligned.update(range(r0, r1))
            diagonal.append((i, q0, q1))
            paired += q1 - q0
            for r, q in zip(range(r0, r1), range(q0, q1)):
                if reference.sequence[r] == sequence[q] and sequence[q] in "ACGT":
                    matched.add(r)
                else:
                    mismatches += 1
    if not diagonal:
        return None
    gaps = 0
    for i in range(diagonal[0][0], diagonal[-1][0] + 1):
        dr = int(coordinates[0, i + 1] - coordinates[0, i])
        dq = int(coordinates[1, i + 1] - coordinates[1, i])
        if not dr or not dq:
            gaps += max(dr, dq)
    return Match(reference, sequence, float(alignment.score), paired / len(sequence),
                 (mismatches + gaps) / (paired + gaps), matched, aligned)


def reassign_sequences(counts, references, candidates, aligner, *, min_score_margin=10.0,
                       min_query_fraction=0.8, max_error_fraction=0.15,
                       min_anchor=3, both_orientations=False):
    """Return audit rows and confident assignments, aligning each sequence once.

    Each count belongs to exactly one status. Different amplicons of the same
    allele compete as a group, never as distinct runner-up alleles. A best score
    that fails acceptance is rejected, not replaced with a weaker passing allele.
    """
    audit = []
    assignments = []
    cached = {}
    for (sequence, strand), count in sorted(counts.items()):
        if sequence in cached:
            row, match = cached[sequence]
            audit.append(dict(row, strand=strand, count=count))
            if row["status"] == "assigned":
                assignments.append((match, strand, count))
            continue
        orientations = (sequence, str(Seq(sequence).reverse_complement())) if both_orientations else (sequence,)
        scored = []
        for reference in references:
            score, oriented = max((aligner.score(reference.sequence, seq), seq)
                                  for seq in orientations)
            scored.append((score, reference, oriented))
        best_by_allele = {}
        for score, reference, oriented in scored:
            previous = best_by_allele.get(reference.name)
            if previous is None or score > previous[0]:
                best_by_allele[reference.name] = (score, reference, oriented)
        ranked = sorted(best_by_allele.values(), key=lambda item: (-item[0], item[1].name))
        score, reference, oriented = ranked[0]
        runner = ranked[1][0] if len(ranked) > 1 else None
        margin = score - runner if runner is not None else None
        match = alignment_match(reference, oriented, aligner)
        reason = "assigned"
        if match is None or match.query_fraction < min_query_fraction:
            reason = "insufficient_query_coverage"
        elif match.error_fraction > max_error_fraction:
            reason = "excess_alignment_errors"
        elif margin is not None and (margin <= 1e-6 or margin < min_score_margin):
            reason = "ambiguous_score"
        elif reference.edit is not None and not match.covers(reference.edit, min_anchor):
            reason = "insufficient_variant_anchors"
        elif reference.name == "WT" and candidates and not any(
                match.covers(edit, min_anchor) for edit in candidates.values()):
            reason = "non_informative"
        status = "assigned" if reason == "assigned" else (
            "ambiguous" if reason in {"ambiguous_score", "non_informative", "insufficient_variant_anchors"}
            else "rejected")
        audit.append(dict(sequence=sequence, oriented_sequence=oriented, strand=strand, count=count,
                          status=status, reason=reason, allele=reference.name,
                          amplicon=reference.amplicon, score=score,
                          runner_up_allele=ranked[1][1].name if len(ranked) > 1 else "",
                          runner_up_score=runner, score_margin=margin,
                          query_fraction=match.query_fraction if match else 0,
                          error_fraction=match.error_fraction if match else 1))
        if status == "assigned":
            assignments.append((match, strand, count))
        cached[sequence] = (audit[-1], match)
    return audit, assignments


def quantify(candidates, assignments, min_anchor):
    """Candidate-specific denominators: confident fragments spanning both WT flanks."""
    support = Counter()
    coverage = Counter()
    for match, strand, count in assignments:
        if match.reference.name != "WT":
            support[match.reference.name] += count
        for name, edit in candidates.items():
            if match.covers(edit, min_anchor):
                coverage[name] += count
    return support, coverage


def recount_sequences(buckets, references, candidates, aligner, baseline, *, mode="balanced",
                      min_query_fraction=None, max_error_fraction=None, min_anchor=None,
                      min_score_margin=None, both_orientations=False):
    """Separate allele confidence from coverage; return sequence and coverage audits.

    Buckets are (sequence, strand, amplicon) -> count. An empty amplicon permits
    rescue against all references. Baseline maps (amplicon, sequence) to accepted
    full-WT alignment spans. Assigned reads can move only through a confident
    boundary rescue into a containing WT interval, retaining original provenance.
    """
    q, error, anchor, margin = RECOUNT_PRESETS[mode]
    q = q if min_query_fraction is None else min_query_fraction
    error = error if max_error_fraction is None else max_error_fraction
    anchor = anchor if min_anchor is None else min_anchor
    margin = margin if min_score_margin is None else min_score_margin
    audit, coverage_audit, assignments = [], [], []
    cache = {}
    wild_types = {r.amplicon: r for r in references if r.name == "WT"}
    for (sequence, strand, label), count in sorted(buckets.items()):
        cache_key = (sequence, label)
        if cache_key not in cache:
            orientations = {sequence}
            if both_orientations:
                orientations.add(str(Seq(sequence).reverse_complement()))
            rescue_labels = set()
            original_wt = wild_types.get(label)
            if original_wt:
                start, end = original_wt.offset, original_wt.offset + len(original_wt.sequence)
                for ref in references:
                    edit = ref.edit
                    if edit is None or ref.amplicon == label:
                        continue
                    containing = wild_types[ref.amplicon]
                    contains = containing.offset <= start and containing.offset + len(containing.sequence) >= end
                    truncated = (start - anchor <= edit.start <= end + anchor and
                                 edit.end >= start - anchor and
                                 (edit.start - anchor < start or edit.end + anchor > end))
                    if contains and truncated and len(containing.sequence) > len(original_wt.sequence):
                        rescue_labels.add(ref.amplicon)
            scored = []
            for ref in references:
                if label and ref.amplicon != label and ref.amplicon not in rescue_labels:
                    continue
                score, oriented = max((aligner.score(ref.sequence, seq), seq) for seq in orientations)
                scored.append((score, ref, oriented))
            if not scored:
                raise ValueError(f"No compatible references for amplicon {label!r}")
            scored.sort(key=lambda item: (-item[0], item[1].name, item[1].amplicon))
            boundary_rescue = False
            if label and rescue_labels:
                original_scores = [item for item in scored if item[1].amplicon == label]
                top_score, top_ref, top_sequence = scored[0]
                if top_ref.amplicon != label and top_ref.edit is not None:
                    runner = next((item[0] for item in scored if item[1].name != top_ref.name), None)
                    match = alignment_match(top_ref, top_sequence, aligner)
                    source_truncates = (top_ref.edit.start - anchor < original_wt.offset or
                                        top_ref.edit.end + anchor > original_wt.offset + len(original_wt.sequence))
                    boundary_rescue = bool(source_truncates and match and
                        match.query_fraction >= q and match.error_fraction <= error and
                        match.covers(top_ref.edit, anchor) and
                        top_score - original_scores[0][0] >= margin and
                        top_score - original_scores[0][0] > 1e-6 and
                        (runner is None or (top_score - runner >= margin and top_score - runner > 1e-6)))
                if not boundary_rescue:
                    scored = original_scores
            # Materialize near-optimal alignments and weaker acceptable alternatives
            # when the top score fails quality. Never force a failed alignment.
            matches = []
            for score, ref, oriented in scored:
                if matches and score < matches[0].score - max(margin, 1e-6):
                    break
                match = alignment_match(ref, oriented, aligner)
                if match and match.query_fraction >= q and match.error_fraction <= error:
                    matches.append(match)
            cache[cache_key] = (scored, matches, boundary_rescue)
        scored, matches, boundary_rescue = cache[cache_key]
        selected_label = (scored[0][1].amplicon if boundary_rescue else label) or (
            matches[0].reference.amplicon if matches else scored[0][1].amplicon)
        # A rescue may choose one amplicon only; it must not multiply coverage.
        plausible = [m for m in matches if m.reference.amplicon == selected_label]
        best = plausible[0] if plausible else None
        best_score, best_ref, oriented = scored[0]
        runner = next((item for item in scored if item[1].name != best_ref.name), None)
        delta = best_score - runner[0] if runner else None
        reason = "assigned"
        if best is None:
            reason = "poor_alignment_all_alternatives"
        elif best.reference.name != best_ref.name:
            # A weaker valid alignment can establish coverage, not confident support.
            reason = "best_alignment_failed_quality"
        elif delta is not None and (delta <= 1e-6 or delta < margin):
            reason = "ambiguous_score"
        elif best.reference.edit and not best.covers(best.reference.edit, anchor):
            reason = "insufficient_variant_anchors"
        status = "assigned" if reason == "assigned" else ("rejected" if best is None else "ambiguous")
        if status == "assigned":
            assignments.append((best, strand, count))
        audit.append(dict(sequence=sequence, oriented_sequence=best.sequence if best else oriented,
                          strand=strand, count=count, status=status, reason=reason,
                          allele=best_ref.name, amplicon=selected_label, score=best_score,
                          runner_up_allele=runner[1].name if runner else "",
                          runner_up_score=runner[0] if runner else None, score_margin=delta,
                          query_fraction=best.query_fraction if best else 0,
                          error_fraction=best.error_fraction if best else 1, mode=mode,
                          original_amplicon=label, amplicon_boundary_rescue=boundary_rescue))
        original = baseline.get((selected_label, sequence)) if selected_label == label else None
        if original is None and both_orientations and selected_label == label:
            original = baseline.get((selected_label, str(Seq(sequence).reverse_complement())))
        for name, edit in candidates.items():
            if not any(r.name == name and r.amplicon == selected_label for r in references):
                continue
            # Baseline uses the same normalized boundary for a read-by-read
            # reconciliation. Raw legacy coverage can differ after 3' shifting.
            was_covered = bool(original and original[0] <= edit.start < original[1])
            eligible = False
            coverage_reason = "insufficient_site_coverage"
            if not plausible:
                coverage_reason = "poor_alignment_all_alternatives"
            elif mode == "lenient":
                eligible = was_covered or all(m.spans(edit) for m in plausible)
            elif mode == "balanced":
                eligible = all(m.spans(edit) for m in plausible)
            else:
                eligible = status == "assigned" and best.covers(edit, anchor)
                if status != "assigned":
                    coverage_reason = "ambiguous_assignment"
                elif not eligible:
                    coverage_reason = "insufficient_exact_anchors"
            supporting = status == "assigned" and best.reference.name == name
            # Confident support must always be represented in its denominator.
            if supporting:
                eligible = True
            if eligible:
                coverage_reason = "retained" if was_covered else "rescued"
            unresolved = eligible and status != "assigned" and any(m.reference.name == name for m in plausible)
            coverage_audit.append(dict(sequence=sequence, strand=strand, count=count, amplicon=selected_label,
                HGVS=name, mode=mode, baseline_covered=was_covered, coverage_eligible=eligible,
                coverage_reason="amplicon_boundary_rescue" if boundary_rescue and eligible else coverage_reason,
                supporting=supporting, unresolved=unresolved, original_amplicon=label,
                amplicon_boundary_rescue=boundary_rescue))
        if boundary_rescue:
            # Preserve source coverage accounting without counting the moved
            # fragment in two denominators for the same candidate.
            source_span = baseline.get((label, sequence))
            if source_span is None and both_orientations:
                source_span = baseline.get((label, str(Seq(sequence).reverse_complement())))
            for name, edit in candidates.items():
                if not any(r.name == name and r.amplicon == label for r in references):
                    continue
                source_covered = bool(source_span and source_span[0] <= edit.start < source_span[1])
                coverage_audit.append(dict(sequence=sequence, strand=strand, count=count, amplicon=label,
                    HGVS=name, mode=mode, baseline_covered=source_covered, coverage_eligible=False,
                    coverage_reason="transferred_to_containing_amplicon", supporting=False, unresolved=False,
                    original_amplicon=label, amplicon_boundary_rescue=True))
    return audit, coverage_audit, assignments