"""HGVS-oriented normalization helpers for mergeITD mutation calls.

Coordinates in this module are 0-based, half-open indices on a reference stored
in transcript 5'->3' orientation.  HGVS coordinate strings are generated only
after sequence-level 3' normalization.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class ReferenceEdit:
    """Replace reference[start:end] with alt."""

    start: int
    end: int
    alt: str

    @property
    def deleted_length(self) -> int:
        return self.end - self.start

    @property
    def inserted_length(self) -> int:
        return len(self.alt)

    @property
    def net_insert(self) -> int:
        return self.inserted_length - self.deleted_length


@dataclass(frozen=True)
class NormalizedMutation:
    edit: ReferenceEdit
    kind: str
    hgvs: str

    @property
    def site_index(self) -> int:
        return self.edit.start


_MUTATION_RE = re.compile(
    r"^(?P<start>\d+)(?:_(?P<end>\d+))?"
    r"(?P<op>delins\[(?P<del_len>\d+),(?P<delins_len>\d+)\]"
    r"|ins\[(?P<ins_len>\d+)\]|dup|del)"
    r"(?P<seq>[ACGTNacgtn]*)$"
)


def _upper(sequence: str) -> str:
    return sequence.upper()


def _validate_edit(reference: str, edit: ReferenceEdit) -> None:
    if edit.start < 0 or edit.end < edit.start or edit.end > len(reference):
        raise ValueError(
            f"Edit {edit.start}:{edit.end} lies outside reference length {len(reference)}"
        )


def apply_edit(reference: str, edit: ReferenceEdit) -> str:
    """Apply an edit and return the resulting sequence."""
    reference = _upper(reference)
    edit = ReferenceEdit(edit.start, edit.end, _upper(edit.alt))
    _validate_edit(reference, edit)
    return reference[: edit.start] + edit.alt + reference[edit.end :]


def minimize_edit(reference: str, edit: ReferenceEdit) -> ReferenceEdit:
    """Trim identical replacement prefix/suffix while preserving the haplotype."""
    reference = _upper(reference)
    edit = ReferenceEdit(edit.start, edit.end, _upper(edit.alt))
    _validate_edit(reference, edit)
    original_mutant = apply_edit(reference, edit)

    start, end, alt = edit.start, edit.end, edit.alt
    while start < end and alt and reference[start] == alt[0]:
        start += 1
        alt = alt[1:]
    while start < end and alt and reference[end - 1] == alt[-1]:
        end -= 1
        alt = alt[:-1]

    minimized = ReferenceEdit(start, end, alt)
    if apply_edit(reference, minimized) != original_mutant:
        raise AssertionError("Edit minimization changed the resulting sequence")
    return minimized


def _rightmost_equivalent_same_lengths(
    reference: str, edit: ReferenceEdit
) -> ReferenceEdit:
    """Move an edit to the right-most equivalent start with unchanged allele lengths.

    Prefix and suffix equality are scanned once, avoiding repeated construction of
    full candidate haplotypes.  The replacement sequence may rotate while moving
    through a tandem repeat.
    """
    reference = _upper(reference)
    edit = ReferenceEdit(edit.start, edit.end, _upper(edit.alt))
    _validate_edit(reference, edit)
    mutant = apply_edit(reference, edit)
    deleted_length = edit.deleted_length
    inserted_length = edit.inserted_length

    prefix_equal_end = 0
    prefix_limit = min(len(reference), len(mutant))
    while (
        prefix_equal_end < prefix_limit
        and reference[prefix_equal_end] == mutant[prefix_equal_end]
    ):
        prefix_equal_end += 1

    ref_i = len(reference) - 1
    mut_i = len(mutant) - 1
    while ref_i >= 0 and mut_i >= 0 and reference[ref_i] == mutant[mut_i]:
        ref_i -= 1
        mut_i -= 1
    suffix_ref_start = ref_i + 1

    max_start = min(prefix_equal_end, len(reference) - deleted_length)
    min_start = max(edit.start, suffix_ref_start - deleted_length)
    if max_start < min_start:
        return edit

    new_start = max_start
    new_end = new_start + deleted_length
    new_alt = mutant[new_start : new_start + inserted_length]
    shifted = ReferenceEdit(new_start, new_end, new_alt)
    if apply_edit(reference, shifted) != mutant:
        raise AssertionError("3' normalization changed the resulting sequence")
    return shifted


def normalize_3prime(reference: str, edit: ReferenceEdit) -> ReferenceEdit:
    """Return the most 3' sequence-equivalent minimal edit.

    The reference must be supplied in the coordinate system's 5'->3' orientation.
    For FLT3 this is transcript orientation, not increasing genomic coordinate.
    """
    reference = _upper(reference)
    original = ReferenceEdit(edit.start, edit.end, _upper(edit.alt))
    _validate_edit(reference, original)
    original_mutant = apply_edit(reference, original)

    current = minimize_edit(reference, original)
    for _ in range(8):
        shifted = _rightmost_equivalent_same_lengths(reference, current)
        normalized = minimize_edit(reference, shifted)
        if normalized == current:
            break
        current = normalized
    else:  # pragma: no cover - defensive guard against an implementation loop
        raise AssertionError("3' normalization did not converge")

    if apply_edit(reference, current) != original_mutant:
        raise AssertionError("Normalized edit is not sequence-equivalent to raw edit")
    return current


def classify_edit(reference: str, edit: ReferenceEdit) -> str:
    """Classify a normalized edit using HGVS prioritization."""
    reference = _upper(reference)
    edit = ReferenceEdit(edit.start, edit.end, _upper(edit.alt))
    _validate_edit(reference, edit)

    if edit.deleted_length == 0 and edit.inserted_length > 0:
        length = edit.inserted_length
        if edit.start >= length and reference[edit.start - length : edit.start] == edit.alt:
            return "dup"
        return "ins"
    if edit.deleted_length > 0 and edit.inserted_length == 0:
        return "del"
    if edit.deleted_length > 0 and edit.inserted_length > 0:
        return "delins"
    return "identity"


def _coord_range(
    start: int,
    end_exclusive: int,
    coord_for_index: Callable[[int], str],
) -> str:
    if end_exclusive <= start:
        raise ValueError("HGVS base range must contain at least one reference base")
    first = coord_for_index(start)
    last = coord_for_index(end_exclusive - 1)
    return first if end_exclusive - start == 1 else f"{first}_{last}"


def format_hgvs(
    reference: str,
    edit: ReferenceEdit,
    coord_for_index: Callable[[int], str],
) -> tuple[str, str]:
    """Format a normalized edit using mergeITD's existing sequence-length style."""
    reference = _upper(reference)
    edit = ReferenceEdit(edit.start, edit.end, _upper(edit.alt))
    kind = classify_edit(reference, edit)

    if kind == "dup":
        source_start = edit.start - edit.inserted_length
        source_end = edit.start
        return kind, f"c.{_coord_range(source_start, source_end, coord_for_index)}dup"
    if kind == "ins":
        if edit.start <= 0 or edit.start >= len(reference):
            raise ValueError(
                f"Insertion boundary {edit.start} cannot be represented by two flanking bases"
            )
        left = coord_for_index(edit.start - 1)
        right = coord_for_index(edit.start)
        return kind, f"c.{left}_{right}ins[{edit.inserted_length}]{edit.alt}"
    if kind == "del":
        return kind, f"c.{_coord_range(edit.start, edit.end, coord_for_index)}del"
    if kind == "delins":
        ref_range = _coord_range(edit.start, edit.end, coord_for_index)
        return (
            kind,
            f"c.{ref_range}delins[{edit.deleted_length},{edit.inserted_length}]{edit.alt}",
        )
    raise ValueError("Identity edit cannot be formatted as a mutation")


def parse_mergeitd_edit(
    mutation_name: str,
    amplicon_reference: str,
    full_reference_offset: int,
) -> ReferenceEdit:
    """Convert mergeITD's amplicon-local mutation name to a full-reference edit."""
    match = _MUTATION_RE.fullmatch(mutation_name)
    if not match:
        raise ValueError(f"Unsupported mergeITD mutation name: {mutation_name!r}")

    amplicon_reference = _upper(amplicon_reference)
    start = int(match.group("start"))
    end_text = match.group("end")
    end = int(end_text) if end_text is not None else start
    op = match.group("op")
    seq = _upper(match.group("seq"))

    if op.startswith("ins["):
        declared = int(match.group("ins_len"))
        if declared != len(seq):
            raise ValueError(
                f"mergeITD insertion length {declared} does not match sequence length {len(seq)}"
            )
        if not 0 <= start <= len(amplicon_reference):
            raise ValueError(f"Insertion boundary {start} lies outside inferred amplicon")
        return ReferenceEdit(full_reference_offset + start, full_reference_offset + start, seq)

    if op == "dup":
        if end < start or start < 0 or end >= len(amplicon_reference):
            raise ValueError(f"Duplication {start}_{end} lies outside inferred amplicon")
        duplicated = amplicon_reference[start : end + 1]
        boundary = end + 1
        return ReferenceEdit(
            full_reference_offset + boundary,
            full_reference_offset + boundary,
            duplicated,
        )

    if end < start or start < 0 or end >= len(amplicon_reference):
        raise ValueError(f"Reference edit {start}_{end} lies outside inferred amplicon")

    if op == "del":
        return ReferenceEdit(full_reference_offset + start, full_reference_offset + end + 1, "")

    if op.startswith("delins["):
        declared_deleted = int(match.group("del_len"))
        declared_inserted = int(match.group("delins_len"))
        observed_deleted = end - start + 1
        if declared_deleted != observed_deleted or declared_inserted != len(seq):
            raise ValueError(
                "mergeITD delins length annotation does not match mutation coordinates/sequence"
            )
        return ReferenceEdit(
            full_reference_offset + start,
            full_reference_offset + end + 1,
            seq,
        )

    raise ValueError(f"Unsupported mergeITD operation in {mutation_name!r}")


def normalize_mergeitd_name(
    mutation_name: str,
    amplicon_reference: str,
    full_reference: str,
    full_reference_offset: int,
    coord_for_index: Callable[[int], str],
) -> NormalizedMutation:
    raw = parse_mergeitd_edit(mutation_name, amplicon_reference, full_reference_offset)
    normalized = normalize_3prime(full_reference, raw)
    kind, hgvs = format_hgvs(full_reference, normalized, coord_for_index)
    if apply_edit(full_reference, raw) != apply_edit(full_reference, normalized):
        raise AssertionError("Raw and normalized mutation are not sequence-equivalent")
    return NormalizedMutation(edit=normalized, kind=kind, hgvs=hgvs)
