"""Load and validate ITDmapper TOML settings."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Optional

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib


PACKAGED_SETTINGS_PATH = Path(__file__).resolve().parent / "config" / "default.toml"

_SCHEMA: dict[str, dict[str, tuple[type | tuple[type, ...], Optional[float], Optional[float]]]] = {
    "reads": {
        "min_mapq": (int, 0, None),
        "min_mean_base_quality": ((int, float), 0, None),
        "exclude_duplicates": (bool, None, None),
        "exclude_qcfail": (bool, None, None),
        "min_read_copies": (int, 1, None),
    },
    "amplicon": {
        "endpoint_tolerance": (int, 0, None),
        "min_cluster_fraction": ((int, float), 0, 1),
        "min_cluster_support": (int, 1, None),
        "target_fetch_flank": (int, 0, None),
        "end_anchor_tolerance": (int, 0, None),
        "min_overlap_fraction": ((int, float), 0, 1),
    },
    "fastq": {
        "discovery_top_unique": (int, 1, None),
    },
    "alignment": {
        "match": ((int, float), None, None),
        "mismatch": ((int, float), None, None),
        "gap_open": ((int, float), None, None),
        "gap_extend": ((int, float), None, None),
        "min_aligned_block": (int, 1, None),
        "min_reference_fraction": ((int, float), 0, 1),
        "min_query_fraction": ((int, float), 0, 1),
        "max_indel_fraction": ((int, float), 0, 1),
    },
    "calling": {
        "min_insert_length": (int, 0, None),
        "min_supporting_reads": (int, 1, None),
        "min_vaf_percent": ((int, float), 0, None),
    },
}


def _read_toml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"Settings file not found: {path}")
    try:
        with path.open("rb") as handle:
            data = tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"Invalid TOML in settings file {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"Settings file must contain TOML tables: {path}")
    return data


def _validate_known_keys(data: dict[str, Any], path: Path, require_all: bool) -> None:
    unknown_sections = sorted(set(data) - set(_SCHEMA))
    if unknown_sections:
        raise ValueError(
            f"Unknown settings section(s) in {path}: {', '.join(unknown_sections)}"
        )

    for section, spec in _SCHEMA.items():
        values = data.get(section)
        if values is None:
            if require_all:
                raise ValueError(f"Missing settings section [{section}] in {path}")
            continue
        if not isinstance(values, dict):
            raise ValueError(f"Settings section [{section}] must be a TOML table in {path}")

        unknown_keys = sorted(set(values) - set(spec))
        if unknown_keys:
            raise ValueError(
                f"Unknown setting(s) in [{section}] in {path}: {', '.join(unknown_keys)}"
            )
        if require_all:
            missing = sorted(set(spec) - set(values))
            if missing:
                raise ValueError(
                    f"Missing setting(s) in [{section}] in {path}: {', '.join(missing)}"
                )

        for key, value in values.items():
            expected_type, minimum, maximum = spec[key]
            if expected_type is not bool and isinstance(value, bool):
                valid_type = False
            else:
                valid_type = isinstance(value, expected_type)
            if not valid_type:
                raise ValueError(
                    f"Invalid type for [{section}].{key} in {path}: {type(value).__name__}"
                )
            if key in {
                "min_cluster_fraction",
                "min_overlap_fraction",
                "min_reference_fraction",
                "min_query_fraction",
            } and value <= 0:
                raise ValueError(f"[{section}].{key} must be > 0 in {path}")
            if minimum is not None and value < minimum:
                raise ValueError(f"[{section}].{key} must be >= {minimum} in {path}")
            if maximum is not None and value > maximum:
                raise ValueError(f"[{section}].{key} must be <= {maximum} in {path}")


def validate_settings(data: dict[str, Any], source: str = "settings") -> None:
    """Validate a complete resolved settings dictionary."""
    _validate_known_keys(data, Path(source), require_all=True)


def _merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for section, values in override.items():
        merged[section].update(values)
    return merged


def load_settings(
    custom_path: Optional[Path] = None,
    packaged_path: Path = PACKAGED_SETTINGS_PATH,
) -> dict[str, Any]:
    """Load packaged defaults, then overlay an optional partial custom TOML file."""
    packaged = _read_toml(packaged_path)
    _validate_known_keys(packaged, packaged_path, require_all=True)
    if custom_path is None:
        return packaged

    custom_path = Path(custom_path)
    custom = _read_toml(custom_path)
    _validate_known_keys(custom, custom_path, require_all=False)
    merged = _merge(packaged, custom)
    _validate_known_keys(merged, custom_path, require_all=True)
    return merged
