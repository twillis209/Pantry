"""Pheast config handling: legacy shim, defaults merge and loaders.

Pure Python (no Snakemake import) so the unit tests can import it directly.
setup.smk loads it via workflow.source_path, which also works when the
workflow is used as a remote module.
"""

import copy
import re
from pathlib import Path

# Top-level settings shared by all datasets of a module instance.
_TOP_LEVEL_KEYS = ("fusion_dir", "gcta")

CROSS_MODALITY_MESSAGE = (
    'The "cross_modality" modality should only be used with analysis "qtl" '
    "for cross-modality cis-QTL mapping."
)


def _legacy_pathvars(config: dict) -> dict:
    """Pathvar overrides from the legacy intermediate_dir key."""
    if "intermediate_dir" not in config:
        return {}
    intermediate = str(Path(config["intermediate_dir"]).expanduser())
    return {
        "pantry_intermediate": intermediate,
        "pantry_resources": f"{intermediate}/resources",
    }


def _shim(config: dict) -> dict:
    """Wrap a single-dataset (legacy) config as dataset 'default'."""
    dataset = copy.deepcopy(config)
    if "samples_file" in dataset:
        dataset["samples"] = dataset.pop("samples_file")
    dataset.pop("intermediate_dir", None)
    out = {"datasets": {"default": dataset}, "legacy_pathvars": _legacy_pathvars(config)}
    for key in _TOP_LEVEL_KEYS:
        if key in dataset:
            out[key] = dataset.pop(key)
    return out


def _check_name(kind: str, name):
    """Reject names that cannot be used as a single path component."""
    if (
        not isinstance(name, str)
        or not re.fullmatch(r"[^/]+", name)
        or name in (".", "..")
    ):
        raise ValueError(
            f"{kind} name {name!r} must be a non-empty string without '/' "
            "and not '.' or '..'; it is used as a directory name"
        )


def _deep_merge(base: dict, override: dict) -> dict:
    """Merge override into a copy of base; nested dicts merge, other values replace."""
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _check_dataset(name: str, ds: dict):
    """Cross-field checks a JSON schema cannot express.

    These run before schema validation, so they guard the types they touch.
    Absent required keys are left to the schema.
    """
    modalities = ds.get("modalities", {})
    if not isinstance(modalities, dict):
        raise ValueError(f"datasets.{name}.modalities must be a mapping")
    analyses = ds.get("analyses", {})
    if not isinstance(analyses, dict) or any(
        not isinstance(a, dict) or not isinstance(a.get("files"), list)
        for a in analyses.values()
    ):
        raise ValueError(
            f"datasets.{name}.analyses must map analysis names to mappings "
            "with a 'files' list"
        )

    # 'cross_modality' pseudomodality should only be used for cross-modality cis-QTL mapping
    if "cross_modality" in modalities:
        if any(x != "qtl" for x in analyses) or "{modality}.trans_qtl.txt.gz" in (
            analyses.get("qtl", {}).get("files", [])
        ):
            raise ValueError(f"datasets.{name}: {CROSS_MODALITY_MESSAGE}")


def normalise(config: dict) -> dict:
    """Return a multi-dataset config built from config, which is not modified.

    A config without 'datasets' is a legacy single-dataset config and is wrapped
    as dataset 'default' (samples_file becomes samples). 'defaults' is
    deep-merged under every dataset. The result carries 'legacy_pathvars', the
    pathvar overrides implied by a legacy intermediate_dir key.
    """
    if not config:
        raise ValueError("No config file provided, or it is empty.")
    if "datasets" in config:
        out = copy.deepcopy(config)
        out["legacy_pathvars"] = {}
    else:
        out = _shim(config)
    for key in ("datasets", "defaults"):
        if not isinstance(out.get(key, {}), dict):
            raise ValueError(f"{key} must be a mapping")

    defaults = out.get("defaults", {})
    for name, ds in list(out["datasets"].items()):
        _check_name("Dataset", name)
        if not isinstance(ds, dict):
            raise ValueError(f"datasets.{name} must be a mapping")
        ds = _deep_merge(defaults, ds)
        out["datasets"][name] = ds
        _check_dataset(name, ds)
    return out


def load_samples(spec) -> list:
    """Sample IDs: a list as given, or a path to a one-column file."""
    if isinstance(spec, list):
        return spec
    import pandas as pd

    return pd.read_csv(Path(spec).expanduser(), sep="\t", header=None, dtype=str)[
        0
    ].tolist()


def qtl_filter_args(ds: dict, window: bool = True) -> str:
    """tensorQTL wrapper genotype filter arguments for a dataset.

    Matches upstream; --window is only added when qtl_cis_window is set, so
    run_tensorqtl.py's own default applies otherwise. window=False leaves it
    out for genotype_qc.py, which has no --window option.
    """
    args = f"--min_allele_count {int(ds['qtl_min_allele_count'])}"
    if ds.get("qtl_maf_threshold") is not None:
        args += f" --maf_threshold {float(ds['qtl_maf_threshold'])}"
    if window and ds.get("qtl_cis_window") is not None:
        args += f" --window {int(ds['qtl_cis_window'])}"
    return args
