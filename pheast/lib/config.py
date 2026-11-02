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


def dataset_targets(ds: dict) -> list:
    """Output file names for a dataset, relative to <pantry_output>/{dataset}/."""
    targets = []
    for analysis, params in ds["analyses"].items():
        if analysis == "qtl":
            targets.append(f"{analysis}/genotype_qc.tsv")
        for f in params["files"]:
            targets.extend(
                f"{analysis}/" + f.format(modality=modality)
                for modality in ds["modalities"]
            )
    return targets


def check_genotype_samples(samples: list, fam_path) -> None:
    """Raise if any sample is absent from the PLINK .fam file (upstream message)."""
    import pandas as pd

    geno_samples = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str)[
        1
    ].tolist()
    missing_samples = [s for s in samples if s not in geno_samples]
    if len(missing_samples) > 0:
        raise Exception(f"Samples missing from genotypes: {missing_samples}")


def phenotype_summary(bed_path) -> tuple:
    """Chromosomes (unique, in file order) and phenotype count of a BED.

    The two reads are upstream's: the parse-time chromosome read for
    heritability and twas_model_count for TWAS batches.
    """
    import gzip

    import pandas as pd

    chroms = pd.read_csv(bed_path, sep="\t", usecols=[0], dtype=str)
    with gzip.open(bed_path, "rt") as handle:
        n_models = sum(1 for _ in handle) - 1
    return chroms.iloc[:, 0].unique().tolist(), n_models


def write_phenotype_summary(bed_path, out_path) -> None:
    """Write a BED's summary as a TSV with columns chroms (comma-joined), n_models."""
    chroms, n_models = phenotype_summary(bed_path)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(f"chroms\tn_models\n{','.join(chroms)}\t{n_models}\n")


def read_phenotype_summary(path) -> tuple:
    """(chroms, n_models) from a file written by write_phenotype_summary."""
    with open(path) as f:
        header, row = f.read().splitlines()
    if header != "chroms\tn_models":
        raise ValueError(f"Unexpected phenotype summary header in {path}")
    chroms, n_models = row.split("\t")
    return (chroms.split(",") if chroms else []), int(n_models)


def fusion_paths(config: dict) -> dict:
    """FUSION files used by TWAS: the install_fusion copy unless overridden.

    A configured fusion_dir must hold FUSION.compute_weights.R and
    utils/{plink_utils.R,FUSION.profile_wgt.R} (compute_weights sources
    plink_utils.R relative to itself); a configured gcta is a binary path.
    """
    installed = "<pantry_resources>/fusion"
    fusion_dir = config.get("fusion_dir")
    base = str(Path(fusion_dir).expanduser()) if fusion_dir else installed
    gcta = config.get("gcta")
    return {
        "compute_weights": f"{base}/FUSION.compute_weights.R",
        "plink_utils": f"{base}/utils/plink_utils.R",
        "profile_wgt": f"{base}/utils/FUSION.profile_wgt.R",
        "gcta": str(Path(gcta).expanduser()) if gcta else f"{installed}/gcta_nr_robust",
    }
