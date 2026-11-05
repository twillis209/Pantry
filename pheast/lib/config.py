"""Pheast config handling: legacy shim, defaults merge and loaders.

Pure Python (no Snakemake import) so the unit tests can import it directly.
setup.smk loads it via workflow.source_path, which also works when the
workflow is used as a remote module.
"""

import copy
import re
from pathlib import Path

# FUSION override keys of module-v0.1.0/v0.1.1, removed once FUSION moved into
# the container image; rejected wherever they appear, as they would be ignored.
_REMOVED_KEYS = ("fusion_dir", "gcta", "gcta_url", "gcta_sha256")
# Top-level keys of the legacy single-dataset format. qtl_seed is new (the
# module's), but a top-level qtl_seed is still read by the shim, so it is
# rejected beside 'datasets' like the others (qtl_cis_window is not listed).
_LEGACY_TOP_LEVEL_KEYS = (
    "phenotype_dir",
    "geno_prefix",
    "samples_file",
    "modalities",
    "analyses",
    "qtl_min_allele_count",
    "qtl_maf_threshold",
    "qtl_seed",
    "twas_snps",
    "twas_models_per_job",
    "twas_threads",
    "intermediate_dir",
)

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
    return {"datasets": {"default": dataset}, "legacy_pathvars": _legacy_pathvars(config)}


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


def _reject_legacy_keys(config: dict, legacy_keys) -> None:
    """Raise if a multi-dataset config also has legacy single-dataset keys.

    They would be ignored. This happens e.g. when an old-format --configfile
    is merged over a new-format config.yml (Snakemake merges both).
    """
    found = [key for key in legacy_keys if key in config]
    if found:
        raise ValueError(
            "Config has 'datasets' and also legacy single-dataset top-level keys, "
            f"which would be ignored: {', '.join(found)}. Use one format; note "
            "that Snakemake merges --configfile with the workflow's config.yml, "
            "so an old-format config should replace config.yml (or be used from "
            "a directory without one, via --directory)."
        )


def _reject_removed_keys(config: dict) -> None:
    """Raise if the config sets a removed FUSION override key, at the top
    level, under defaults or under a dataset."""
    places = {"": config}
    if isinstance(config.get("defaults"), dict):
        places["defaults."] = config["defaults"]
    if isinstance(config.get("datasets"), dict):
        for name, ds in config["datasets"].items():
            if isinstance(ds, dict):
                places[f"datasets.{name}."] = ds
    found = [p + k for p, d in places.items() for k in _REMOVED_KEYS if k in d]
    if found:
        raise ValueError(
            f"Config keys {', '.join(found)} are no longer supported: the FUSION "
            "scripts and gcta_nr_robust now ship in the TWAS rules' container "
            "image (ghcr.io/twillis209/pantry-fusion). Remove them."
        )


def normalise(config: dict) -> dict:
    """Return a multi-dataset config built from config, which is not modified.

    A config without 'datasets' is a legacy single-dataset config and is wrapped
    as dataset 'default' (samples_file becomes samples). 'defaults' is
    deep-merged under every dataset. The result carries 'legacy_pathvars', the
    pathvar overrides implied by a legacy intermediate_dir key.
    """
    if not config:
        raise ValueError("No config file provided, or it is empty.")
    _reject_removed_keys(config)
    if "datasets" in config:
        _reject_legacy_keys(config, _LEGACY_TOP_LEVEL_KEYS)
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


def qtl_seed_arg(ds: dict) -> str:
    """run_tensorqtl.py --seed for the cis and cis_independent rules.

    Empty when qtl_seed is null (unseeded, as upstream).
    """
    if ds.get("qtl_seed") is None:
        return ""
    return f"--seed {int(ds['qtl_seed'])}"


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
