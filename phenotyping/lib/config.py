"""Phenotyping config handling: legacy shim, defaults merge and loaders.

Pure Python (no Snakemake import) so the unit tests can import it directly.
setup.smk loads it via workflow.source_path, which also works when the
workflow is used as a remote module.
"""

import copy
import math
import re
from pathlib import Path

# Legacy keys that leave the dataset: reference files and intermediate dirs.
_LEGACY_REFERENCE_KEYS = {"ref_genome": "genome", "ref_anno": "annotation"}
_LEGACY_GROUP_REFERENCE_KEYS = {
    "RNA_editing": "edit_sites_bed",
    "intron_retention": "annotation_gff3",
}
_LEGACY_DIR_KEYS = ("intermediate_dir", "intermediate_ref_dir")
# Top-level keys of the legacy single-dataset format.
_LEGACY_TOP_LEVEL_KEYS = (
    "paired_end",
    "read_length",
    "fragment_length_mean",
    "fragment_length_sd",
    "fastq_dir",
    "fastq_map",
    "ref_genome",
    "ref_anno",
    "samples_file",
    "modality_groups",
    "laddr",
) + _LEGACY_DIR_KEYS

LADDR_TARGETS = [
    "latent_residual.bed.gz",
    "latent_residual.bed.gz.tbi",
    "latent_residual.phenotype_groups.txt",
]
CROSS_MODALITY_TARGETS = [
    "cross_modality.bed.gz",
    "cross_modality.bed.gz.tbi",
    "cross_modality.phenotype_groups.txt",
]


def _legacy_pathvars(config: dict) -> dict:
    """Pathvar overrides from the legacy intermediate_dir/intermediate_ref_dir keys."""
    if not any(key in config for key in _LEGACY_DIR_KEYS):
        return {}
    intermediate = str(Path(config.get("intermediate_dir", "intermediate")).expanduser())
    reference = str(
        Path(config.get("intermediate_ref_dir", f"{intermediate}/reference")).expanduser()
    )
    return {"pantry_intermediate": intermediate, "pantry_reference": reference}


def _shim(config: dict) -> dict:
    """Wrap a single-dataset (legacy) config as reference 'default', dataset 'default'."""
    dataset = copy.deepcopy(config)
    reference = {}
    for old, new in _LEGACY_REFERENCE_KEYS.items():
        if old in dataset:
            reference[new] = dataset.pop(old)
    groups = dataset.get("modality_groups")
    if isinstance(groups, dict):
        for group, key in _LEGACY_GROUP_REFERENCE_KEYS.items():
            if isinstance(groups.get(group), dict) and key in groups[group]:
                reference[key] = groups[group].pop(key)
    if "samples_file" in dataset:
        dataset["samples"] = dataset.pop("samples_file")
    for key in _LEGACY_DIR_KEYS:
        dataset.pop(key, None)
    dataset["reference"] = "default"
    out = {
        "references": {"default": reference},
        "datasets": {"default": dataset},
        "legacy_pathvars": _legacy_pathvars(config),
    }
    # Top-level settings shared by all datasets stay top-level.
    if "majiq_env" in dataset:
        out["majiq_env"] = dataset.pop("majiq_env")
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


def _check_dataset(name: str, ds: dict, references: dict):
    """Cross-field checks a JSON schema cannot express.

    These run before schema validation, so they guard the types they touch.
    """
    if not isinstance(ds.get("reference"), str) or ds["reference"] not in references:
        raise ValueError(
            f"datasets.{name}.reference {ds.get('reference')!r} is not one of the "
            f"configured references: {', '.join(references) or '(none)'}"
        )
    reference = references[ds["reference"]]
    if not isinstance(reference, dict):
        raise ValueError(f"references.{ds['reference']} must be a mapping")

    groups = ds.get("modality_groups", {})
    if not isinstance(groups, dict) or any(
        not isinstance(g, dict) or not isinstance(g.get("files"), list)
        for g in groups.values()
    ):
        raise ValueError(
            f"datasets.{name}.modality_groups must map group names to mappings "
            "with a 'files' list"
        )

    # Fragment length parameters are only required for single-end samples.
    if ds.get("paired_end") is False:
        for key in ("fragment_length_mean", "fragment_length_sd"):
            if key not in ds:
                raise ValueError(f"datasets.{name}: {key} required for single-end data")

    if "RNA_editing" in groups and "edit_sites_bed" not in reference:
        raise ValueError(
            f"datasets.{name} uses RNA_editing, so references.{ds['reference']} "
            "needs edit_sites_bed"
        )

    cross = ds.get("cross_modality")
    if cross is not None:
        mods = cross.get("modalities") if isinstance(cross, dict) else None
        if not isinstance(mods, list) or any(not isinstance(x, str) for x in mods):
            raise ValueError(
                f"datasets.{name}.cross_modality.modalities must be a list of strings"
            )
        duplicates = sorted({x for x in mods if mods.count(x) > 1})
        if duplicates:
            raise ValueError(
                f"datasets.{name}.cross_modality.modalities lists modalities more "
                "than once: " + ", ".join(duplicates)
            )
        invalid = sorted(set(mods) - set(available_modalities(groups)))
        if invalid:
            raise ValueError(
                f"datasets.{name}.cross_modality.modalities contains modalities "
                "without a BED output: " + ", ".join(invalid)
            )

    laddr = ds.get("laddr")
    if isinstance(laddr, dict):
        if "regress_modalities" not in laddr:
            laddr["regress_modalities"] = available_modalities(groups)
        regress = laddr["regress_modalities"]
        if not isinstance(regress, list) or any(not isinstance(x, str) for x in regress):
            raise ValueError("laddr.regress_modalities must be a list of strings")
        invalid = sorted(set(regress) - set(available_modalities(groups)))
        if invalid:
            raise ValueError(
                "laddr.regress_modalities contains unknown modality names: "
                + ", ".join(invalid)
            )
        if "dataset_name" in laddr:
            _check_name("LaDDR dataset", laddr["dataset_name"])
        if "project_dir" in laddr:
            # LaDDR always works in <pantry_intermediate>/{dataset}/laddr now.
            if laddr.get("enabled"):
                print(
                    f"Note: laddr.project_dir (dataset {name!r}) is ignored; LaDDR "
                    f"uses <pantry_intermediate>/{name}/laddr."
                )
            del laddr["project_dir"]
        if laddr.get("use_existing_bins"):
            if "existing_info_dir" not in laddr or "existing_bins_dir" not in laddr:
                raise ValueError(
                    "laddr.existing_info_dir and laddr.existing_bins_dir are required "
                    "when laddr.use_existing_bins is true"
                )


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


def normalise(config: dict) -> dict:
    """Return a multi-dataset config built from config, which is not modified.

    A config without 'datasets' is a legacy single-dataset config and is wrapped
    as reference 'default' and dataset 'default'. 'defaults' is deep-merged
    under every dataset. The result carries 'legacy_pathvars', the pathvar
    overrides implied by legacy intermediate_dir/intermediate_ref_dir keys.
    """
    if not config:
        raise ValueError("No config file provided, or it is empty.")
    if "datasets" in config:
        _reject_legacy_keys(config, _LEGACY_TOP_LEVEL_KEYS)
        out = copy.deepcopy(config)
        out["legacy_pathvars"] = {}
    else:
        out = _shim(config)
    out.setdefault("references", {})
    for key in ("references", "datasets", "defaults"):
        if not isinstance(out.get(key, {}), dict):
            raise ValueError(f"{key} must be a mapping")

    for name, reference in out["references"].items():
        _check_name("Reference", name)
        if not isinstance(reference, dict):
            raise ValueError(f"references.{name} must be a mapping")
        annotation = reference.get("annotation")
        if isinstance(annotation, str) and not annotation.lower().endswith(".gtf"):
            what = "ref_anno" if "datasets" not in config else f"references.{name}.annotation"
            raise ValueError(f"{what} must be an uncompressed GTF file")

    defaults = out.get("defaults", {})
    for name, ds in list(out["datasets"].items()):
        _check_name("Dataset", name)
        if not isinstance(ds, dict):
            raise ValueError(f"datasets.{name} must be a mapping")
        ds = _deep_merge(defaults, ds)
        out["datasets"][name] = ds
        _check_dataset(name, ds, out["references"])
    return out


def load_samples(spec) -> list:
    """Sample IDs: a list as given, or a path to a one-column file."""
    if isinstance(spec, list):
        return spec
    import pandas as pd

    return pd.read_csv(Path(spec).expanduser(), sep="\t", header=None, dtype=str)[
        0
    ].tolist()


def _check_endedness(paths: dict, paired_end):
    if paired_end is not None:
        for sample_id, (_, is_paired) in paths.items():
            if is_paired != paired_end:
                raise ValueError(
                    f"Sample {sample_id} has {'paired' if is_paired else 'single'}-end "
                    f"reads, but config specifies {'paired' if paired_end else 'single'}-end"
                )


def check_fragment_length(dataset: str, ds: dict, paths: dict) -> None:
    """Require fragment_length_mean/_sd if any sample of the dataset is single-end.

    _check_dataset covers paired_end: false; this covers an omitted paired_end,
    where endedness is known only once the FASTQ map is loaded.
    """
    single = [sample_id for sample_id, (_, is_paired) in paths.items() if not is_paired]
    if not single:
        return
    missing = [k for k in ("fragment_length_mean", "fragment_length_sd") if k not in ds]
    if missing:
        raise ValueError(
            f"datasets.{dataset}: {' and '.join(missing)} required for single-end "
            f"data (sample {single[0]} has single-end reads)"
        )


def load_fastq_map(spec, fastq_dir, paired_end) -> dict:
    """Load the FASTQ file paths for each sample.

    spec is either a path to a tab-separated map file (rows 'R1 R2 sample' or
    'R sample'; paths relative to fastq_dir), or a dict
    {sample: [[R1...], [R2...]]} (paired) / {sample: [[R...]]} (single).

    Returns a dictionary where each value is a tuple (paths, is_paired).
    paths is either a list (single-end) or tuple of two lists (paired-end).
    is_paired is a boolean indicating if the sample is paired-end.
    """
    if isinstance(spec, dict):
        base = Path(fastq_dir).expanduser() if fastq_dir is not None else None

        def resolve(p: str) -> str:
            return str(base / p) if base is not None else p

        paths = {}
        for sample_id, reads in spec.items():
            if len(reads) == 2 and len(reads[0]) == len(reads[1]):
                paths[sample_id] = (
                    ([resolve(p) for p in reads[0]], [resolve(p) for p in reads[1]]),
                    True,
                )
            elif len(reads) == 1:
                paths[sample_id] = ([resolve(p) for p in reads[0]], False)
            else:
                raise ValueError(
                    f"fastq_map entry for {sample_id} must be [[R1...], [R2...]] with "
                    "equal-length lists, or [[R...]]"
                )
        _check_endedness(paths, paired_end)
        return paths

    if fastq_dir is None:
        raise ValueError("fastq_dir is required when fastq_map is a file")
    fastq_dir = Path(fastq_dir).expanduser()
    paths = {}
    with open(Path(spec).expanduser(), "r") as f:
        for line in f.read().splitlines():
            fields = line.split("\t")
            sample_id = fields[-1]
            is_paired = len(fields) == 3

            if sample_id not in paths:
                paths[sample_id] = (
                    (([], []), is_paired) if is_paired else ([], is_paired)
                )
            elif paths[sample_id][1] != is_paired:
                raise ValueError(
                    f"Sample {sample_id} has mixed paired-end and single-end reads"
                )

            if is_paired:
                paths[sample_id][0][0].append(str(fastq_dir / fields[0]))
                paths[sample_id][0][1].append(str(fastq_dir / fields[1]))
            else:
                paths[sample_id][0].append(str(fastq_dir / fields[0]))

    _check_endedness(paths, paired_end)
    return paths


def check_annotation(ref_anno):
    """Check that the GTF includes transcript entries.

    These are lacking in some annotation files such as RefSeq Select GTFs.
    """
    has_transcript = False
    with open(ref_anno, "r") as fh:
        checked = 0
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if cols[2] == "transcript":
                has_transcript = True
                break
            checked += 1
            if checked >= 1000:
                break
    if not has_transcript:
        raise ValueError(
            "Reference annotation file does not contain transcript entries "
            f"(at least not in the first 1000 lines): {ref_anno}"
        )


def edit_sites_min_samples(configured: int, fraction: float, n: int) -> int:
    """Effective RNA editing sample threshold, capped at the number of samples."""
    return min(max(configured, math.ceil(fraction * n)), n)


def available_modalities(modality_groups: dict) -> list:
    """Get phenotype modality names with normalized BED outputs."""
    modalities = []
    for params in modality_groups.values():
        for f in params["files"]:
            if f.endswith(".bed.gz"):
                modalities.append(f.removesuffix(".bed.gz"))
    return sorted(set(modalities))


def dataset_targets(ds: dict) -> list:
    """Output file names for a dataset, relative to <pantry_output>/{dataset}/."""
    targets = [f for params in ds["modality_groups"].values() for f in params["files"]]
    if ds.get("laddr", {}).get("enabled", False):
        targets.extend(LADDR_TARGETS)
    if ds.get("cross_modality"):
        targets.extend(CROSS_MODALITY_TARGETS)
    return targets
