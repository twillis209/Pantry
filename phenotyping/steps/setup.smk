"""This script is loaded near the beginning of the main phenotyping Snakefile.

It normalises and validates the configuration and defines config helpers.
"""

import functools
from copy import deepcopy
import importlib.util

from snakemake.utils import validate


def _load_config_lib():
    """Import lib/config.py by path, so it also works for a remote module."""
    spec = importlib.util.spec_from_file_location(
        "pantry_phenotyping_config", str(workflow.source_path("../lib/config.py"))
    )
    lib = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lib)
    return lib


pantry_config = _load_config_lib()

CONFIG = pantry_config.normalise(config)
validate(CONFIG, workflow.source_path("../schemas/config.schema.yaml"))


# Legacy intermediate_dir/intermediate_ref_dir override the Snakefile's
# default pathvars; this registers after them, so it wins.
if CONFIG["legacy_pathvars"]:

    pathvars:
        pantry_intermediate=CONFIG["legacy_pathvars"]["pantry_intermediate"],
        pantry_reference=CONFIG["legacy_pathvars"]["pantry_reference"],


def script(name: str):
    """Params callable giving the path to a bundled script.

    source_path returns a per-process cache path, so it is wrapped in a
    callable taking input: such derived params are excluded from the params
    rerun trigger.
    """
    path = workflow.source_path(f"../scripts/{name}")
    return lambda wildcards, input: path


def parent_dir(w, output) -> str:
    """Directory of a rule's first output, e.g. for mkdir -p in params."""
    return str(Path(output[0]).parent)


def bgzip_source(w, output) -> str:
    """Uncompressed path that bgzip turns into a rule's first (.gz) output."""
    return output[0].removesuffix(".gz")


def ds(w) -> dict:
    """Config of the dataset named by the {dataset} wildcard."""
    return CONFIG["datasets"][w.dataset]


def ref(w) -> dict:
    """Reference set of the {dataset} wildcard, else named by {reference}."""
    name = ds(w)["reference"] if "dataset" in w.keys() else w.reference
    return CONFIG["references"][name]


@functools.cache
def dataset_samples(dataset: str) -> list:
    """Sample IDs of a dataset, read on first use."""
    return pantry_config.load_samples(CONFIG["datasets"][dataset]["samples"])


@functools.cache
def fastq(dataset: str) -> dict:
    """FASTQ map of a dataset ({sample: (paths, is_paired)}), read on first use."""
    d = CONFIG["datasets"][dataset]
    return pantry_config.load_fastq_map(
        d["fastq_map"], d.get("fastq_dir"), d.get("paired_end")
    )


# Temporary: the old global names, taken from dataset 'default', until the
# rules take a {dataset} wildcard.
if "default" not in CONFIG["datasets"]:
    raise ValueError(
        "Until rules take a {dataset} wildcard, the phenotyping config needs a "
        "dataset named 'default' (or a single-dataset config without 'datasets')."
    )
_ds = CONFIG["datasets"]["default"]
_ref = CONFIG["references"][_ds["reference"]]

if not isinstance(_ds["samples"], str):
    raise ValueError(
        "Inline sample lists are not supported yet: until rules take a {dataset} "
        "wildcard, 'samples' must be a path to a one-column file."
    )
samples_file = Path(_ds["samples"]).expanduser()
samples = dataset_samples("default")

paired_end = _ds.get("paired_end", None)
read_length = _ds["read_length"]
fastq_map = fastq("default")
if any(not is_paired for paths, is_paired in fastq_map.values()):
    fragment_length_mean = _ds["fragment_length_mean"]
    fragment_length_sd = _ds["fragment_length_sd"]

ref_genome = Path(_ref["genome"]).expanduser()
ref_anno = Path(_ref["annotation"]).expanduser()
ref_cdna = "<pantry_reference>/cDNA.fa.gz"

modality_groups = deepcopy(_ds["modality_groups"])
# The shim moved annotation_gff3 to the reference; intron_retention.smk still
# reads it from the modality group.
if "intron_retention" in modality_groups and "annotation_gff3" in _ref:
    modality_groups["intron_retention"]["annotation_gff3"] = _ref["annotation_gff3"]
laddr_config = dict(_ds.get("laddr", {}))
for _key in ("project_dir", "existing_info_dir", "existing_bins_dir"):
    if _key in laddr_config:
        laddr_config[_key] = Path(laddr_config[_key]).expanduser()
laddr_enabled = laddr_config.get("enabled", False)

if "RNA_editing" in modality_groups:
    rna_editing_params = modality_groups["RNA_editing"]
    edit_sites_bed = Path(_ref["edit_sites_bed"]).expanduser()
    edit_sites_min_coverage = rna_editing_params["edit_sites_min_coverage"]
    edit_sites_min_samples = pantry_config.edit_sites_min_samples(
        rna_editing_params["edit_sites_min_samples"],
        rna_editing_params["edit_sites_min_samples_fraction"],
        len(samples),
    )
    edit_sites_cluster = rna_editing_params["edit_sites_cluster"]
    edit_sites_correlation_threshold = rna_editing_params[
        "edit_sites_correlation_threshold"
    ]
    if max(
        rna_editing_params["edit_sites_min_samples"],
        math.ceil(rna_editing_params["edit_sites_min_samples_fraction"] * len(samples)),
    ) > len(samples):
        print(
            f"Note: RNA editing site sample threshold is greater than the number of samples ({len(samples)}). Setting to {len(samples)}."
        )

outputs = [f"<pantry_output>/{f}" for f in pantry_config.dataset_targets(_ds)]
