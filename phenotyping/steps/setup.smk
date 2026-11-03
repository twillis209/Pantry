"""This script is loaded near the beginning of the main phenotyping Snakefile.

It normalises and validates the configuration and defines config helpers.
"""

import functools
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
    paths = pantry_config.load_fastq_map(
        d["fastq_map"], d.get("fastq_dir"), d.get("paired_end")
    )
    pantry_config.check_fragment_length(dataset, d, paths)
    return paths


def ref_file(key: str):
    """Input function: a file of the reference set (genome, annotation, ...)."""

    def path(w) -> str:
        return str(Path(ref(w)[key]).expanduser())

    return path


def ref_path(name: str):
    """Input function: a generated file under the {dataset}'s reference directory."""

    def path(w) -> str:
        return f"<pantry_reference>/{ds(w)['reference']}/{name}"

    return path


def per_sample(pattern: str):
    """Input function: pattern expanded over the {dataset}'s samples."""

    def paths(w) -> list:
        return expand(pattern, dataset=w.dataset, sample_id=dataset_samples(w.dataset))

    return paths


def is_paired(w) -> bool:
    """Whether the {sample_id} of the {dataset} has paired-end reads."""
    return fastq(w.dataset)[w.sample_id][1]


def rna_editing(w) -> dict:
    """RNA_editing modality group settings of the {dataset}."""
    return ds(w)["modality_groups"]["RNA_editing"]


def rna_editing_min_samples(w) -> int:
    """Effective RNA editing sample threshold of the {dataset}."""
    params = rna_editing(w)
    n = len(dataset_samples(w.dataset))
    configured = params["edit_sites_min_samples"]
    fraction = params["edit_sites_min_samples_fraction"]
    pantry_config.note_edit_sites_threshold(w.dataset, configured, fraction, n)
    min_samples = pantry_config.edit_sites_min_samples(configured, fraction, n)
    return min_samples
