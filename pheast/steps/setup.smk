"""This script is loaded near the beginning of the main Pheast Snakefile.

It normalises and validates the configuration and defines config helpers.
"""

import functools
import importlib.util

from snakemake.utils import validate


def _load_source(name: str, rel_path: str):
    """Import a bundled Python file by path, so it also works for a remote module."""
    spec = importlib.util.spec_from_file_location(
        name, str(workflow.source_path(rel_path))
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pantry_config = _load_source("pantry_pheast_config", "../lib/config.py")
twas_batch_ranges = _load_source(
    "pantry_pheast_twas_batches", "../scripts/twas_batches.py"
).batch_ranges


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


CONFIG = pantry_config.normalise(config)
validate(CONFIG, workflow.source_path("../schemas/config.schema.yaml"))


# The legacy intermediate_dir key overrides the Snakefile's default pathvars;
# this registers after them, so it wins.
if CONFIG["legacy_pathvars"]:

    pathvars:
        pantry_intermediate=CONFIG["legacy_pathvars"]["pantry_intermediate"],
        pantry_resources=CONFIG["legacy_pathvars"]["pantry_resources"],


def ds(w) -> dict:
    """Config of the dataset named by the {dataset} wildcard."""
    return CONFIG["datasets"][w.dataset]


@functools.cache
def dataset_samples(dataset: str) -> list:
    """Sample IDs of a dataset, read on first use."""
    return pantry_config.load_samples(CONFIG["datasets"][dataset]["samples"])


def geno_prefix(dataset: str) -> str:
    """PLINK 1 genotype prefix of a dataset."""
    return str(Path(CONFIG["datasets"][dataset]["geno_prefix"]).expanduser())


def geno_files(w) -> list:
    """Input function: the {dataset}'s PLINK .bed/.bim/.fam files."""
    return [geno_prefix(w.dataset) + ext for ext in (".bed", ".bim", ".fam")]


def geno_fam(w) -> str:
    """Input function: the {dataset}'s PLINK .fam file."""
    return geno_prefix(w.dataset) + ".fam"


def geno_input_prefix(w, input) -> str:
    """PLINK prefix of a rule's geno input, for params."""
    return input.geno[0].removesuffix(".bed")


def pheno_file(suffix: str):
    """Input function: {modality}{suffix} in the {dataset}'s phenotype_dir."""

    def path(w) -> str:
        return str(Path(ds(w)["phenotype_dir"]).expanduser() / f"{w.modality}{suffix}")

    return path


pheno_bed = pheno_file(".bed.gz")


def phenotype_summary(w) -> tuple:
    """(chroms, n_models) of the {dataset}'s {modality} BED, from the checkpoint."""
    summary = checkpoints.phenotype_summary.get(
        dataset=w.dataset, modality=w.modality
    ).output[0]
    return pantry_config.read_phenotype_summary(summary)
