"""This script is loaded near the beginning of the main Pheast Snakefile.

It normalises and validates the configuration and sets global variables.
"""

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


# Transitional: rules do not take {dataset} yet, so they use dataset 'default'.
if "default" not in CONFIG["datasets"]:
    raise ValueError(
        "Pheast rules do not yet run per dataset: the config needs a dataset "
        "named 'default' (or use the legacy single-dataset format)."
    )
DS = CONFIG["datasets"]["default"]

pheno_dir = Path(DS["phenotype_dir"]).expanduser()
geno_prefix = DS["geno_prefix"]

samples = pantry_config.load_samples(DS["samples"])
geno_samples = pd.read_csv(geno_prefix + ".fam", sep=r"\s+", header=None, dtype=str)[
    1
].tolist()
missing_samples = [s for s in samples if s not in geno_samples]
if len(missing_samples) > 0:
    raise Exception(f"Samples missing from genotypes: {missing_samples}")

modalities = DS["modalities"]
modality_chroms = {}
for modality in modalities:
    chroms = pd.read_csv(
        pheno_dir / f"{modality}.bed.gz", sep="\t", usecols=[0], dtype=str
    )
    modality_chroms[modality] = chroms.iloc[:, 0].unique().tolist()

analyses = DS["analyses"]
outputs = []
for analysis, params in analyses.items():
    if analysis == "qtl":
        outputs.append(f"<pantry_output>/{analysis}/genotype_qc.tsv")
    for f in params["files"]:
        outputs.append(
            expand(f"<pantry_output>/{analysis}/{f}", modality=modalities.keys())
        )
