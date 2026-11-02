"""This script is loaded near the beginning of the main Pheast Snakefile.

It validates and processes the configuration and sets global variables.
"""

import importlib.util


def _load_source(name: str, rel_path: str):
    """Import a bundled Python file by path, so it also works for a remote module."""
    spec = importlib.util.spec_from_file_location(
        name, str(workflow.source_path(rel_path))
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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


def validate_config(config: dict):
    """Validate the configuration"""
    if len(config.keys()) == 0:
        raise Exception("No config file provided.")
    fields = [
        "phenotype_dir",
        "geno_prefix",
        "samples_file",
        "modalities",
        "analyses",
    ]
    for field in fields:
        if field not in config.keys():
            raise Exception(f"{field} not in config file.")

    # 'cross_modality' pseudomodality should only be used for cross-modality cis-QTL mapping
    if "cross_modality" in config["modalities"]:
        if (
            any([x != "qtl" for x in config["analyses"]])
            or "{modality}.trans_qtl.txt.gz" in config["analyses"]["qtl"]["files"]
        ):
            raise ValueError(
                'The "cross_modality" modality should only be used with analysis "qtl" for cross-modality cis-QTL mapping.'
            )


def process_config(config: dict):
    """Prepare user config for use in the pipeline.

    Updates some values in place to expand paths and parse numbers.
    """
    paths = [
        "phenotype_dir",
        "samples_file",
        "twas_snps",
        "intermediate_dir",
    ]
    for path in paths:
        if path in config.keys():
            config[path] = str(Path(config[path]).expanduser())

    config["samples"] = pd.read_csv(
        config["samples_file"], sep="\t", header=None, dtype=str
    )[0].tolist()

    if "qtl_min_allele_count" not in config:
        config["qtl_min_allele_count"] = 10
    config["qtl_min_allele_count"] = int(config["qtl_min_allele_count"])

    if "qtl_maf_threshold" in config and config["qtl_maf_threshold"] is not None:
        config["qtl_maf_threshold"] = float(config["qtl_maf_threshold"])

    if "twas_models_per_job" not in config:
        config["twas_models_per_job"] = 200
    config["twas_models_per_job"] = int(config["twas_models_per_job"])
    if config["twas_models_per_job"] < 1:
        raise ValueError("twas_models_per_job must be at least 1")

    if "twas_threads" not in config:
        config["twas_threads"] = 4
    config["twas_threads"] = int(config["twas_threads"])
    if config["twas_threads"] < 1:
        raise ValueError("twas_threads must be at least 1")


validate_config(config)
process_config(config)


# The legacy intermediate_dir key overrides the Snakefile's default pathvars;
# this registers after them, so it wins.
if "intermediate_dir" in config:

    pathvars:
        pantry_intermediate=config["intermediate_dir"],
        pantry_resources=f"{config['intermediate_dir']}/resources",


pheno_dir = Path(config["phenotype_dir"])
geno_prefix = config["geno_prefix"]

samples = config["samples"]
geno_samples = pd.read_csv(geno_prefix + ".fam", sep=r"\s+", header=None, dtype=str)[
    1
].tolist()
missing_samples = [s for s in samples if s not in geno_samples]
if len(missing_samples) > 0:
    raise Exception(f"Samples missing from genotypes: {missing_samples}")

modalities = config["modalities"]
for modality, params in modalities.items():
    chroms = pd.read_csv(
        pheno_dir / f"{modality}.bed.gz", sep="\t", usecols=[0], dtype=str
    )
    params["chroms"] = chroms.iloc[:, 0].unique().tolist()

analyses = config["analyses"]
outputs = []
for analysis, params in analyses.items():
    if analysis == "qtl":
        outputs.append(f"<pantry_output>/{analysis}/genotype_qc.tsv")
    for f in params["files"]:
        outputs.append(
            expand(f"<pantry_output>/{analysis}/{f}", modality=modalities.keys())
        )

qtl_min_allele_count = config["qtl_min_allele_count"]
qtl_maf_threshold = config.get("qtl_maf_threshold")
