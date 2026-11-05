# tensorQTL rules run in this container (under --sdm apptainer) rather than a
# conda env: containers/tensorqtl (tensorqtl 1.0.10, torch 2.3.1+cu121, pandas
# 2.3.3, Python 3.12, as in the baseline env), built by
# .github/workflows/tensorqtl-image.yml and pinned by digest (tag 1.0.10).
# GPU use needs apptainer's --nv (see tests/profile).
TENSORQTL_CONTAINER = "docker://ghcr.io/twillis209/pantry-tensorqtl@sha256:65a42f9fd74a4bc5f0f6a1a2a88a19e1140b62bbef5b3f2e6fc9de4330bde9b8"


def grouped(wildcards) -> bool:
    """Whether the {dataset}'s {modality} has a phenotype groups file."""
    return ds(wildcards)["modalities"][wildcards.modality]["grouped"]


def groups_arg(wildcards, input):
    """Pass the phenotype groups file as a Python arg if applicable"""
    if grouped(wildcards):
        return f"--groups {input.groups}"
    else:
        return ""


def groups_input(wildcards):
    """Include the phenotype groups file as an input if applicable"""
    if grouped(wildcards):
        return pheno_file(".phenotype_groups.txt")(wildcards)
    else:
        return []


def qtl_filter_args(wildcards):
    """Pass the {dataset}'s QTL genotype filter settings to tensorQTL wrapper scripts."""
    return pantry_config.qtl_filter_args(ds(wildcards))


def qtl_seed_arg(wildcards):
    """run_tensorqtl.py --seed for the {dataset}, or nothing if qtl_seed is null.

    A wildcards-only function, so its value is a plain param and a changed
    seed triggers a rerun.
    """
    return pantry_config.qtl_seed_arg(ds(wildcards))


def genotype_qc_filter_args(wildcards):
    """The same filter settings for genotype_qc.py, which takes no --window."""
    return pantry_config.qtl_filter_args(ds(wildcards), window=False)


rule qtl_genotype_qc:
    """Report genotype variants below the QTL allele-count/MAF threshold."""
    input:
        geno=geno_files,
        samples_ok="<pantry_intermediate>/{dataset}/genotype_samples.ok",
    output:
        "<pantry_output>/{dataset}/qtl/genotype_qc.tsv",
    conda:
        "../envs/python.yaml"
    params:
        script=script("genotype_qc.py"),
        geno_prefix=geno_input_prefix,
        qtl_dir=parent_dir,
        filter_args=genotype_qc_filter_args,
    shell:
        """
        mkdir -p {params.qtl_dir}
        python3 {params.script} \
            {params.geno_prefix} \
            {output} \
            {params.filter_args}
        """


rule tensorqtl_cis:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        geno=geno_files,
        bed=pheno_bed,
        bedi=pheno_file(".bed.gz.tbi"),
        covar="<pantry_intermediate>/{dataset}/covar/{modality}.covar.tsv",
        groups=groups_input,
    output:
        "<pantry_output>/{dataset}/qtl/{modality}.cis_qtl.txt.gz",
    container:
        TENSORQTL_CONTAINER
    resources:
        mem_mb=32000,
        runtime="12h",
    params:
        script=script("run_tensorqtl.py"),
        geno_prefix=geno_input_prefix,
        qtl_dir=parent_dir,
        groups_arg=groups_arg,
        filter_args=qtl_filter_args,
        seed_arg=qtl_seed_arg,
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 {params.script} \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.groups_arg} \
            {params.filter_args} \
            {params.seed_arg} \
            --mode cis
        """


rule tensorqtl_cis_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        geno=geno_files,
        bed=pheno_bed,
        bedi=pheno_file(".bed.gz.tbi"),
        covar="<pantry_intermediate>/{dataset}/covar/{modality}.covar.tsv",
        groups=groups_input,
        cis="<pantry_output>/{dataset}/qtl/{modality}.cis_qtl.txt.gz",
    output:
        "<pantry_output>/{dataset}/qtl/{modality}.cis_independent_qtl.txt.gz",
    container:
        TENSORQTL_CONTAINER
    resources:
        mem_mb=32000,
        runtime="20h",
    params:
        script=script("run_tensorqtl.py"),
        geno_prefix=geno_input_prefix,
        groups_arg=groups_arg,
        filter_args=qtl_filter_args,
        seed_arg=qtl_seed_arg,
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        python3 {params.script} \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            {params.groups_arg} \
            {params.filter_args} \
            {params.seed_arg} \
            --mode cis_independent
        """


rule tensorqtl_trans:
    """Map trans-QTL associations without permutation testing.
    Outputs all variant pairs with TSS distance > 5 Mb, configured MAF threshold, and p < 1e-5.

    """
    input:
        geno=geno_files,
        bed=pheno_bed,
        bedi=pheno_file(".bed.gz.tbi"),
        covar="<pantry_intermediate>/{dataset}/covar/{modality}.covar.tsv",
    output:
        "<pantry_output>/{dataset}/qtl/{modality}.trans_qtl.txt.gz",
    container:
        TENSORQTL_CONTAINER
    resources:
        mem_mb=32000,
        runtime="12h",
    params:
        script=script("run_tensorqtl.py"),
        geno_prefix=geno_input_prefix,
        qtl_dir=parent_dir,
        filter_args=qtl_filter_args,
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 {params.script} \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.filter_args} \
            --mode trans
        """
