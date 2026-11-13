# FUSION weight fitting and profiling run in this container (under --sdm
# apptainer) rather than a conda env: containers/fusion (the FUSION R scripts
# and gcta_nr_robust of pheast/scripts/fusion_twas under /opt/fusion, R 4.5.3,
# plink 1.90, gemma 0.98.5, Python 3.12, as in the baseline env), built by
# .github/workflows/fusion-image.yml and pinned by digest. The /opt/fusion
# paths are params, not inputs: they exist only inside the container.
FUSION_CONTAINER = "docker://ghcr.io/twillis209/pantry-fusion@sha256:87335916d0395fe07cb6abe01e6cfa52838a85cd8bf538f32bfa11c17fdb4009"
FUSION_DIR = "/opt/fusion"


def twas_snps(w) -> str:
    """Input function: the {dataset}'s twas_snps file."""
    return str(Path(ds(w)["twas_snps"]).expanduser())


def twas_geno_input(w) -> list:
    """Genotypes for TWAS: subset to twas_snps if configured, else the dataset's."""
    if "twas_snps" in ds(w):
        return expand(
            "<pantry_intermediate>/{dataset}/twas/geno{ext}",
            dataset=w.dataset,
            ext=[".bed", ".bim", ".fam"],
        )
    return geno_files(w)


def twas_threads(w) -> int:
    """The {dataset}'s twas_threads setting."""
    return ds(w)["twas_threads"]


rule twas_geno:
    """Subset genotypes to SNPs used in TWAS LD reference panel."""
    input:
        geno=geno_files,
        snps=twas_snps,
    output:
        multiext("<pantry_intermediate>/{dataset}/twas/geno", ".bed", ".bim", ".fam"),
    conda:
        "../envs/plink.yaml"
    params:
        in_geno_prefix=geno_input_prefix,
        out_geno_prefix=lambda w, output: output[0].removesuffix(".bed"),
    shell:
        """
        plink \
            --bfile {params.in_geno_prefix} \
            --extract {input.snps} \
            --make-bed \
            --out {params.out_geno_prefix}
        """


rule twas_compute_weights_batch:
    """Use FUSION to compute TWAS weights from expression and genotypes.

    Outputs also include the actual weights
    ('<pantry_intermediate>/{dataset}/twas/{modality}/{gene}.wgt.RDat'), but the genes for
    which we end up with outputs are not known until the rule is run.
    """
    input:
        geno=twas_geno_input,
        bed=pheno_bed,
        covar="<pantry_intermediate>/{dataset}/covar/{modality}.covar.plink.tsv",
    output:
        "<pantry_intermediate>/{dataset}/twas/status_{modality}/{batch_start}_{batch_end}.tsv",
    container:
        FUSION_CONTAINER
    threads: twas_threads
    resources:
        runtime="8h",
    params:
        script=script("fit_twas_weights.py"),
        twas_geno_prefix=lambda w, input: input.geno[0].removesuffix(".bed"),
        twas_interm_dir=lambda w, output: str(Path(output[0]).parent.parent),
        fusion_script=f"{FUSION_DIR}/FUSION.compute_weights.R",
        gcta=f"{FUSION_DIR}/gcta_nr_robust",
    shell:
        """
        python {params.script} \
            --geno {params.twas_geno_prefix} \
            --bed {input.bed} \
            --covar {input.covar} \
            --modality {wildcards.modality} \
            --batch-start {wildcards.batch_start} \
            --batch-end {wildcards.batch_end} \
            --output-dir {params.twas_interm_dir} \
            --status {output} \
            --threads {threads} \
            --fusion-script {params.fusion_script} \
            --gcta {params.gcta}
        """


def twas_batch_status_input(wildcards):
    """Get the status files for all model-sized TWAS batches.

    None if the modality has no phenotypes: twas_assemble_summary then writes
    empty outputs (R33).
    """
    _, n_models = phenotype_summary(wildcards)
    for start, end in twas_batch_ranges(n_models, ds(wildcards)["twas_models_per_job"]):
        yield f"<pantry_intermediate>/{wildcards.dataset}/twas/status_{wildcards.modality}/{start}_{end}.tsv"


rule twas_assemble_summary:
    """Summarize weights from all batches/genes.

    Batch status tables define which phenotypes produced weights and retain
    compact diagnostics for phenotypes that FUSION intentionally skipped.
    If no phenotype has weights, FUSION.profile_wgt.R is not run: the script
    writes a header-only profile and, in .profile.err, the skip reasons (R33).
    """
    input:
        status=twas_batch_status_input,
    output:
        file_list="<pantry_intermediate>/{dataset}/twas/{modality}.list",
        profile="<pantry_intermediate>/{dataset}/twas/{modality}.profile",
        summary="<pantry_intermediate>/{dataset}/twas/{modality}.profile.err",
    container:
        FUSION_CONTAINER
    params:
        script=script("assemble_twas_weight_list.py"),
        # An absolute path, as the shell changes into the intermediate dir:
        profile_script=f"{FUSION_DIR}/utils/FUSION.profile_wgt.R",
        twas_interm_dir=lambda w, output: str(Path(output.file_list).parent),
    shell:
        """
        python {params.script} \
            --status {input.status} \
            --weights-dir {params.twas_interm_dir} \
            --output {output.file_list} \
            --profile {output.profile} \
            --profile-err {output.summary}
        if [ -s {output.file_list} ]; then
            cd {params.twas_interm_dir}
            Rscript {params.profile_script} \
                {wildcards.modality}.list \
                >{wildcards.modality}.profile \
                2>{wildcards.modality}.profile.err
        fi
        """


rule twas_pos_file:
    input:
        file_list="<pantry_intermediate>/{dataset}/twas/{modality}.list",
        bed=pheno_bed,
    output:
        pos_file="<pantry_intermediate>/{dataset}/twas/{modality}.pos",
    localrule: True
    conda:
        "../envs/plink.yaml"
    params:
        n_samples=lambda w: len(dataset_samples(w.dataset)),
    shell:
        """
        echo 'WGT\tID\tCHR\tP0\tP1\tN' >{output.pos_file}
        cut -d'/' -f2 {input.file_list} \
            | sed 's/.wgt.RDat//' \
            | paste {input.file_list} - \
            | sort -k 2b,2 \
            | join -1 2 -2 4 - <(zcat {input.bed} | cut -f1-4 | sort -k4) \
            | awk '{{OFS="\t"; print $2, $1, $3, $4, $5, {params.n_samples}}}' \
                >>{output.pos_file}
        """


rule twas_compress_output:
    """Combine TWAS weights and summary files into a single archive.

    Uses the same format as these:
    http://gusevlab.org/projects/fusion/#single-tissue-gene-expression
    """
    input:
        file_list="<pantry_intermediate>/{dataset}/twas/{modality}.list",
        profile="<pantry_intermediate>/{dataset}/twas/{modality}.profile",
        summary="<pantry_intermediate>/{dataset}/twas/{modality}.profile.err",
        pos_file="<pantry_intermediate>/{dataset}/twas/{modality}.pos",
    output:
        "<pantry_output>/{dataset}/twas/{modality}.tar.bz2",
    conda:
        "../envs/fusion.yaml"
    params:
        twas_interm_dir=lambda w, input: str(Path(input.file_list).parent),
        twas_output_dir=parent_dir,
    shell:
        """
        mkdir -p {params.twas_output_dir}
        tar -cjf {output} \
            -C {params.twas_interm_dir} \
            {wildcards.modality}.list \
            {wildcards.modality}.profile \
            {wildcards.modality}.profile.err \
            {wildcards.modality}.pos \
            --files-from={input.file_list}
        """
