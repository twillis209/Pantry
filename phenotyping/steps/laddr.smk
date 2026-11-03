def laddr(w) -> dict:
    """LaDDR settings of the {dataset}."""
    return ds(w)["laddr"]


def laddr_dataset_name(w) -> str:
    """Dataset ID used in LaDDR's manifest and outputs; defaults to the dataset key."""
    return laddr(w).get("dataset_name", w.dataset)


def laddr_setting(key: str):
    """Params function: a LaDDR setting of the {dataset}."""

    def setting(w):
        return laddr(w)[key]

    return setting


def laddr_existing_dir(key: str):
    """Params function: an existing LaDDR directory of the {dataset}."""

    def path(w) -> str:
        return str(Path(laddr(w)[key]).expanduser())

    return path


def laddr_use_existing_bins_flag(w) -> str:
    return "--use-existing-bins" if laddr(w)["use_existing_bins"] else ""


def laddr_project_dir(w, input) -> str:
    """LaDDR project directory, derived from the LaDDR config input."""
    return str(Path(input.config).parent)


def laddr_bigwigs(w) -> list:
    """bigWig coverage of the {dataset}'s samples, laid out as LaDDR expects."""
    return expand(
        "<pantry_intermediate>/{dataset}/laddr/covg_bigwig/{laddr_name}/{sample_id}.bw",
        dataset=w.dataset,
        laddr_name=laddr_dataset_name(w),
        sample_id=dataset_samples(w.dataset),
    )


def laddr_bins_ready(w) -> str:
    """Marker of linked existing bins or of bins built by LaDDR."""
    kind = "existing" if laddr(w)["use_existing_bins"] else "built"
    return f"<pantry_intermediate>/{w.dataset}/laddr/bins_ready.{kind}"


def laddr_coverage_done(w) -> str:
    return f"<pantry_intermediate>/{w.dataset}/laddr/covg_norm/{laddr_dataset_name(w)}/.coverage_done"


def laddr_pheno_inputs(w) -> list:
    return [
        f"<pantry_output>/{w.dataset}/{name}.bed.gz"
        for name in laddr(w)["regress_modalities"]
    ]


rule bam_to_bigwig:
    """Convert BAM alignments to base-resolution bigWig coverage."""
    input:
        bam="<pantry_intermediate>/{dataset}/bam/{sample_id}.bam",
        bai="<pantry_intermediate>/{dataset}/bam/{sample_id}.bam.bai",
    output:
        bw="<pantry_intermediate>/{dataset}/laddr/covg_bigwig/{laddr_name}/{sample_id}.bw",
    conda:
        "../envs/deeptools.yaml"
    threads: laddr_setting("bam_coverage_threads")
    params:
        outdir=parent_dir,
        bin_size=laddr_setting("bam_coverage_bin_size"),
    shell:
        """
        mkdir -p {params.outdir}
        bamCoverage \
            -b {input.bam} \
            -o {output.bw} \
            -of bigwig \
            --binSize {params.bin_size} \
            -p {threads}
        """


rule laddr_manifest:
    """Write coverage manifest for LaDDR."""
    input:
        laddr_bigwigs,
    output:
        manifest="<pantry_intermediate>/{dataset}/laddr/coverage_manifest.tsv",
    localrule: True
    params:
        dataset=laddr_dataset_name,
        samples=lambda w: dataset_samples(w.dataset),
    run:
        Path(str(output.manifest)).parent.mkdir(parents=True, exist_ok=True)
        with open(output.manifest, "w") as f:
            for sample_id in params.samples:
                f.write(
                    f"{params.dataset}\t{sample_id}\t{params.dataset}/{sample_id}.bw\n"
                )


rule laddr_config:
    """Write config used for LaDDR commands."""
    input:
        manifest=rules.laddr_manifest.output.manifest,
        phenos=laddr_pheno_inputs,
        ref_anno=ref_file("annotation"),
    output:
        config="<pantry_intermediate>/{dataset}/laddr/config.yaml",
    localrule: True
    conda:
        "../envs/python.yaml"
    params:
        script=script("write_laddr_config.py"),
        project_dir=parent_dir,
        manifest_rel="coverage_manifest.tsv",
        min_samples_expressed=laddr_setting("min_samples_expressed"),
        binning_method=laddr_setting("binning_method"),
        batch_size=laddr_setting("batch_size"),
        max_bin_width=laddr_setting("max_bin_width"),
        adaptive_max_samples=laddr_setting("adaptive_max_samples"),
        adaptive_bins_per_gene=laddr_setting("adaptive_bins_per_gene"),
        adaptive_min_mean_total_covg=laddr_setting("adaptive_min_mean_total_covg"),
        adaptive_max_corr=laddr_setting("adaptive_max_corr"),
        model_var_explained_max=laddr_setting("model_var_explained_max"),
        model_n_pcs_max=laddr_setting("model_n_pcs_max"),
        use_existing_bins_flag=laddr_use_existing_bins_flag,
    shell:
        """
        mkdir -p {params.project_dir}
        python3 {params.script} \
            --output {output.config} \
            --gtf {input.ref_anno} \
            --manifest {params.manifest_rel} \
            --coverage-directory covg_bigwig \
            --pheno-files {input.phenos:q} \
            {params.use_existing_bins_flag} \
            --min-samples-expressed {params.min_samples_expressed} \
            --binning-method {params.binning_method} \
            --batch-size {params.batch_size} \
            --max-bin-width {params.max_bin_width} \
            --adaptive-max-samples {params.adaptive_max_samples} \
            --adaptive-bins-per-gene {params.adaptive_bins_per_gene} \
            --adaptive-min-mean-total-covg {params.adaptive_min_mean_total_covg} \
            --adaptive-max-corr {params.adaptive_max_corr} \
            --model-var-explained-max {params.model_var_explained_max} \
            --model-n-pcs-max {params.model_n_pcs_max}
        """


rule laddr_prepare_existing_bins:
    """Attach existing LaDDR info and bins directories."""
    input:
        config=rules.laddr_config.output.config,
    output:
        ready="<pantry_intermediate>/{dataset}/laddr/bins_ready.existing",
    conda:
        "../envs/laddr.yaml"
    params:
        project_dir=laddr_project_dir,
        info_dir=laddr_existing_dir("existing_info_dir"),
        bins_dir=laddr_existing_dir("existing_bins_dir"),
    shell:
        """
        mkdir -p {params.project_dir}
        ln -sfn {params.info_dir:q} {params.project_dir}/info
        ln -sfn {params.bins_dir:q} {params.project_dir}/gene_bins
        test -s {params.project_dir}/info/n_batches.txt
        test -d {params.project_dir}/gene_bins
        touch {output.ready}
        """


rule laddr_setup:
    """Prepare gene metadata and batch layout for LaDDR."""
    input:
        config=rules.laddr_config.output.config,
        manifest=rules.laddr_manifest.output.manifest,
        bw=laddr_bigwigs,
    output:
        genes="<pantry_intermediate>/{dataset}/laddr/info/genes.tsv",
        exons="<pantry_intermediate>/{dataset}/laddr/info/exons.tsv.gz",
        n_batches="<pantry_intermediate>/{dataset}/laddr/info/n_batches.txt",
        median_covg="<pantry_intermediate>/{dataset}/laddr/info/median_coverage.txt",
    conda:
        "../envs/laddr.yaml"
    params:
        project_dir=laddr_project_dir,
    shell:
        """
        laddr setup -p {params.project_dir} -c {input.config}
        """


rule laddr_binning:
    """Generate bins for all LaDDR batches."""
    input:
        config=rules.laddr_config.output.config,
        genes=rules.laddr_setup.output.genes,
        exons=rules.laddr_setup.output.exons,
        n_batches=rules.laddr_setup.output.n_batches,
        bw=laddr_bigwigs,
    output:
        ready="<pantry_intermediate>/{dataset}/laddr/bins_ready.built",
    conda:
        "../envs/laddr.yaml"
    params:
        project_dir=laddr_project_dir,
    shell:
        """
        laddr binning -p {params.project_dir} -c {input.config}
        test -d {params.project_dir}/gene_bins
        touch {output.ready}
        """


rule laddr_coverage:
    """Generate normalized, binned coverage and regress out explicit phenotypes."""
    input:
        config=rules.laddr_config.output.config,
        manifest=rules.laddr_manifest.output.manifest,
        bins_ready=laddr_bins_ready,
        bw=laddr_bigwigs,
        phenos=laddr_pheno_inputs,
    output:
        done="<pantry_intermediate>/{dataset}/laddr/covg_norm/{laddr_name}/.coverage_done",
    conda:
        "../envs/laddr.yaml"
    params:
        project_dir=laddr_project_dir,
        dataset=laddr_dataset_name,
    shell:
        """
        laddr coverage -p {params.project_dir} -c {input.config} -d {params.dataset}
        touch {output.done}
        """


rule laddr_fit:
    """Fit latent phenotype models from normalized coverage."""
    input:
        config=rules.laddr_config.output.config,
        coverage_done=laddr_coverage_done,
    output:
        done="<pantry_intermediate>/{dataset}/laddr/models/.fit_done",
    conda:
        "../envs/laddr.yaml"
    params:
        project_dir=laddr_project_dir,
    shell:
        """
        laddr fit -p {params.project_dir} -c {input.config}
        touch {output.done}
        """


rule laddr_transform:
    """Generate latent phenotypes from fitted models.

    LaDDR names its output after the LaDDR dataset name; it is moved to a
    path that does not depend on it.
    """
    input:
        config=rules.laddr_config.output.config,
        coverage_done=laddr_coverage_done,
        fit_done=rules.laddr_fit.output.done,
    output:
        phenos="<pantry_intermediate>/{dataset}/laddr/latent_phenos.tsv.gz",
    conda:
        "../envs/laddr.yaml"
    params:
        project_dir=laddr_project_dir,
        dataset=laddr_dataset_name,
    shell:
        """
        laddr transform -p {params.project_dir} -c {input.config} -d {params.dataset}
        mv {params.project_dir}/phenotypes/latent_phenos.{params.dataset}.tsv.gz {output.phenos}
        """


rule assemble_latent_bed:
    """Convert latent phenotypes to BED format."""
    input:
        latent=rules.laddr_transform.output.phenos,
        ref_anno=ref_file("annotation"),
    output:
        bed="<pantry_output>/{dataset}/unnorm/latent_residual.bed",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    params:
        script=script("assemble_bed.py"),
        unnorm_dir=parent_dir,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.script} latent \
            --input {input.latent} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """


rule normalize_latent:
    """Quantile-normalize latent phenotypes for QTL mapping."""
    input:
        bed="<pantry_output>/{dataset}/unnorm/latent_residual.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/latent_residual.bed.gz",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    params:
        script=script("normalize_phenotypes.py"),
        bed=bgzip_source,
        outdir=parent_dir,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.script} \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """


rule pheno_groups_latent:
    """Group latent phenotypes by gene for tensorQTL."""
    input:
        "<pantry_output>/{dataset}/latent_residual.bed.gz",
    output:
        "<pantry_output>/{dataset}/latent_residual.phenotype_groups.txt",
    localrule: True
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
