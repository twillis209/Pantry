from math import ceil

majiq_dir = "<pantry_intermediate>/default/intron_retention"
majiq_build_dir = f"{majiq_dir}/build"
majiq_psicov_dir = f"{majiq_dir}/psicov"
majiq_params = modality_groups["intron_retention"]
majiq_annotation_gff3 = (
    Path(majiq_params["annotation_gff3"]).expanduser()
    if "annotation_gff3" in majiq_params
    else f"<pantry_reference>/{default_reference}/majiq/annotation.gff3"
)
majiq_license = (
    str(Path(majiq_params["license"]).expanduser())
    if "license" in majiq_params
    else script("intron_retention/majiq_license_academic_official.lic")
)
majiq_threads = majiq_params.get("threads", 8)
majiq_psicov_batch_size = majiq_params.get("psicov_batch_size", 25)
majiq_min_experiments = majiq_params.get("min_experiments", 0.01)
majiq_batches = list(range(ceil(len(samples) / majiq_psicov_batch_size)))
# MAJIQ names junction coverage after the resolved BAM, including symlink targets.
# Resolving needs the real path, so the pathvar is applied here at parse time.
majiq_sj_by_sample = {
    sample: f"{majiq_build_dir}/"
    + Path(workflow.pathvars.apply(f"<pantry_intermediate>/default/bam/{sample}.bam"))
    .resolve()
    .stem
    + ".sj"
    for sample in samples
}
majiq_sj = [directory(majiq_sj_by_sample[sample]) for sample in samples]


def majiq_batch_samples(batch):
    start = int(batch) * majiq_psicov_batch_size
    return samples[start : start + majiq_psicov_batch_size]


def majiq_batch_sj(wildcards):
    return [
        str(majiq_sj_by_sample[sample])
        for sample in majiq_batch_samples(wildcards.batch)
    ]


def majiq_batch_prefixes(wildcards):
    return " ".join(majiq_batch_samples(wildcards.batch))


if "annotation_gff3" not in majiq_params:

    rule majiq_gff3_from_gtf:
        """Convert the reference GTF to GFF3 for MAJIQ."""
        input:
            ref_anno,
        output:
            majiq_annotation_gff3,
        resources:
            mem_mb=32000,
        params:
            script=script("gtf_to_majiq_gff3.py"),
            outdir=parent_dir,
        shell:
            """
            mkdir -p {params.outdir}
            python3 {params.script} \
                --input {input} \
                --output {output}
            """


rule make_majiq_experiments_tsv:
    """Write MAJIQ experiment groups and BAM paths."""
    input:
        bams=expand(
            "<pantry_intermediate>/default/bam/{sample_id}.bam", sample_id=samples
        ),
        bai=expand(
            "<pantry_intermediate>/default/bam/{sample_id}.bam.bai", sample_id=samples
        ),
    output:
        f"{majiq_dir}/experiments.tsv",
    localrule: True
    params:
        group=majiq_params.get("group", "pantry"),
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as out:
            out.write("group\tpath\n")
            for bam in input.bams:
                out.write(f"{params.group}\t{bam}\n")


rule majiq_build:
    """Build MAJIQ splice graph and per-sample junction coverage."""
    input:
        gff3=majiq_annotation_gff3,
        experiments=rules.make_majiq_experiments_tsv.output,
    output:
        splicegraph=directory(f"{majiq_build_dir}/splicegraph.zarr"),
        sj=majiq_sj,
    conda:
        "majiq3"
    threads: majiq_threads
    resources:
        mem_mb=64000,
        runtime="12h",
    params:
        outdir=lambda w, output: str(Path(output.splicegraph).parent),
        min_experiments=majiq_min_experiments,
        license=majiq_license,
    shell:
        """
        mkdir -p {params.outdir}
        rm -rf {output.splicegraph} {output.sj}
        majiq build \
            --license {params.license} \
            --all-introns \
            --min-experiments {params.min_experiments} \
            --nthreads {threads} \
            --overwrite \
            {input.gff3} \
            {input.experiments} \
            {params.outdir}
        """


rule majiq_psi_coverage:
    """Calculate MAJIQ PSI coverage in sample batches."""
    input:
        splicegraph=rules.majiq_build.output.splicegraph,
        sj=majiq_batch_sj,
    output:
        psicov=directory(f"{majiq_psicov_dir}/batch_{{batch}}.psicov.zarr"),
    conda:
        "majiq3"
    threads: majiq_threads
    resources:
        mem_mb=32000,
    params:
        prefixes=majiq_batch_prefixes,
        license=majiq_license,
        outdir=parent_dir,
    shell:
        """
        mkdir -p {params.outdir}
        rm -rf {output.psicov}
        majiq psi-coverage \
            --license {params.license} \
            --nthreads {threads} \
            --prefixes {params.prefixes} \
            --overwrite \
            {input.splicegraph} \
            {output.psicov} \
            {input.sj}
        """


rule majiq_psi:
    """Estimate MAJIQ PSI values across all batches."""
    input:
        splicegraph=rules.majiq_build.output.splicegraph,
        psicov=expand(
            f"{majiq_psicov_dir}/batch_{{batch}}.psicov.zarr", batch=majiq_batches
        ),
    output:
        f"{majiq_dir}/majiq.psi.tsv",
    conda:
        "majiq3"
    threads: majiq_threads
    resources:
        mem_mb=64000,
    params:
        license=majiq_license,
    shell:
        """
        majiq psi \
            --license {params.license} \
            --splicegraph {input.splicegraph} \
            --output-tsv {output} \
            --quantiles 0.025 0.5 0.975 \
            --nthreads {threads} \
            --overwrite \
            {input.psicov}
        """


rule extract_retained_introns:
    """Extract retained-intron rows from MAJIQ PSI output."""
    input:
        f"{majiq_dir}/majiq.psi.tsv",
    output:
        f"{majiq_dir}/retained_intron_psi.tsv.gz",
    resources:
        mem_mb=32000,
    params:
        script=script("intron_retention/extract_ir_psi.py"),
    shell:
        """
        python3 {params.script} \
            --input {input} \
            --output {output}
        """


rule assemble_intron_retention_bed:
    """Convert MAJIQ retained-intron PSI values into BED format."""
    input:
        ir=rules.extract_retained_introns.output,
        ref_anno=ref_anno,
    output:
        bed="<pantry_output>/default/unnorm/intron_retention.bed",
    resources:
        mem_mb=32000,
    params:
        script=script("assemble_bed.py"),
        unnorm_dir=parent_dir,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.script} intron-retention \
            --input {input.ir} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """


rule normalize_intron_retention:
    """Quantile-normalize retained-intron PSI values for QTL mapping."""
    input:
        bed="<pantry_output>/default/unnorm/intron_retention.bed",
        samples=samples_file,
    output:
        "<pantry_output>/default/intron_retention.bed.gz",
    resources:
        mem_mb=32000,
    params:
        script=script("normalize_phenotypes.py"),
        bed=bgzip_source,
    shell:
        """
        python3 {params.script} \
            --input {input.bed} \
            --samples {input.samples} \
            --output {params.bed}
        bgzip {params.bed}
        """


rule pheno_groups_intron_retention:
    """Group retained-intron phenotypes by gene for tensorQTL."""
    input:
        "<pantry_output>/default/intron_retention.bed.gz",
    output:
        "<pantry_output>/default/intron_retention.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
