from math import ceil

# One MAJIQ conda env (by name) per module instance: a top-level setting.
MAJIQ_ENV = CONFIG.get("majiq_env", "majiq3")


def majiq(w) -> dict:
    """intron_retention modality group settings of the {dataset}."""
    return ds(w)["modality_groups"]["intron_retention"]


def majiq_gff3(w) -> str:
    """The reference's annotation_gff3 if given, else the GFF3 converted from its GTF."""
    reference = ref(w)
    if "annotation_gff3" in reference:
        return str(Path(reference["annotation_gff3"]).expanduser())
    return f"<pantry_reference>/{ds(w)['reference']}/majiq/annotation.gff3"


_majiq_default_license = script("intron_retention/majiq_license_academic_official.lic")


def majiq_license(w, input) -> str:
    """Configured MAJIQ license file, else the bundled academic license."""
    params = majiq(w)
    if "license" in params:
        return str(Path(params["license"]).expanduser())
    return _majiq_default_license(w, input)


def majiq_threads(w) -> int:
    return majiq(w)["threads"]


def majiq_batch_samples(w) -> list:
    """Samples of the {dataset} in psi-coverage batch {batch}."""
    size = majiq(w)["psicov_batch_size"]
    start = int(w.batch) * size
    return dataset_samples(w.dataset)[start : start + size]


def majiq_batch_bams(w) -> list:
    return expand(
        "<pantry_intermediate>/{dataset}/bam/{sample_id}.bam",
        dataset=w.dataset,
        sample_id=majiq_batch_samples(w),
    )


def majiq_batch_prefixes(w) -> str:
    return " ".join(majiq_batch_samples(w))


def majiq_psicov(w) -> list:
    """psi-coverage outputs of all batches of the {dataset}."""
    n = ceil(len(dataset_samples(w.dataset)) / majiq(w)["psicov_batch_size"])
    return expand(
        "<pantry_intermediate>/{dataset}/intron_retention/psicov/batch_{batch}.psicov.zarr",
        dataset=w.dataset,
        batch=range(n),
    )


def _majiq_sj(build_dir, bams) -> str:
    """Junction coverage paths written by majiq build to build_dir.

    MAJIQ names them after the resolved BAM, including symlink targets, so
    they are computed from the BAM paths when the job runs.
    """
    return " ".join(
        str(Path(build_dir) / f"{Path(bam).resolve().stem}.sj") for bam in bams
    )


def majiq_build_sj(w, input, output) -> str:
    """Junction coverage of all the {dataset}'s samples (majiq_build)."""
    return _majiq_sj(output.build, input.bams)


def majiq_batch_sj(w, input) -> str:
    """Junction coverage of the samples in batch {batch} (majiq_psi_coverage)."""
    return _majiq_sj(Path(input.splicegraph).parent, input.bams)


rule majiq_gff3_from_gtf:
    """Convert the reference GTF to GFF3 for MAJIQ."""
    input:
        ref_anno=ref_file("annotation"),
        anno_ok="<pantry_reference>/{reference}/annotation.ok",
    output:
        "<pantry_reference>/{reference}/majiq/annotation.gff3",
    resources:
        mem_mb=32000,
    params:
        script=script("gtf_to_majiq_gff3.py"),
        outdir=parent_dir,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.script} \
            --input {input.ref_anno} \
            --output {output}
        """


rule make_majiq_experiments_tsv:
    """Write MAJIQ experiment groups and BAM paths."""
    input:
        bams=per_sample("<pantry_intermediate>/{dataset}/bam/{sample_id}.bam"),
        bai=per_sample("<pantry_intermediate>/{dataset}/bam/{sample_id}.bam.bai"),
    output:
        "<pantry_intermediate>/{dataset}/intron_retention/experiments.tsv",
    localrule: True
    params:
        group=lambda w: majiq(w).get("group", "pantry"),
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as out:
            out.write("group\tpath\n")
            for bam in input.bams:
                out.write(f"{params.group}\t{bam}\n")


rule majiq_build:
    """Build MAJIQ splice graph and per-sample junction coverage."""
    input:
        gff3=majiq_gff3,
        experiments=rules.make_majiq_experiments_tsv.output,
        # Only used to name the junction coverage files (see _majiq_sj).
        bams=per_sample("<pantry_intermediate>/{dataset}/bam/{sample_id}.bam"),
    output:
        build=directory("<pantry_intermediate>/{dataset}/intron_retention/build"),
        splicegraph=directory(
            "<pantry_intermediate>/{dataset}/intron_retention/build/splicegraph.zarr"
        ),
    conda:
        MAJIQ_ENV
    threads: majiq_threads
    resources:
        mem_mb=64000,
        runtime="12h",
    params:
        outdir=lambda w, output: output.build,
        sj=majiq_build_sj,
        min_experiments=lambda w: majiq(w)["min_experiments"],
        license=majiq_license,
    shell:
        """
        mkdir -p {params.outdir}
        rm -rf {output.splicegraph} {params.sj}
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
        # Only used to name the junction coverage files (see _majiq_sj).
        bams=majiq_batch_bams,
    output:
        psicov=directory(
            "<pantry_intermediate>/{dataset}/intron_retention/psicov/batch_{batch}.psicov.zarr"
        ),
    conda:
        MAJIQ_ENV
    threads: majiq_threads
    resources:
        mem_mb=32000,
    params:
        sj=majiq_batch_sj,
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
            {params.sj}
        """


rule majiq_psi:
    """Estimate MAJIQ PSI values across all batches."""
    input:
        splicegraph=rules.majiq_build.output.splicegraph,
        psicov=majiq_psicov,
    output:
        "<pantry_intermediate>/{dataset}/intron_retention/majiq.psi.tsv",
    conda:
        MAJIQ_ENV
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
        "<pantry_intermediate>/{dataset}/intron_retention/majiq.psi.tsv",
    output:
        "<pantry_intermediate>/{dataset}/intron_retention/retained_intron_psi.tsv.gz",
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
        ref_anno=ref_file("annotation"),
    output:
        bed="<pantry_output>/{dataset}/unnorm/intron_retention.bed",
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
        bed="<pantry_output>/{dataset}/unnorm/intron_retention.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/intron_retention.bed.gz",
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
        "<pantry_output>/{dataset}/intron_retention.bed.gz",
    output:
        "<pantry_output>/{dataset}/intron_retention.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
