TXREVISE_N_BATCHES = 100


rule preprocess_gtf_for_txrevise:
    """Preprocess the GTF file for txrevise"""
    input:
        ref_anno,
    output:
        "<pantry_reference>/txrevise/gene_annotations.gtf",
    resources:
        mem_mb=32000,
    params:
        script=script("txrevise/preprocess_gtf.py"),
    shell:
        """
        python3 {params.script} \
            --input {input} \
            --output {output}
        """


rule txrevise_extract_tags:
    """Extract transcript tags from the GTF file

    txrevise rules are adapted from https://github.com/kauralasoo/txrevise/blob/master/scripts/Snakefile

    """
    input:
        "<pantry_reference>/txrevise/gene_annotations.gtf",
    output:
        tags="<pantry_reference>/txrevise/transcript_tags.txt",
    params:
        script=script("txrevise/extractTranscriptTags.py"),
    shell:
        """
        python3 {params.script} \
            --gtf <(gzip -c {input}) \
            >{output.tags}
        """


rule txrevise_prepare_annotations:
    """Prepare annotations"""
    input:
        gtf="<pantry_reference>/txrevise/gene_annotations.gtf",
        tags="<pantry_reference>/txrevise/transcript_tags.txt",
    output:
        "<pantry_reference>/txrevise/txrevise_annotations.rds",
    params:
        script=script("txrevise/prepareAnnotations.R"),
    shell:
        """
        Rscript {params.script} \
            --gtf {input.gtf} \
            --tags {input.tags} \
            --out {output}
        """


rule txrevise_construct_events:
    """Construct events for annotations"""
    input:
        "<pantry_reference>/txrevise/txrevise_annotations.rds",
    output:
        "<pantry_reference>/txrevise/batch/txrevise.grp_1.upstream.{batch}_{n_batches}.gff3",
        "<pantry_reference>/txrevise/batch/txrevise.grp_2.upstream.{batch}_{n_batches}.gff3",
        "<pantry_reference>/txrevise/batch/txrevise.grp_1.contained.{batch}_{n_batches}.gff3",
        "<pantry_reference>/txrevise/batch/txrevise.grp_2.contained.{batch}_{n_batches}.gff3",
        "<pantry_reference>/txrevise/batch/txrevise.grp_1.downstream.{batch}_{n_batches}.gff3",
        "<pantry_reference>/txrevise/batch/txrevise.grp_2.downstream.{batch}_{n_batches}.gff3",
    resources:
        runtime="16h",
    params:
        script=script("txrevise/constructEvents.R"),
        batch_str="'{batch} {n_batches}'",
        outdir=parent_dir,
    shell:
        """
        Rscript {params.script} \
            --annot {input} \
            --batch {params.batch_str} \
            --out {params.outdir}
        """


rule txrevise_merge_gff_files:
    """Merge txrevise output files"""
    input:
        gff=expand(
            "<pantry_reference>/txrevise/batch/txrevise.{{group}}.{{position}}.{batch}_{n_batches}.gff3",
            batch=[i for i in range(1, TXREVISE_N_BATCHES + 1)],
            n_batches=TXREVISE_N_BATCHES,
        ),
    output:
        gff="<pantry_reference>/txrevise/txrevise.{group}.{position}.gff3",
    shell:
        'cat {input.gff} | grep -v "^#" > {output.gff}'


rule txrevise_make_one:
    """Iterate over groups and positions"""
    input:
        gff=expand(
            "<pantry_reference>/txrevise/txrevise.{group}.{position}.gff3",
            group=["grp_1", "grp_2"],
            position=["upstream", "contained", "downstream"],
        ),
        # metadata = 'processed/{annotation}_{kind}/txrevise_{kind}_phenotype_metadata.tsv.gz'
    output:
        "<pantry_reference>/txrevise/completed.txt",
    shell:
        "echo 'Done!' > {output}"


rule alt_TSS_polyA_transcript_seqs:
    """Generate transcript sequences from txrevise annotations for kallisto."""
    input:
        ref_genome=ref_genome,
        gff3="<pantry_reference>/txrevise/txrevise.{group}.{position}.gff3",
    output:
        "<pantry_reference>/txrevise/txrevise.{group}.{position}.fa.gz",
    shell:
        "gffread {input.gff3} -g {input.ref_genome} -w - | bgzip -c > {output}"


rule alt_TSS_polyA_kallisto_index:
    """Generate the index for kallisto using annotations from txrevise."""
    input:
        "<pantry_reference>/txrevise/txrevise.{group}.{position}.fa.gz",
    output:
        "<pantry_reference>/txrevise/txrevise.{group}.{position}.idx",
    resources:
        mem_mb=32000,
    shell:
        "kallisto index --index {output} {input}"


rule alt_TSS_polyA_kallisto:
    """Quantify expression using annotations from txrevise."""
    input:
        index="<pantry_reference>/txrevise/txrevise.{group}.{position}.idx",
        fastq=fastq_inputs,
    output:
        h5="<pantry_intermediate>/alt_TSS_polyA/{group}.{position}/{sample_id}/abundance.h5",
        tsv="<pantry_intermediate>/alt_TSS_polyA/{group}.{position}/{sample_id}/abundance.tsv",
        json="<pantry_intermediate>/alt_TSS_polyA/{group}.{position}/{sample_id}/run_info.json",
    threads: 16
    resources:
        runtime="16h",
    params:
        alt_group_pos_dir=lambda w, output: str(Path(output.h5).parent.parent),
        out_dir=parent_dir,
        single_end_flag=kallisto_single_end_flag,
        # TODO add strandedness parameter
    shell:
        """
        mkdir -p {params.alt_group_pos_dir}
        kallisto quant \
            --index {input.index} \
            --output-dir {params.out_dir} \
            {params.single_end_flag} \
            --threads {threads} \
            {input.fastq}
        """


rule assemble_alt_TSS_polyA_bed:
    """Assemble kallisto log2 and tpm outputs into alt TSS/polyA BED files"""
    input:
        kallisto=lambda w: expand(
            "<pantry_intermediate>/alt_TSS_polyA/grp_{grp}.{position}/{sample_id}/abundance.tsv",
            grp=[1, 2],
            position=dict(TSS="upstream", polyA="downstream")[w.type],
            sample_id=samples,
        ),
        samples=samples_file,
        ref_anno=ref_anno,
    output:
        "<pantry_output>/unnorm/alt_{type}.bed",
    resources:
        mem_mb=32000,
    params:
        script=script("assemble_bed.py"),
        unnorm_dir=parent_dir,
        # input.kallisto lists all grp_1 samples, then all grp_2 samples.
        alt_group1_pos_dir=lambda w, input: str(Path(input.kallisto[0]).parent.parent),
        alt_group2_pos_dir=lambda w, input: str(Path(input.kallisto[-1]).parent.parent),
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.script} alt-tss-polya \
            --group1-dir {params.alt_group1_pos_dir} \
            --group2-dir {params.alt_group2_pos_dir} \
            --samples {input.samples} \
            --ref-anno {input.ref_anno} \
            --output {output}
        """


rule normalize_alt_TSS_polyA:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed="<pantry_output>/unnorm/alt_{type}.bed",
        samples=samples_file,
    output:
        "<pantry_output>/alt_{type}.bed.gz",
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


rule pheno_groups_alt_TSS_polyA:
    """Group phenotypes by gene for tensorQTL"""
    input:
        "<pantry_output>/alt_{type}.bed.gz",
    output:
        "<pantry_output>/alt_{type}.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\\t" g }}' \
                >{output}
        """
