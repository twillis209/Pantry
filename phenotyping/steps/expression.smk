rule kallisto_index:
    """Generate the index for kallisto."""
    input:
        ref_cdna,
    output:
        "<pantry_reference>/kallisto.idx",
    resources:
        mem_mb=32000,
    params:
        ref_dir=parent_dir,
    shell:
        """
        mkdir -p {params.ref_dir}
        kallisto index --index {output} {input}
        """


rule kallisto:
    """Quantify expression from FASTQ files using kallisto."""
    input:
        index="<pantry_reference>/kallisto.idx",
        fastq=fastq_inputs,
    output:
        h5="<pantry_intermediate>/expression/{sample_id}/abundance.h5",
        tsv="<pantry_intermediate>/expression/{sample_id}/abundance.tsv",
        json="<pantry_intermediate>/expression/{sample_id}/run_info.json",
    threads: 16
    resources:
        runtime="16h",
    params:
        expr_dir=lambda w, output: str(Path(output.h5).parent.parent),
        out_dir=parent_dir,
        single_end_flag=kallisto_single_end_flag,
        # TODO add strandedness parameter
    shell:
        """
        mkdir -p {params.expr_dir}
        kallisto quant \
            --index {input.index} \
            --output-dir {params.out_dir} \
            {params.single_end_flag} \
            --threads {threads} \
            {input.fastq}
        """


rule assemble_expression_bed:
    """Assemble kallisto outputs into expression BED files

    Currently only tpm is used for downstream analysis.
    """
    input:
        kallisto=expand(
            "<pantry_intermediate>/expression/{sample_id}/abundance.tsv",
            sample_id=samples,
        ),
        samples=samples_file,
        ref_anno=ref_anno,
    output:
        iso="<pantry_output>/unnorm/isoforms.bed",
        gene="<pantry_output>/unnorm/expression.bed",
    resources:
        mem_mb=32000,
    params:
        script=script("assemble_bed.py"),
        unnorm_dir=parent_dir,
        expr_dir=lambda w, input: str(Path(input.kallisto[0]).parent.parent),
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.script} expression \
            --input-dir {params.expr_dir} \
            --samples {input.samples} \
            --ref-anno {input.ref_anno} \
            --output-isoforms {output.iso} \
            --output-expression {output.gene}
        """


rule normalize_expression:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed="<pantry_output>/unnorm/expression.bed",
        samples=samples_file,
    output:
        "<pantry_output>/expression.bed.gz",
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


rule normalize_isoforms:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed="<pantry_output>/unnorm/isoforms.bed",
        samples=samples_file,
    output:
        "<pantry_output>/isoforms.bed.gz",
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


rule pheno_groups_isoforms:
    """Group phenotypes by gene for tensorQTL"""
    input:
        "<pantry_output>/isoforms.bed.gz",
    output:
        "<pantry_output>/isoforms.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
