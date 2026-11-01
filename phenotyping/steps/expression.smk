rule kallisto_index:
    """Generate the index for kallisto."""
    input:
        "<pantry_reference>/{reference}/cDNA.fa.gz",
    output:
        "<pantry_reference>/{reference}/kallisto.idx",
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
        index=ref_path("kallisto.idx"),
        fastq=fastq_inputs,
    output:
        h5="<pantry_intermediate>/{dataset}/expression/{sample_id}/abundance.h5",
        tsv="<pantry_intermediate>/{dataset}/expression/{sample_id}/abundance.tsv",
        json="<pantry_intermediate>/{dataset}/expression/{sample_id}/run_info.json",
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
        kallisto=per_sample(
            "<pantry_intermediate>/{dataset}/expression/{sample_id}/abundance.tsv"
        ),
        samples="<pantry_intermediate>/{dataset}/samples.txt",
        ref_anno=ref_file("annotation"),
    output:
        iso="<pantry_output>/{dataset}/unnorm/isoforms.bed",
        gene="<pantry_output>/{dataset}/unnorm/expression.bed",
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
        bed="<pantry_output>/{dataset}/unnorm/expression.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/expression.bed.gz",
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
        bed="<pantry_output>/{dataset}/unnorm/isoforms.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/isoforms.bed.gz",
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
        "<pantry_output>/{dataset}/isoforms.bed.gz",
    output:
        "<pantry_output>/{dataset}/isoforms.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
