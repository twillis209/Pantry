rule regtools_junctions:
    """Run regtools junctions for LeafCutter"""
    input:
        bam="<pantry_intermediate>/{dataset}/bam/{sample_id}.bam",
        bai="<pantry_intermediate>/{dataset}/bam/{sample_id}.bam.bai",
    output:
        "<pantry_intermediate>/{dataset}/splicing/{sample_id}.junc",
    conda:
        "../envs/regtools.yaml"
    params:
        splice_dir=parent_dir,
        min_anchor_len=8,
        min_intron_len=50,
        max_intron_len=500000,
        strandedness="XS",
    shell:
        """
        mkdir -p {params.splice_dir}
        regtools junctions extract \
            -a {params.min_anchor_len} \
            -m {params.min_intron_len} \
            -M {params.max_intron_len} \
            -s {params.strandedness} \
            -o {output} \
            {input.bam}
        """


rule cluster_junctions:
    """Cluster splice junctions with common boundaries"""
    input:
        per_sample("<pantry_intermediate>/{dataset}/splicing/{sample_id}.junc"),
    output:
        "<pantry_intermediate>/{dataset}/splicing/leafcutter_perind_numers.counts.gz",
    conda:
        "../envs/python.yaml"
    params:
        script=script("leafcutter_cluster_regtools_py3.py"),
        juncfile_list=lambda w, output: str(Path(output[0]).parent / "juncfiles.txt"),
        splice_dir=parent_dir,
        max_intron_len=100000,
        min_clust_reads=30,
        min_clust_ratio=0.001,
    shell:
        """
        printf '%s\\n' {input} >{params.juncfile_list}
        python3 {params.script} \
            --juncfiles {params.juncfile_list} \
            --rundir {params.splice_dir} \
            --maxintronlen {params.max_intron_len} \
            --minclureads {params.min_clust_reads} \
            --mincluratio {params.min_clust_ratio}
        """


rule assemble_splicing_bed:
    """Convert leafcutter output into splicing BED file"""
    input:
        counts="<pantry_intermediate>/{dataset}/splicing/leafcutter_perind_numers.counts.gz",
        ref_anno=ref_file("annotation"),
    output:
        bed="<pantry_output>/{dataset}/unnorm/splicing.bed",
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
        python3 {params.script} splicing \
            --input {input.counts} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """


rule normalize_splicing:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed="<pantry_output>/{dataset}/unnorm/splicing.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/splicing.bed.gz",
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


rule pheno_groups_splicing:
    """Group phenotypes by gene for tensorQTL"""
    input:
        "<pantry_output>/{dataset}/splicing.bed.gz",
    output:
        "<pantry_output>/{dataset}/splicing.phenotype_groups.txt",
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
