rule regtools_junctions:
    """Run regtools junctions for LeafCutter"""
    input:
        bam="<pantry_intermediate>/bam/{sample_id}.bam",
        bai="<pantry_intermediate>/bam/{sample_id}.bam.bai",
    output:
        "<pantry_intermediate>/splicing/{sample_id}.junc",
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
        expand("<pantry_intermediate>/splicing/{sample_id}.junc", sample_id=samples),
    output:
        "<pantry_intermediate>/splicing/leafcutter_perind_numers.counts.gz",
    params:
        script=script("leafcutter_cluster_regtools_py3.py"),
        juncfile_list=lambda w, output: Path(output[0]).parent / "juncfiles.txt",
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
        counts="<pantry_intermediate>/splicing/leafcutter_perind_numers.counts.gz",
        ref_anno=ref_anno,
    output:
        bed="<pantry_output>/unnorm/splicing.bed",
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
        bed="<pantry_output>/unnorm/splicing.bed",
        samples=samples_file,
    output:
        "<pantry_output>/splicing.bed.gz",
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
        "<pantry_output>/splicing.bed.gz",
    output:
        "<pantry_output>/splicing.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
