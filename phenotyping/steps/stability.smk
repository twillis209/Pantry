def feature_gtf(w) -> str:
    """The {feature_type} GTF of the {dataset}'s reference."""
    return ref_path(f"{w.feature_type}.gtf")(w)


rule exonic_intronic_from_GTF:
    """Extract exonic and intronic regions from GTF"""
    input:
        ref_anno=ref_file("annotation"),
        anno_ok="<pantry_reference>/{reference}/annotation.ok",
    output:
        exon_gtf="<pantry_reference>/{reference}/exonic.gtf",
        intron_gtf="<pantry_reference>/{reference}/intronic.gtf",
    conda:
        "../envs/python.yaml"
    params:
        script=script("exonic_intronic_from_gtf.py"),
    shell:
        """
        python3 {params.script} \
            {input.ref_anno} \
            {output.exon_gtf} \
            {output.intron_gtf}
        """


rule run_featureCounts:
    """Run featureCounts from subread to get exon or intron read counts

    Intron counts represent reads that could only come from pre-mRNA, so they
    need only overlap an intron by one base. Exon counts represent reads that
    could come from pre-mRNA or mature mRNA, so they must be fully within an
    exon.
    """
    input:
        bam="<pantry_intermediate>/{dataset}/bam/{sample_id}.bam",
        gtf=feature_gtf,
    output:
        "<pantry_intermediate>/{dataset}/stability/{sample_id}.{feature_type}.counts.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    resources:
        runtime="16h",
    params:
        stab_dir=parent_dir,
        paired_flag=lambda w: "-p" if is_paired(w) else "",
        feature_id=lambda w: {"exonic": "exon", "intronic": "intron"}[w.feature_type],
        frac_overlap=lambda w: {"exonic": 1, "intronic": 0}[w.feature_type],
        # TODO add strandedness parameter
    shell:
        """
        mkdir -p {params.stab_dir}
        featureCounts \
            {input.bam} \
            {params.paired_flag} \
            -a {input.gtf} \
            -t {params.feature_id} \
            --fracOverlap {params.frac_overlap} \
            -T {threads} \
            -o {output}
        """


rule assemble_stability_bed:
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    input:
        exon=per_sample(
            "<pantry_intermediate>/{dataset}/stability/{sample_id}.exonic.counts.txt"
        ),
        intron=per_sample(
            "<pantry_intermediate>/{dataset}/stability/{sample_id}.intronic.counts.txt"
        ),
        samples="<pantry_intermediate>/{dataset}/samples.txt",
        ref_anno=ref_file("annotation"),
    output:
        bed="<pantry_output>/{dataset}/unnorm/stability.bed",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    params:
        script=script("assemble_bed.py"),
        unnorm_dir=parent_dir,
        stab_dir=lambda w, input: str(Path(input.exon[0]).parent),
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.script} stability \
            --input-dir {params.stab_dir} \
            --samples {input.samples} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """


rule normalize_stability:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed="<pantry_output>/{dataset}/unnorm/stability.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/stability.bed.gz",
    conda:
        "../envs/python.yaml"
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
