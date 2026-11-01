rule exonic_intronic_from_GTF:
    """Extract exonic and intronic regions from GTF"""
    input:
        ref_anno=ref_anno,
    output:
        exon_gtf="<pantry_reference>/exonic.gtf",
        intron_gtf="<pantry_reference>/intronic.gtf",
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
        bam="<pantry_intermediate>/bam/{sample_id}.bam",
        gtf="<pantry_reference>/{feature_type}.gtf",
    output:
        "<pantry_intermediate>/stability/{sample_id}.{feature_type}.counts.txt",
    threads: 8
    resources:
        runtime="16h",
    params:
        stab_dir=parent_dir,
        paired_flag=lambda w: "-p" if fastq_map[w.sample_id][1] else "",
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
        exon=expand(
            "<pantry_intermediate>/stability/{sample_id}.exonic.counts.txt",
            sample_id=samples,
        ),
        intron=expand(
            "<pantry_intermediate>/stability/{sample_id}.intronic.counts.txt",
            sample_id=samples,
        ),
        samples=samples_file,
        ref_anno=ref_anno,
    output:
        bed="<pantry_output>/unnorm/stability.bed",
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
        bed="<pantry_output>/unnorm/stability.bed",
        samples=samples_file,
    output:
        "<pantry_output>/stability.bed.gz",
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
