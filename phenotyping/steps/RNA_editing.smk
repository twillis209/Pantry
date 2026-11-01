rule query_editing_level:
    input:
        edit_sites=ref_file("edit_sites_bed"),
        ref_genome=ref_file("genome"),
        bam="<pantry_intermediate>/{dataset}/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
    output:
        levels="<pantry_intermediate>/{dataset}/RNA_editing/edit_levels/{sample_id}.rnaeditlevel.tsv.gz",
    resources:
        runtime="8h",
    params:
        script=script("RNA_editing/query_editing_level.py"),
        edit_levels_dir=parent_dir,
    shell:
        """
        mkdir -p {params.edit_levels_dir}
        python3 {params.script} \
            --edit_sites {input.edit_sites} \
            --ref {input.ref_genome} \
            --bam {input.bam} \
            --output {output.levels}
        """


rule shared_sample_sites:
    input:
        levels=per_sample(
            "<pantry_intermediate>/{dataset}/RNA_editing/edit_levels/{sample_id}.rnaeditlevel.tsv.gz"
        ),
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        matrix="<pantry_intermediate>/{dataset}/RNA_editing/edit_site_matrix.tsv",
    params:
        script=script("RNA_editing/shared_samples_sites_matrix.py"),
        edit_levels_dir=lambda w, input: str(Path(input.levels[0]).parent),
        min_cov=lambda w: rna_editing(w)["edit_sites_min_coverage"],
        min_sam=rna_editing_min_samples,
    shell:
        """
        python3 {params.script} \
            --path_to_edit_files {params.edit_levels_dir} \
            --samples_file {input.samples} \
            --output_file {output.matrix} \
            --min_coverage {params.min_cov} \
            --min_samples {params.min_sam}
        """


rule map_edit_sites_to_genes:
    input:
        edit_sites=ref_file("edit_sites_bed"),
        ref_anno=ref_file("annotation"),
        anno_ok="<pantry_reference>/{reference}/annotation.ok",
    output:
        tsv="<pantry_reference>/{reference}/edit_sites_to_genes.tsv",
    localrule: True
    params:
        script=script("RNA_editing/map_edit_sites_to_genes.py"),
    shell:
        """
        python3 {params.script} {input.edit_sites} {input.ref_anno} {output.tsv}
        """


rule prepare_RNA_editing_phenotypes:
    input:
        matrix="<pantry_intermediate>/{dataset}/RNA_editing/edit_site_matrix.tsv",
        edit_sites_to_genes=ref_path("edit_sites_to_genes.tsv"),
    output:
        matrix="<pantry_intermediate>/{dataset}/RNA_editing/phenotype_matrix.tsv",
        site_map="<pantry_output>/{dataset}/RNA_editing.site_to_phenotype.tsv",
    localrule: True
    params:
        script=script("RNA_editing/prepare_rna_editing_phenotypes.py"),
        cluster=lambda w: "--cluster" if rna_editing(w)["edit_sites_cluster"] else "",
        correlation_threshold=lambda w: rna_editing(w)[
            "edit_sites_correlation_threshold"
        ],
        outdir=lambda w, output: str(Path(output.site_map).parent),
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.script} \
            --edit-matrix {input.matrix} \
            --site-to-gene {input.edit_sites_to_genes} \
            --output-matrix {output.matrix} \
            --output-site-map {output.site_map} \
            --correlation-threshold {params.correlation_threshold} \
            {params.cluster}
        """


rule assemble_RNA_editing_bed:
    """Convert RNA editing matrix into BED file"""
    input:
        matrix="<pantry_intermediate>/{dataset}/RNA_editing/phenotype_matrix.tsv",
        ref_anno=ref_file("annotation"),
    output:
        bed="<pantry_output>/{dataset}/unnorm/RNA_editing.bed",
    resources:
        mem_mb=32000,
    params:
        script=script("assemble_bed.py"),
        unnorm_dir=parent_dir,
    shell:
        """
        mkdir -p {params.unnorm_dir}
        python3 {params.script} rna-editing \
            --input {input.matrix} \
            --ref-anno {input.ref_anno} \
            --output {output.bed}
        """


rule normalize_RNA_editing:
    """Quantile-normalize values for QTL mapping"""
    input:
        bed="<pantry_output>/{dataset}/unnorm/RNA_editing.bed",
        samples="<pantry_intermediate>/{dataset}/samples.txt",
    output:
        "<pantry_output>/{dataset}/RNA_editing.bed.gz",
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


rule pheno_groups_RNA_editing:
    """Group phenotypes by gene for tensorQTL"""
    input:
        "<pantry_output>/{dataset}/RNA_editing.bed.gz",
    output:
        "<pantry_output>/{dataset}/RNA_editing.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
                >{output}
        """
