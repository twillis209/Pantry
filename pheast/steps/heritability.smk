rule calculate_heritability_chr_batch:
    """Get heritability for phenotypes in one batch using plink and gcta64.

    Runs one batch of phenotypes from one chromosome. Batches are
    assigned using the last digit of the TSS position (i.e. the `end` column)
    so that phenotypes for the same gene are in the same batch.
    """
    input:
        bed=pheno_bed,
        geno=geno_files,
        covar="<pantry_intermediate>/{dataset}/covar/{modality}.covar.plink.tsv",
    output:
        hsq="<pantry_intermediate>/{dataset}/heritability/{modality}.hsq.{chrom}.batch{batch}.tsv",
    conda:
        "../envs/gcta.yaml"
    resources:
        mem_mb=8000,
        runtime="24h",
    params:
        script=script("heritability.py"),
        geno_prefix=geno_input_prefix,
        grm_dir=lambda w, output: str(
            Path(output.hsq).parent / f"grm_{w.modality}_{w.chrom}"
        ),
        tmp_dir=lambda w, output: str(
            Path(output.hsq).parent / f"tmp_{w.modality}_{w.chrom}"
        ),
    shell:
        """
        mkdir -p {params.grm_dir}
        mkdir -p {params.tmp_dir}
        python3 {params.script} \
            --bed {input.bed} \
            --geno {params.geno_prefix} \
            --covar {input.covar} \
            --chrom {wildcards.chrom} \
            --batch {wildcards.batch} \
            --grm-dir {params.grm_dir} \
            --tmp-dir {params.tmp_dir} \
            --output {output.hsq}
        """


def heritability_batches(w) -> list:
    """Per-chromosome, per-batch heritability files of the {dataset}'s {modality}."""
    chroms, _ = phenotype_summary(w)
    return expand(
        "<pantry_intermediate>/{dataset}/heritability/{modality}.hsq.{chrom}.batch{batch}.tsv",
        dataset=w.dataset,
        modality=w.modality,
        chrom=chroms,
        batch=range(10),
    )


rule combine_heritability_chr:
    """Combine per-chromosome heritability stats into one file.

    A modality with no phenotypes has no batches: the header is written alone (R33).
    """
    input:
        heritability_batches,
    output:
        hsq="<pantry_output>/{dataset}/heritability/{modality}.hsq.tsv",
    conda:
        "../envs/gcta.yaml"
    params:
        hsq_dir=parent_dir,
    shell:
        """
        mkdir -p {params.hsq_dir}
        if [ -z "{input}" ]; then
            echo "NOTE: no phenotypes; writing a header-only heritability table"
            printf 'phenotype_id\thsq\tse\tpval\n' >{output.hsq}
            exit 0
        fi
        FILES=({input})
        head -n 1 ${{FILES[0]}} >{output.hsq}
        for FILE in {input}; do
            tail -n +2 $FILE >>{output.hsq}
        done
        """
