## Cross-modality assembly. Mirrors scripts/combine_modalities.sh per dataset.


rule combine_modalities:
    """Concatenate all listed modalities into one BED, with '{modality}:' prefixed phenotype names."""
    input:
        lambda w: [
            f"<pantry_output>/{w.dataset}/{m}.bed.gz"
            for m in ds(w)["cross_modality"]["modalities"]
        ],
    output:
        bed="<pantry_output>/{dataset}/cross_modality.bed.gz",
        tbi="<pantry_output>/{dataset}/cross_modality.bed.gz.tbi",
        tmp=temp("<pantry_intermediate>/{dataset}/cross_modality_tmp.bed"),
    localrule: True
    params:
        modalities=lambda w: ds(w)["cross_modality"]["modalities"],
        bed=lambda w, output: output.bed.removesuffix(".gz"),
        indir=lambda w, output: str(Path(output.bed).parent),
    shell:
        """
        # zcat | head would fail on SIGPIPE under pipefail (as in the script, off)
        set +o pipefail
        rm -f {output.tmp}
        first=$(echo {params.modalities} | cut -d' ' -f1)
        # First confirm that all bed files have the same header:
        for modality in {params.modalities}; do
            diff <(zcat {params.indir}/$first.bed.gz | head -n 1) \
                <(zcat {params.indir}/$modality.bed.gz | head -n 1)
        done
        for modality in {params.modalities}; do
            echo $modality
            # Prepend modality name to phenotype names:
            zcat {params.indir}/$modality.bed.gz | tail -n+2 \
                | awk -v OFS='\\t' -v m=$modality '{{$4 = m":"$4; print}}' \
                    >>{output.tmp}
        done
        zcat {params.indir}/$first.bed.gz | head -n1 >{params.bed}
        sort -k1,1 -k2,2n {output.tmp} >>{params.bed}
        bgzip -f {params.bed}
        tabix -p bed {output.bed}
        """


rule pheno_groups_cross_modality:
    """Group cross-modality phenotypes by gene for tensorQTL."""
    input:
        "<pantry_output>/{dataset}/cross_modality.bed.gz",
    output:
        "<pantry_output>/{dataset}/cross_modality.phenotype_groups.txt",
    localrule: True
    shell:
        """
        zcat <{input} \
            | tail -n +2 \
            | cut -f4 \
            | awk -F':' '{{ g=$2; sub(/__.*$/, "", g); print $0 "\t" g }}' \
                >{output}
        """
