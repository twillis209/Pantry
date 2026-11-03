def genome_size(ref_genome: Path) -> int:
    """Estimate genome size from a FASTA file."""
    assert ref_genome.suffix != ".gz", "Reference genome must be uncompressed"
    return ref_genome.stat().st_size


def star_index_path(w) -> str:
    """The STAR index for the {dataset}'s reference and read length."""
    d = ds(w)
    return f"<pantry_reference>/{d['reference']}/star_index_{d['read_length']}/SAindex"


rule star_index:
    """Generate the index for STAR."""
    input:
        ref_genome=ref_file("genome"),
        ref_anno=ref_file("annotation"),
        anno_ok="<pantry_reference>/{reference}/annotation.ok",
    output:
        # Among other generated files:
        "<pantry_reference>/{reference}/star_index_{read_length}/SAindex",
        "<pantry_reference>/{reference}/star_index_{read_length}/chrNameLength.txt",
    conda:
        "../envs/star.yaml"
    threads: 16
    resources:
        mem_mb=60000,
        runtime="4h",
    params:
        index_dir=parent_dir,
        overhang=lambda w: int(w.read_length) - 1,
        genomeSAindexNbases=lambda w, input: int(
            np.log2(genome_size(Path(input.ref_genome))) / 2 - 1
        ),
    shell:
        """
        mkdir -p {params.index_dir}
        STAR --runMode genomeGenerate \
            --outTmpDir {params.index_dir}/tmp \
            --outFileNamePrefix {params.index_dir}/STAR_ \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.ref_genome} \
            --sjdbGTFfile {input.ref_anno} \
            --sjdbOverhang {params.overhang} \
            --genomeSAindexNbases {params.genomeSAindexNbases} \
            --runThreadN {threads}
        """


def fastq_star_param(wildcards):
    """Get a string listing the fastq input files, with two lists if paired_end.

    This is the string supplied directly to the STAR command.
    """
    paths, is_paired = fastq(wildcards.dataset)[wildcards.sample_id]
    if is_paired:
        return " ".join([",".join(paths[0]), ",".join(paths[1])])
    else:
        return ",".join(paths)


rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        fastq=fastq_inputs,
        index=star_index_path,
    output:
        "<pantry_intermediate>/{dataset}/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
    conda:
        "../envs/star.yaml"
    threads: 16
    resources:
        mem_mb=60000,
        runtime="8h",
    params:
        fastq_list=fastq_star_param,
        index_dir=lambda w, input: str(Path(input.index).parent),
        star_out_dir=parent_dir,
        prefix=lambda w, output: output[0].removesuffix("Aligned.sortedByCoord.out.bam"),
    shell:
        """
        mkdir -p {params.star_out_dir}
        STAR --runMode alignReads \
            --genomeDir {params.index_dir} \
            --readFilesIn {params.fastq_list} \
            --readFilesCommand 'zcat < ' \
            --twopassMode Basic \
            --outSAMstrandField intronMotif \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {threads}
        """


rule shrink_bam:
    """Remove SEQ and QUAL fields from BAM file to reduce size."""
    input:
        "<pantry_intermediate>/{dataset}/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
    output:
        "<pantry_intermediate>/{dataset}/bam/{sample_id}.bam",
    conda:
        "../envs/samtools.yaml"
    params:
        bam_dir=parent_dir,
    shell:
        """
        mkdir -p {params.bam_dir}
        samtools view -h {input} \
            | awk -v OFS="\t" '{{if (substr($0, 1, 1) != "@") {{$10="*"; $11="*"}}; print ;}}' \
            | samtools view -h -b \
                >{output}
        """
