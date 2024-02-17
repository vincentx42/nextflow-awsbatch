/*
* Define channels
*/
def raw_reads = Channel
                    .fromFilePairs("s3://vx-nf-bucket/nf_rna/data/*_{1,2}.fastq.gz")
                    .map {item -> [item[0], item[1]]}

/*
 * Quality control and trimming - trim_galore
 */
process trim_galore {
    container 'nfcore/rnaseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir}/qc" , mode: 'copy', pattern: "SRR*_1_val_1.fq.gz"
    publishDir "${params.outdir}/qc" , mode: 'copy', pattern: "SRR*_2_val_2.fq.gz"


    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("SRR*_1_val_1.fq.gz"), path("SRR*_2_val_2.fq.gz")

    """
    trim_galore \
        -j 8 \
        -q 20 \
        --phred33 \
        --fastqc \
        --stringency 3 \
        --length 20 \
        -e 0.1 \
        --paired ${fastq[0]} ${fastq[1]} \
        --gzip
    """
}

/*
 * Read alignment - hisat2
 */
process hisat2 {
    container 'nfcore/rnaseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir}/mapping", mode: 'copy', pattern: "*.sam"
    publishDir "${params.outdir}/mapping/hisat2_logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(name), path("trimmed_1_val_1.fq.gz"), path("trimmed_2_val_2.fq.gz")
    path index

    output:
    tuple val(name), path("${name}.sam"), emit: samfiles
    path "${name}.log", emit: hisat2logs

    """
    hisat2 \
        -p 8 \
        -x ${index}/genome_tran \
        -1 trimmed_1_val_1.fq.gz \
        -2 trimmed_2_val_2.fq.gz \
        -S ${name}.sam \
        --summary-file ${name}.log
    """
}

/*
 * Convert to .bam files and sort - samtools
 */
process samtools {
    container 'nfcore/rnaseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir}/mapping", mode: 'copy', pattern: "*_sorted.bam"

    input:
    tuple val(name), path("query.sam")

    output:
    tuple val(name), path("${name}_sorted.bam")

    """
    samtools view -bS query.sam | samtools sort -o ${name}_sorted.bam
    """
}

/*
 * Read summarization - featureCounts
 */
process feature_counts {
    container 'nfcore/rnaseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir}/expr", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(name), path("query.bam")
    path annotation

    output:
    path "${name}_counts.txt"


    """
    featureCounts \
            -a ${annotation} \
            -o ${name}_counts.txt \
            query.bam
    """
}

workflow {
    trim_galore(raw_reads)
    hisat2(trim_galore.out, params.hisat2_index)
    samtools(hisat2.out.samfiles)
    feature_counts(samtools.out, params.annotation)
}



