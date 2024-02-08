/*
* Define channels
*/
def raw_reads = Channel
                    .fromFilePairs("s3://vx-nf-bucket/nf_chipseq/data/*_{1,2}.fastq.gz")
                    .map {item -> [item[0], item[1]]}

/*
 * Quality control and trimming - trim_galore
 */
process trim_galore {
    container 'nfcore/chipseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir2}/qc" , mode: 'copy', pattern: "SRR*_1_val_1.fq.gz"
    publishDir "${params.outdir2}/qc" , mode: 'copy', pattern: "SRR*_2_val_2.fq.gz"


    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("SRR*_1_val_1.fq.gz"), path("SRR*_2_val_2.fq.gz")

    """
    trim_galore \
        -j 10 \
        -q 20 \
        --phred33 \
        --fastqc \
        --stringency 3 \
        --length 20 \
        -e 0.1 \
        --paired ${fastq[0]} ${fastq[1]} \
        --gzip \
        --fastqc
    """
}

/*
 * Read alignment - bowtie2
 */
process bowtie2 {
    container 'staphb/bowtie2:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir2}/mapping", mode: 'copy', pattern: "*.sam"

    input:
    tuple val(name), path("trimmed_1_val_1.fq.gz"), path("trimmed_2_val_2.fq.gz")
    path index

    output:
    tuple val(name), path("${name}.sam")
 
    """
    bowtie2 \
        -p 10 \
        -x ${index}/grch38_1kgmaj \
        -1 trimmed_1_val_1.fq.gz \
        -2 trimmed_2_val_2.fq.gz \
        -S ${name}.sam
    """
}

/*
 * Sort bam - samtools
 */
process sortbam {
    container 'nfcore/chipseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir2}/mapping/sorted_bam", mode: 'copy', pattern: "*_sorted.bam"

    input:
    tuple val(name), path("query.sam")

    output:
    tuple val(name), path("${name}_sorted.bam")

    """
    samtools view -bS query.sam | samtools sort -o ${name}_sorted.bam
    """
}

/*
 * Remove duplicates - picard
 */
process mark_duplicates {
    container 'nfcore/chipseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir2}/mapping/picard_markdup", mode: 'copy', pattern: "*_sorted_marked.bam"
    publishDir "${params.outdir2}/mapping/picard_markdup", mode: 'copy', pattern: "*_picard_markdup_metrics.txt"

    input:
    tuple val(name), path("query.bam")

    output:
    tuple val(name), path("${name}_sorted_marked.bam"), emit: bamfiles
    path("${name}_picard_markdup_metrics.txt"), emit: markdup_metrics


    """
    picard MarkDuplicates \
        I=query.bam \
        O=${name}_sorted_marked.bam \
        M=${name}_picard_markdup_metrics.txt
    """
}

/*
 * Peak calling - macs2
 */
process peak_calling {
    container 'nfcore/chipseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir2}/peak_calling", mode: 'copy'

    input:
    tuple val(name), path("query.bam")

    output:
    tuple val(name), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(name), path("*.xls")                   , emit: xls
    tuple val(name), path("*.r")                     , emit: model
    tuple val(name), path("*.gappedPeak"), optional:true, emit: gapped
    tuple val(name), path("*.bed")       , optional:true, emit: bed
    tuple val(name), path("*.bdg")       , optional:true, emit: bdg

    """
    macs2 callpeak \
        -f BAM \
        -t query.bam \
        -n ${name} \
        -g hs \
        --bdg \
        -q 0.01
    """
}

/*
 * Convert bam to bigwig - bamCoverage 
 */
process convertTobw {
    container 'nfcore/chipseq:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir2}/bigwig", mode: 'copy'

    input:
    tuple val(name), path("query.bam")

    output:
    tuple val(name), path("${name}.bw")

    """
    samtools index -b query.bam 
    bamCoverage -bs 10 -b query.bam -o ${name}.bw 
    """
}

workflow {
    trim_galore(raw_reads)
    bowtie2(trim_galore.out, params.bowtie2_index) | sortbam | mark_duplicates
    convertTobw(mark_duplicates.out.bamfiles)
    peak_calling(mark_duplicates.out.bamfiles)

}



