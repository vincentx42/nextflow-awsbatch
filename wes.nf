/*
* Define channels
*/
def raw_reads = Channel
                    .fromFilePairs("s3://vx-nf-bucket/nf_wes/data/*_{1,2}.fastq.gz")
                    .map {item -> [item[0], item[1]]}

def bwa_index_dir = Channel
                        .fromPath(params.bwa_idx)
                        .collect()

/*
 * Quality control and trimming - trim_galore
 */
process trim_galore {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/qc" , mode: 'copy', pattern: "SRR*_1_val_1.fq.gz"
    publishDir "${params.outdir3}/qc" , mode: 'copy', pattern: "SRR*_2_val_2.fq.gz"


    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("SRR*_1_val_1.fq.gz"), path("SRR*_2_val_2.fq.gz")

    """
    trim_galore \
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
 * Read alignment - bwa mem
 */
process bwa_mem {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/mapping", mode: 'copy', pattern: "*.sam"
    publishDir "${params.outdir3}/mapping/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(name), path("trimmed_1_val_1.fq.gz"), path("trimmed_2_val_2.fq.gz")
    path index
    path ref

    output:
    tuple val(name), path("${name}.sam"), emit: samfiles
    path "${name}.log", emit: bwa_mem_logs

    """
    bwa mem \
        -t 10 \
        -M \
        -R "@RG\\tID:${name}\\tSM:${name}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1" \
        ${ref} \
        trimmed_1_val_1.fq.gz \
        trimmed_2_val_2.fq.gz \
        >${name}.sam \
        2>${name}.log
    """
}

/*
 * Sort bam - samtools
 */
process sortbam {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/mapping/sorted_bam", mode: 'copy', pattern: "*_sorted.bam"

    input:
    tuple val(name), path("query.sam")

    output:
    tuple val(name), path("${name}_sorted.bam")

    """
    samtools view -bS query.sam | samtools sort -o ${name}_sorted.bam
    """
}

/*
 * Mark duplicates - gatk
 */
process mark_duplicates {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/mapping/gatk_markdup", mode: 'copy', pattern: "*_sorted_marked.bam"
    publishDir "${params.outdir3}/mapping/gatk_markdup", mode: 'copy', pattern: "*_markdup_metrics.txt"

    input:
    tuple val(name), path("query.bam")

    output:
    tuple val(name), path("${name}_sorted_marked.bam"), emit: bamfiles
    path "${name}_markdup_metrics.txt", emit: markdup_metrics


    """
    gatk MarkDuplicates \
        -I query.bam \
        -O ${name}_sorted_marked.bam \
        -M ${name}_markdup_metrics.txt
    """
}

/*
 * Variant calling - gatk
 */
process base_recalibration {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/gatk_variant_calling", mode: 'copy', pattern: "*.table"

    input:
    tuple val(name), path("query.bam")
    path ref_genome
    path ref_vcf

    output:
    path "${name}_recal_data.table"



    """
    samtools faidx ${ref_genome}

    gatk CreateSequenceDictionary \
	    -R ${ref_genome} \
	    -O grch38.dict
    
    gatk IndexFeatureFile \
        -I ${ref_vcf}
    
    gatk BaseRecalibrator \
        -I query.bam \
        -R ${ref_genome} \
        --known-sites ${ref_vcf} \
        -O ${name}_recal_data.table
    """
}

process apply_BQSR {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/gatk_variant_calling", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(name), path("query.bam")
    path ref_genome
    path recalibration_table

    output:
    tuple val(name), path("${name}_sorted_marked_BQSR.bam")   

    """
    samtools faidx ${ref_genome}

    gatk CreateSequenceDictionary \
	    -R ${ref_genome} \
	    -O grch38.dict
    
    gatk ApplyBQSR \
        -R ${ref_genome} \
        -I query.bam \
        --bqsr-recal-file ${recalibration_table} \
        -O ${name}_sorted_marked_BQSR.bam
    """
}

process variant_calling {
    container 'nfcore/sarek:latest'
    queue 'nextflow_excess'
    publishDir "${params.outdir3}/gatk_variant_calling", mode: 'copy', pattern: "*.vcf"

    input:
    tuple val(name), path("query.bam")
    path ref_genome

    output:
    tuple val(name), path("${name}.vcf")   

    """
    samtools faidx ${ref_genome}

    samtools index query.bam
    
    gatk CreateSequenceDictionary \
	    -R ${ref_genome} \
	    -O grch38.dict
    
    gatk HaplotypeCaller \
        -R ${ref_genome} \
        -I query.bam \
        -O ${name}.vcf
    """
}

workflow {           
    trim_galore(raw_reads)
    bwa_mem(trim_galore.out, bwa_index_dir, params.hg38_genome_unzip).samfiles | sortbam | mark_duplicates
    base_recalibration(mark_duplicates.out.bamfiles, params.hg38_genome_unzip, params.ref_variant)
    apply_BQSR(mark_duplicates.out.bamfiles, params.hg38_genome_unzip, base_recalibration.out)
    variant_calling(apply_BQSR.out, params.hg38_genome_unzip)
}



