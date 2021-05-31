nextflow.enable.dsl = 2

process estimate_genome_size {
    /* Estimate the input genome size if not given. */
    tag "${sample}"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    tuple val(sample), val(sample_type), val(single_end), file(fq), file(extra)

    output:
    file "${sample}-genome-size-error.txt" optional true
    file("${sample}-genome-size.txt") optional true
    tuple val(sample), val(sample_type), val(single_end), 
        file("fastqs/${sample}*.fastq.gz"), file(extra), file("${sample}-genome-size.txt"),emit: QUALITY_CONTROL, optional: true
    file "${task.process}/*" optional true

    shell:
    genome_size = SPECIES_GENOME_SIZE
    
    template "estimate_genome_size.sh"

    stub:
    """
    mkdir fastqs
    mkdir ${task.process}
    touch ${sample}-genome-size-error.txt
    touch ${sample}-genome-size.txt
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.sample_type, 
        params.single_end,
        file(params.fq),
        file(params.extra)             
        ])

    estimate_genome_size(TEST_PARAMS_CH)
}