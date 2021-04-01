nextflow.enable.dsl = 2

process fastq_status {
    /* Determine if FASTQs are PE or SE, and if they meet minimum basepair/read counts. */
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    tuple val(sample), val(sample_type), val(single_end), file(fq), file(extra)
    output:
    file "*-error.txt" optional true
    tuple val(sample), val(sample_type), val(single_end), 
        file("fastqs/${sample}*.fastq.gz"), file(extra),emit: ESTIMATE_GENOME_SIZE optional true
    file "${task.process}/*" optional true

    shell:
    single_end = fq[1] == null ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    template "fastq_status.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir fastqs
    touch *-error.txt
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    """
}

//###############
//Module testing 
//###############

workflow test{
    
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.sample_type, 
        params.single_end,
        params.fq,
        params.extra             
        ])

    fastq_status(TEST_PARAMS_CH)
}
workflow.onComplete {

    println """

    fastq_status Test Execution Summary
    ---------------------------
    Command Line    : ${workflow.commandLine}
    Resumed         : ${workflow.resume}

    Completed At    : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Exit Code       : ${workflow.exitStatus}
    Error Report    : ${workflow.errorReport ?: '-'}
    """
}
workflow.onError {
    println "This test wasn't successful, Error Message: ${workflow.errorMessage}"
}