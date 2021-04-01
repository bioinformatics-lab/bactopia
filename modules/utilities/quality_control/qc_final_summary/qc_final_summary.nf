nextflow.enable.dsl = 2

process qc_final_summary {
    /* Run FASTQC on the cleaned up FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    tuple val(sample), val(single_end), file(fq), file(genome_size)

    output:
    file "quality-control/*"
    file "${task.process}/*" optional true

    shell:
    template "qc_final_summary.sh"

    stub:
    """
    mkdir quality-control
    mkdir ${task.process}
    touch quality-control/*
    touch ${task.process}/*
    """
}

//###############
//Module testing 
//###############

workflow test{
    
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.single_end,
        params.fq,
        params.genome_size            
        ])

    qc_final_summary(TEST_PARAMS_CH)
}
workflow.onComplete {

    println """

    qc_reads Test Execution Summary
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