nextflow.enable.dsl = 2

process qc_original_summary {
    /* Run FASTQC on the input FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    tuple val(sample), val(sample_type), val(single_end), file(fq), file(extra), file(genome_size)

    output:
    file "quality-control/*"
    file "${task.process}/*" optional true

    shell:
    template "qc_original_summary.sh"

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
        params.sample_type, 
        params.single_end,
        params.fq,
        params.extra, 
        params.genome_size            
        ])

    qc_original_summary(TEST_PARAMS_CH)
}
workflow.onComplete {

    println """

    qc_original_summary Test Execution Summary
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