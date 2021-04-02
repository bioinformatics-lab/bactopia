nextflow.enable.dsl = 2

process minmer_sketch {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{msh,sig}"

    input:
    tuple val(sample), val(single_end), file(fq)

    output:
    file("${sample}*.{msh,sig}")
    tuple val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("${sample}.sig"),emit: MINMER_QUERY
    tuple val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("${sample}-k31.msh"),emit: DOWNLOAD_REFERENCES
    file "${task.process}/*" optional true

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template "minmer_sketch.sh"

    stub:
    """
    mkdir fastqs
    mkdir ${task.process}
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    touch ${sample}.sig
    touch ${sample}-k31.msh

    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.single_end, 
        params.fq         
        ])

    minmer_sketch(TEST_PARAMS_CH)
}
workflow.onComplete {

    println """

    estimate_genome_size Test Execution Summary
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