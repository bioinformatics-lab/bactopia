process mapping_query {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${sample}}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "mapping/*"

    input:
    tuple val(sample), val(single_end), file(fq)
    file(query)// from Channel.from(MAPPING_FASTAS).collect()

    output:
    file "mapping/*"
    file "${task.process}/*" optional true

    when:
    MAPPING_FASTAS.isEmpty() == false

    shell:
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    mkdir mapping
    touch ${task.process}/*
    touch mapping/*
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
        ])
    TEST_PARAMS_CH2 = Channel.of([
        params.query
        ])
    annotate_genome(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}
workflow.onComplete {

    println """

    assemble_genome Test Execution Summary
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
