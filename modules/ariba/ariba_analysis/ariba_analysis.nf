nextflow.enable.dsl = 2

process ariba_analysis {
    /* Run reads against all available (if any) ARIBA datasets */
    tag "${sample} - ${dataset_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/ariba", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${dataset_name}/*"

    input:
    tuple val(sample), val(single_end), file(fq)
    each file(dataset)

    output:
    file "${dataset_name}/*"
    file "${task.process}/*" optional true

    when:
    single_end == false && ARIBA_DATABASES.isEmpty() == false

    shell:
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    noclean = params.ariba_no_clean ? "--noclean" : ""
    template "ariba_analysis.sh"

    stub:
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    """
    mkdir ${dataset_name}
    mkdir ${task.process}
    touch ${dataset_name}/*
    touch ${task.process}/*
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
    TEST_PARAMS_CH2 = Channel.of([
        params.dataset
        ])

    ariba_analysis(TEST_PARAMS_CH,TEST_PARAMS_CH2)
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