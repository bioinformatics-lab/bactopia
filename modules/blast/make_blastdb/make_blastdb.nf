nextflow.enable.dsl = 2

process make_blastdb {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    tuple val(sample), val(single_end), file(fasta)

    output:
    file("blastdb/*")
    tuple val(sample), file("blastdb/*") optional true// emit: BLAST_GENES, BLAST_PRIMERS, BLAST_PROTEINS, optional:true
    file "${task.process}/*" optional true

    shell:
    template "make_blastdb.sh" 

    stub:
    """
    mkdir blastdb
    mkdir ${task.process}
    touch blastdb/*
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
        params.fasta          
        ])

    make_blastdb(TEST_PARAMS_CH)
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

