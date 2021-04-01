nextflow.enable.dsl = 2

process gather_fastqs {
    /* Gather up input FASTQs for analysis. */
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "bactopia.versions"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    tag "${sample}"

    input:
    tuple val(sample), val(sample_type), val(single_end), file(r1: '*???-r1'), file(r2: '*???-r2'), file(extra)// from create_input_channel(run_type)

    output:
    file("*-error.txt") optional true
    tuple val(sample), val(final_sample_type), val(single_end),
        file("fastqs/${sample}*.fastq.gz"), file("extra/*.gz"), emit: FASTQ_PE_STATUS optional true
    file("${task.process}/*") optional true
    file("bactopia.versions") optional true
    file("multiple-read-sets-merged.txt") optional true

    shell:
    bactopia_version = VERSION
    nextflow_version = nextflow.version
    is_assembly = sample_type.startsWith('assembly') ? true : false
    is_compressed = false
    no_cache = params.no_cache ? '-N' : ''
    use_ena = params.use_ena
    if (task.attempt >= 4) {
        if (use_ena) {
            // Try SRA
            use_ena = false 
        } else {
            // Try ENA
            use_ena = true
        }
    }
    if (extra) {
        is_compressed = extra.getName().endsWith('gz') ? true : false
    }
    section = null
    if (sample_type == 'assembly_accession') {
        section = sample.startsWith('GCF') ? 'refseq' : 'genbank'
    }
    fcov = params.coverage.toInteger() == 0 ? 150 : Math.round(params.coverage.toInteger() * 1.5)
    final_sample_type = sample_type
    if (sample_type == 'hybrid-merge-pe') {
        final_sample_type = 'hybrid'
    } else if (sample_type == 'merge-pe') {
        final_sample_type = 'paired-end'
    } else if (sample_type == 'merge-se') {
        final_sample_type = 'single-end'
    }

    template "gather_fastqs.sh"    
    stub:
    """
    mkdir fastqs
    mkdir extra
    mkdir ${task.process}
    touch *-error.txt
    touch fastqs/${sample}.fastq.gz
    touch extra/*.gz
    touch ${task.process}/*
    touch bactopia.versions
    touch multiple-read-sets-merged.txt
    """
}

//###############
//Module testing 
//###############

workflow test{
    
    test_params_input = Channel.of([
        params.sample, 
        params.sample_type, 
        params.single_end,
        params.run_type,
        params.r1,
        params.r2             
        ])

    gather_fastqs(test_params_input)
}
workflow.onComplete {

    println """

    gather_fastqs Test Execution Summary
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