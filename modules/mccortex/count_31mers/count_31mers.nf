nextflow.enable.dsl = 2

process count_31mers {
    /* Count 31mers in the reads using McCortex */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/kmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.ctx"

    input:
    tuple val(sample), val(single_end), file(fq)
    output:
    file "${sample}.ctx"
    file "${task.process}/*" optional true

    shell:
    m = task.memory.toString().split(' ')[0].toInteger() * 1000 - 500
    template "count_31mers.sh"

    stub:
    """
    mkdir ${task.process}
    touch ${sample}.ctx
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing 
//###############

workflow test{
    
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        params.fq
        ])

    count_31mers(TEST_PARAMS_CH)
}
