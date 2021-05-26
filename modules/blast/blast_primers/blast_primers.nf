nextflow.enable.dsl = 2

process blast_primers {
    /*
    Query primer FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "primers/*.{json,json.gz}"

    input:
    tuple val(sample), file(blastdb)
    file(query)

    output:
    file("primers/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PRIMER_FASTAS.isEmpty() == false

    shell:
    template "blast_primers.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir primers
    touch ${task.process}/${sample}
    touch primers/${sample}.json
    touch primers/${sample}.json.gz
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.blastdb,          
        ])
    TEST_PARAMS_CH2 = Channel.of(
        params.query
        )

    blast_primers(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}
