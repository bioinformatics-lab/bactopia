nextflow.enable.dsl = 2

process blast_proteins {
    /*
    Query protein FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "proteins/*.{json,json.gz}"

    input:
    tuple val(sample), file(blastdb)
    file(query)
    output:
    file("proteins/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PROTEIN_FASTAS.isEmpty() == false

    shell:
    
    template "blast_proteins.sh"

    stub:
    """
    mkdir ${task.process}
    mkdir proteins
    touch ${task.process}/${sample}
    touch proteins/${sample}.json
    touch proteins/${sample}.json.gz
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

    blast_proteins(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}