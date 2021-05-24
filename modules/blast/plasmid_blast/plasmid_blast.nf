nextflow.enable.dsl = 2

process plasmid_blast {
    /*
    BLAST a set of predicted genes against the PLSDB BLAST database.
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{json,json.gz}"

    input:
    tuple val(sample), file(genes)
    file(blastdb_files)

    output:
    file("${sample}-plsdb.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    PLASMID_BLASTDB.isEmpty() == false

    shell:
    gunzip_genes = genes.getName().replace('.gz', '')
    blastdb = blastdb_files[0].getBaseName()
    template "plasmid_blast.sh"

    stub:
    """
    mkdir ${task.process}
    touch ${task.process}/${sample}
    touch ${sample}-plsdb.json
    touch ${sample}-plsdb.json.gz
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.genes,          
        ])
    TEST_PARAMS_CH2 = Channel.of([
        params.blastdb_files
        ])

    plasmid_blast(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}
