process blast_proteins {
    /*
    Query protein FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "proteins/*.{json,json.gz}"

    input:
    set val(sample), file(blastdb) from BLAST_PROTEINS
    file(query) from Channel.from(BLAST_PROTEIN_FASTAS).collect()

    output:
    file("proteins/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PROTEIN_FASTAS.isEmpty() == false

    shell:
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    mkdir proteins
    touch ${task.process}/*
    touch proteins/*.json
    touch proteins/*.json.gz
    """
}
