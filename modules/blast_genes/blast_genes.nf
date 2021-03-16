process blast_genes {
    /*
    Query gene FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "genes/*.{json,json.gz}"

    input:
    set val(sample), file(blastdb) from BLAST_GENES
    file(query) from Channel.from(BLAST_GENE_FASTAS).collect()

    output:
    file("genes/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_GENE_FASTAS.isEmpty() == false

    shell:
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    mkdir genes
    touch ${task.process}/*
    touch genes/*.json
    touch genes/*.json.gz
    """
}