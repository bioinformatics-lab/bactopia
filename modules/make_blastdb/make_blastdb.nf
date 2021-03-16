process make_blastdb {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    set val(sample), val(single_end), file(fasta) from MAKE_BLASTDB

    output:
    file("blastdb/*")
    set val(sample), file("blastdb/*") into BLAST_GENES, BLAST_PRIMERS, BLAST_PROTEINS
    file "${task.process}/*" optional true

    shell:
    template(task.ext.template)

    stub:
    """
    mkdir blastdb
    mkdir ${task.process}
    touch blastdb/*
    touch ${task.process}/*
    """
}