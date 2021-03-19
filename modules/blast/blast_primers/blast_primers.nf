process blast_primers {
    /*
    Query primer FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "primers/*.{json,json.gz}"

    input:
    set val(sample), file(blastdb) from BLAST_PRIMERS
    file(query) from Channel.from(BLAST_PRIMER_FASTAS).collect()

    output:
    file("primers/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PRIMER_FASTAS.isEmpty() == false

    shell:
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    mkdir primers
    touch ${task.process}/*
    touch primers/*.json
    touch primers/*.json.gz
    """
}