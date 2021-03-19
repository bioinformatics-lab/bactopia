process plasmid_blast {
    /*
    BLAST a set of predicted genes against the PLSDB BLAST database.
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{json,json.gz}"

    input:
    set val(sample), file(genes) from PLASMID_BLAST
    file(blastdb_files) from Channel.from(PLASMID_BLASTDB).toList()

    output:
    file("${sample}-plsdb.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    PLASMID_BLASTDB.isEmpty() == false

    shell:
    gunzip_genes = genes.getName().replace('.gz', '')
    blastdb = blastdb_files[0].getBaseName()
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    touch ${task.process}/*
    touch ${sample}-plsdb.json
    touch ${sample}-plsdb.json.gz
    """
}