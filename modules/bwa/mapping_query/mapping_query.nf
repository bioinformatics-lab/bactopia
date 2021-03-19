process mapping_query {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${sample}}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "mapping/*"

    input:
    set val(sample), val(single_end), file(fq) from MAPPING_QUERY
    file(query) from Channel.from(MAPPING_FASTAS).collect()

    output:
    file "mapping/*"
    file "${task.process}/*" optional true

    when:
    MAPPING_FASTAS.isEmpty() == false

    shell:
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    mkdir mapping
    touch ${task.process}/*
    touch mapping/*
    """
}