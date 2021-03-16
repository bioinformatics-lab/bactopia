process minmer_sketch {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{msh,sig}"

    input:
    set val(sample), val(single_end), file(fq) from MINMER_SKETCH

    output:
    file("${sample}*.{msh,sig}")
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("${sample}.sig") into MINMER_QUERY
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("${sample}-k31.msh") into DOWNLOAD_REFERENCES
    file "${task.process}/*" optional true

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)

    stub:
    """
    mkdir fastqs
    mkdir ${task.process}
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    touch ${sample}.sig
    touch ${sample}-k31.msh

    """
}