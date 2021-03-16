process qc_final_summary {
    /* Run FASTQC on the cleaned up FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    set val(sample), val(single_end), file(fq), file(genome_size) from QC_FINAL_SUMMARY

    output:
    file "quality-control/*"
    file "${task.process}/*" optional true

    shell:
    template(task.ext.template)

    stub:
    """
    mkdir quality-control
    mkdir ${task.process}
    touch quality-control/*
    touch ${task.process}/*
    """
}