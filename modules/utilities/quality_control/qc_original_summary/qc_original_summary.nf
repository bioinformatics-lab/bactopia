process qc_original_summary {
    /* Run FASTQC on the input FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra), file(genome_size) from QC_ORIGINAL_SUMMARY

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