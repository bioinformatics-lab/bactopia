process estimate_genome_size {
    /* Estimate the input genome size if not given. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra) from ESTIMATE_GENOME_SIZE

    output:
    file "${sample}-genome-size-error.txt" optional true
    file("${sample}-genome-size.txt") optional true
    set val(sample), val(sample_type), val(single_end), 
        file("fastqs/${sample}*.fastq.gz"), file(extra), file("${sample}-genome-size.txt") optional true into QC_READS, QC_ORIGINAL_SUMMARY
    file "${task.process}/*" optional true

    shell:
    genome_size = SPECIES_GENOME_SIZE
    template(task.ext.template)

    stub:
    """
    mkdir fastqs
    mkdir ${task.process}
    touch ${sample}-genome-size-error.txt
    touch ${sample}-genome-size.txt
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    """
}