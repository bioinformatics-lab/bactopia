process fastq_status {
    /* Determine if FASTQs are PE or SE, and if they meet minimum basepair/read counts. */
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra) from FASTQ_PE_STATUS

    output:
    file "*-error.txt" optional true
    set val(sample), val(sample_type), val(single_end), 
        file("fastqs/${sample}*.fastq.gz"), file(extra) optional true into ESTIMATE_GENOME_SIZE
    file "${task.process}/*" optional true

    shell:
    single_end = fq[1] == null ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    template(task.ext.template)

    stub:
    """
    mkdir ${task.process}
    mkdir fastqs
    touch *-error.txt
    touch fastqs/${sample}.fastq.gz
    touch ${task.process}/*
    """
}