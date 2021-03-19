process qc_reads {
    /* Cleanup the reads using Illumina-Cleanup */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*error.txt"

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra), file(genome_size) from QC_READS

    output:
    file "*-error.txt" optional true
    file "quality-control/*"
    set val(sample), val(single_end),
        file("quality-control/${sample}*.fastq.gz") optional true into COUNT_31MERS, ARIBA_ANALYSIS,
                                                                       MINMER_SKETCH, CALL_VARIANTS,
                                                                       MAPPING_QUERY
    set val(sample), val(sample_type), val(single_end),
        file("quality-control/${sample}*.fastq.gz"), file(extra),
        file(genome_size) optional true into ASSEMBLY

    set val(sample), val(single_end),
        file("quality-control/${sample}*.{fastq,error-fq}.gz"),
        file(genome_size) optional true into QC_FINAL_SUMMARY
    file "${task.process}/*" optional true

    shell:
    qc_ram = task.memory.toString().split(' ')[0]
    is_assembly = sample_type.startsWith('assembly') ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapters = params.adapters ? file(params.adapters) : 'adapters'
    phix = params.phix ? file(params.phix) : 'phix'
    template(task.ext.template)

    stub:
    """
    mkdir quality-control
    mkdir ${task.process}
    touch *-error.txt
    touch quality-control/${sample}.fastq.gz
    touch quality-control/${sample}.error-fq.gz
    touch ${task.process}/*
    """
}