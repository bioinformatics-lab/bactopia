process sequence_type {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${sample} - ${schema} - ${method}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/mlst/${schema}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${method}/*"

    input:
    set val(sample), val(single_end), file(fq), file(assembly) from SEQUENCE_TYPE
    each file(dataset) from MLST_DATABASES

    output:
    file "${method}/*"
    file "${task.process}/*" optional true

    when:
    MLST_DATABASES.isEmpty() == false

    shell:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '').split('-')[1]
    schema = dataset_tarball.split('-')[0]
    noclean = params.ariba_no_clean ? "--noclean" : ""
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    template(task.ext.template)

    stub:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = file(dataset).getName()
    schema = dataset_tarball.split('-')[0]
    """
    mkdir ${method}
    mkdir ${task.process}
    touch ${method}/*
    touch ${task.process}/*
    """
}