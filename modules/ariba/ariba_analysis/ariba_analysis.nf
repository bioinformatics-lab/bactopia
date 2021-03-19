process ariba_analysis {
    /* Run reads against all available (if any) ARIBA datasets */
    tag "${sample} - ${dataset_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/ariba", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${dataset_name}/*"

    input:
    set val(sample), val(single_end), file(fq) from ARIBA_ANALYSIS
    each file(dataset) from ARIBA_DATABASES

    output:
    file "${dataset_name}/*"
    file "${task.process}/*" optional true

    when:
    single_end == false && ARIBA_DATABASES.isEmpty() == false

    shell:
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    noclean = params.ariba_no_clean ? "--noclean" : ""
    template(task.ext.template)

    stub:
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    """
    mkdir ${dataset_name}
    mkdir ${task.process}
    touch ${dataset_name}/*
    touch ${task.process}/*
    """
}