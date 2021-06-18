nextflow.enable.dsl = 2

process SEQUENCE_TYPE {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${sample} - ${schema} - ${method}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/mlst/${schema}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${method}/*"

    input:
    tuple val(sample), val(single_end), file(fq), file(assembly)
    each file(dataset)

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

    template "sequence_type.sh"

    stub:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = file(dataset).getName()
    schema = dataset_tarball.split('-')[0]
    """
    mkdir ${method}
    mkdir ${task.process}
    touch ${method}/${sample}
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing
//###############

workflow test{

    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        file(params.fq),
        file(params.assembly)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        file(params.dataset_blast)
        file(params.dataset_ariba))

    sequence_type(TEST_PARAMS_CH,TEST_PARAMS_CH2.collect())
}
