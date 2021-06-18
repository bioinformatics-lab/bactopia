nextflow.enable.dsl = 2

process MINMER_QUERY {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${sample} - ${dataset_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.txt"

    input:
    tuple val(sample), val(single_end), file(fq), file(sourmash)
    each file(dataset)

    output:
    file "*.txt"
    file "${task.process}/*" optional true

    when:
    MINMER_DATABASES.isEmpty() == false

    shell:
    dataset_name = dataset.getName()
    mash_w = params.screen_w ? "-w" : ""
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template "minmer_query.sh"

    stub:
    dataset_name = dataset.getName()
    """
    mkdir ${task.process}
    touch ${sample}.txt
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        file(params.fq),
        file(params.sourmash)
    ])
    TEST_PARAMS_CH2 = Channel.of(file(params.k21),file(params.k31),file(params.k51),file(params.refseqk21))
    minmer_query(TEST_PARAMS_CH,TEST_PARAMS_CH2.collect())
}
