process minmer_query {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${sample} - ${dataset_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.txt"

    input:
    set val(sample), val(single_end), file(fq), file(sourmash) from MINMER_QUERY
    each file(dataset) from MINMER_DATABASES

    output:
    file "*.txt"
    file "${task.process}/*" optional true

    when:
    MINMER_DATABASES.isEmpty() == false

    shell:
    dataset_name = dataset.getName()
    mash_w = params.screen_w ? "-w" : ""
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)

    stub:
    dataset_name = dataset.getName()
    """
    mkdir ${task.process}
    touch *.txt
    touch ${task.process}/*
    """
}