process download_references {
    /*
    Download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${sample} - ${params.max_references} reference(s)"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: 'mash-dist.txt'

    input:
    set val(sample), val(single_end), file(fq), file(sample_sketch) from DOWNLOAD_REFERENCES
    file(refseq_sketch) from REFSEQ_SKETCH

    output:
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("genbank/*.gbk") optional true into CALL_VARIANTS_AUTO
    file("mash-dist.txt")
    file "${task.process}/*" optional true

    when:
    REFSEQ_SKETCH_FOUND == true

    shell:
    no_cache = params.no_cache ? '-N' : ''
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    total = params.max_references
    template(task.ext.template)

    stub:
    """
    mkdir fastqs
    mkdir genbank
    mkdir ${task.process}
    touch fastqs/${sample}.fastq.gz
    touch genbank/*.gbk
    touch ${task.process}/*
    touch mash-dist.txt
    """
}