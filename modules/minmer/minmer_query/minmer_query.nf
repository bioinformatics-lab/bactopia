nextflow.enable.dsl = 2

process minmer_query {
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
    '''
#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}/!{dataset_name}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

# Print captured STDERR incase of exit
function print_stderr {
    cat .command.err 1>&2
    ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
}
trap print_stderr EXIT

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check-staging.py --fq1 !{fq[0]} --extra !{sourmash} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{sourmash}
    fi
fi

if [ "!{dataset_name}" == "refseq-k21-s1000.msh" ]; then
    echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
    mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > !{sample}-refseq-k21.txt
    gzip -cd !{fastq} | \
    mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
    sort -gr >> !{sample}-refseq-k21.txt 2> ${LOG_DIR}/mash.err
elif [ "!{dataset_name}" == "plsdb.msh" ]; then
    echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
    mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > !{sample}-plsdb-k21.txt
    gzip -cd !{fastq} | \
    mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
    sort -gr >> !{sample}-plsdb-k21.txt 2> ${LOG_DIR}/mash.err
elif [ "!{dataset_name}" == "genbank-k21.json.gz" ]; then
    echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
    sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k21.txt 2> ${LOG_DIR}/sourmash.err
elif [ "!{dataset_name}" == "genbank-k31.json.gz" ]; then
    echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
    sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k31.txt 2> ${LOG_DIR}/sourmash.err
else
    echo "# Sourmash Version" >> ${LOG_DIR}/!{task.process}.versions
    sourmash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    sourmash lca gather !{sourmash} !{dataset} > !{sample}-genbank-k51.txt 2> ${LOG_DIR}/sourmash.err
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
    '''

    stub:
    dataset_name = dataset.getName()
    """
    mkdir ${task.process}
    touch *.txt
    touch ${task.process}/*
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.single_end, 
        params.fq,   
        params.sourmash
        ])
    TEST_PARAMS_CH2 = Channel.of([
        params.dataset
        ])
    minmer_query(TEST_PARAMS_CH,TEST_PARAMS_CH2)
}
workflow.onComplete {

    println """

    estimate_genome_size Test Execution Summary
    ---------------------------
    Command Line    : ${workflow.commandLine}
    Resumed         : ${workflow.resume}

    Completed At    : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Exit Code       : ${workflow.exitStatus}
    Error Report    : ${workflow.errorReport ?: '-'}
    """
}
workflow.onError {
    println "This test wasn't successful, Error Message: ${workflow.errorMessage}"
}