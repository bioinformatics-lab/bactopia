nextflow.enable.dsl = 2

process blast_proteins {
    /*
    Query protein FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "proteins/*.{json,json.gz}"

    input:
    tuple val(sample), file(blastdb)
    file(query)
    output:
    file("proteins/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PROTEIN_FASTAS.isEmpty() == false

    shell:
    '''
#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
OUTDIR=proteins
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
echo "# tblastn Version" >> ${LOG_DIR}/!{task.process}.versions
tblastn -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

echo "# Parallel Version" >> ${LOG_DIR}/!{task.process}.versions
parallel --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

for fasta in *.fasta; do
    type=`readlink -f ${fasta}`
    name="${fasta%.*}"
    mkdir -p ${OUTDIR} temp_json
    cat ${fasta} | sed -e 's/<[^>]*>//g' |
    parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
    tblastn -db !{sample} \
            -outfmt 15 \
            -evalue 0.0001 \
            -qcov_hsp_perc !{params.qcov_hsp_perc} \
            -query - \
            -out temp_json/${name}_{#}.json

    merge-blast-json.py temp_json > ${OUTDIR}/${name}.json
    rm -rf temp_json

    if [[ !{params.compress} == "true" ]]; then
        pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.json
    fi
done

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
    """
    mkdir ${task.process}
    mkdir proteins
    touch ${task.process}/*
    touch proteins/*.json
    touch proteins/*.json.gz
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.blastdb,          
        ])
    TEST_PARAMS_CH2 = Channel.of([
        params.query
        ])

    blast_proteins(TEST_PARAMS_CH,TEST_PARAMS_CH2)
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