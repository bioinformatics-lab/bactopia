nextflow.enable.dsl = 2

process plasmid_blast {
    /*
    BLAST a set of predicted genes against the PLSDB BLAST database.
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{json,json.gz}"

    input:
    tuple val(sample), file(genes)
    file(blastdb_files)

    output:
    file("${sample}-plsdb.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    PLASMID_BLASTDB.isEmpty() == false

    shell:
    gunzip_genes = genes.getName().replace('.gz', '')
    blastdb = blastdb_files[0].getBaseName()
    '''
#!/bin/bash
set -e
set -u

LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}

echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

echo "# blastn Version" >> ${LOG_DIR}/!{task.process}.versions
blastn -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

echo "# Parallel Version" >> ${LOG_DIR}/!{task.process}.versions
parallel --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

if [[ !{params.compress} == "true" ]]; then
    gunzip -f !{genes}
fi

file_size=`cat !{gunzip_genes} | wc -c`
block_size=$(( file_size / !{task.cpus} / 2 ))
mkdir -p temp_json
cat !{gunzip_genes} | sed -e 's/<[^>]*>//g' | \
parallel --gnu --plain -j !{task.cpus} --block ${block_size} --recstart '>' --pipe \
blastn -db !{blastdb} \
       -outfmt 15 \
       -task blastn \
       -evalue 1 \
       -max_target_seqs !{params.max_target_seqs} \
       -perc_identity !{params.perc_identity} \
       -qcov_hsp_perc !{params.qcov_hsp_perc} \
       -query - \
       -out temp_json/!{sample}_{#}.json

merge-blast-json.py temp_json > !{sample}-plsdb.json
rm -rf temp_json


if [[ !{params.compress} == "true" ]]; then
    pigz --best -n -p !{task.cpus} !{sample}-plsdb.json
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
    """
    mkdir ${task.process}
    touch ${task.process}/*
    touch ${sample}-plsdb.json
    touch ${sample}-plsdb.json.gz
    """
}

//###############
//Module testing 
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample, 
        params.genes,          
        ])
    TEST_PARAMS_CH2 = Channel.of([
        params.blastdb_files
        ])

    plasmid_blast(TEST_PARAMS_CH,TEST_PARAMS_CH2)
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