nextflow.enable.dsl = 2

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
    tuple val(sample), val(single_end), file(fq), file(sample_sketch)
    file(refseq_sketch)

    output:
    tuple val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("genbank/*.gbk"), emit:CALL_VARIANTS_AUTO, optional: true
    file("mash-dist.txt")
    file "${task.process}/*" optional true

    when:
    REFSEQ_SKETCH_FOUND == true

    shell:
    no_cache = params.no_cache ? '-N' : ''
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    total = params.max_references

    """
#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
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
        check-staging.py --fq1 !{fq[0]} --extra !{sample_sketch} --is_single
    else
        check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{sample_sketch}
    fi
fi

# Get Mash distance
echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
mash dist -t !{sample_sketch} !{refseq_sketch} | grep -v "query" | sort -k 2,2 > distances.txt

# Pick genomes to download
printf "accession\tdistance\tlatest_accession\tupdated\n" > mash-dist.txt
select-references.py distances.txt !{total} !{tie_break} >> mash-dist.txt

# Pick only latest accessions
grep -v distance mash-dist.txt | cut -f3 > download-list.txt

# Download genomes
echo "# ncbi-genome-download Version" >> ${LOG_DIR}/!{task.process}.versions
ncbi-genome-download --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
ncbi-genome-download bacteria -l complete -o ./ -F genbank -p !{task.cpus} -A download-list.txt -r !{params.max_retry} !{no_cache} > ${LOG_DIR}/ncbi-genome-download.out 2> ${LOG_DIR}/ncbi-genome-download.err

# Move and uncompress genomes
mkdir genbank_temp
find refseq -name "*.gbff.gz" | xargs -I {} mv {} genbank_temp/
rename 's/(GC[AF]_\\d+).*/\$1/' genbank_temp/*
mkdir genbank
ls genbank_temp/ | xargs -I {} sh -c 'gzip -cd genbank_temp/{} > genbank/!{sample}-{}.gbk'
rm -rf genbank_temp

if [ "!{params.keep_all_files}" == "false" ]; then
    # Remove intermediate GenBank files
    rm -rf refseq/
fi

# pass the FASTQs along
mkdir -p fastqs
if [[ -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
    fi
else
    if [ "!{single_end}" == "false" ]; then
        # Paired-End Reads
        cp !{fq[0]} fastqs/!{sample}_R1.fastq.gz
        cp !{fq[1]} fastqs/!{sample}_R2.fastq.gz
    else
        # Single-End Reads
        cp  !{fq[0]} fastqs/!{sample}.fastq.gz
    fi
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi

    """

    stub:
    """
    mkdir fastqs
    mkdir genbank
    mkdir ${task.process}
    touch fastqs/${sample}.fastq.gz
    touch genbank/*.gbk
    touch ${task.process}/${sample}
    touch mash-dist.txt
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
        params.sample_sketch
        ])
    TEST_PARAMS_CH2 = Channel.of([
        params.refseq_sketch
        ])
    download_references(TEST_PARAMS_CH,TEST_PARAMS_CH2)
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