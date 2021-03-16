process annotate_genome {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "annotation/${sample}*"

    input:
    set val(sample), val(single_end), file(fq), file(fasta), file(total_contigs) from ANNOTATION
    file prokka_proteins from PROKKA_PROTEINS
    file prodigal_tf from PRODIGAL_TF

    output:
    file "annotation/${sample}*"
    set val(sample), file("annotation/${sample}.{ffn,ffn.gz}") optional true into PLASMID_BLAST
    set val(sample),
        file("annotation/${sample}.{ffn,ffn.gz}"),
        file("annotation/${sample}.{faa,faa.gz}") optional true into ANTIMICROBIAL_RESISTANCE
    file "${task.process}/*" optional true

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    contig_count = total_contigs.getName().replace('total_contigs_', '')
    genus = "Genus"
    species = "species"
    proteins = ""
    if (prokka_proteins.getName() != 'EMPTY_PROTEINS') {
        proteins = "--proteins ${prokka_proteins}"
        if (SPECIES.contains("-")) {
            genus = SPECIES.split('-')[0].capitalize()
            species = SPECIES.split('-')[1]
        } else {
            genus = SPECIES.capitalize()
            species = "spp."
        }
    }

    prodigal = ""
    if (prodigal_tf.getName() != 'EMPTY_TF' && !params.skip_prodigal_tf) {
        prodigal = "--prodigaltf ${prodigal_tf}"
    }

    compliant = params.compliant ? "--compliant" : ""
    locustag = "--locustag ${sample}"
    renamed = false
    // Contig ID must <= 37 characters
    if ("gnl|${params.centre}|${sample}_${contig_count}".length() > 37) {
        locustag = ""
        compliant = "--compliant"
        renamed = true
    }
    addgenes = params.nogenes ? "" : "--addgenes"
    addmrna = params.addmrna ? "--addmrna" : ""
    rawproduct = params.rawproduct ? "--rawproduct" : ""
    cdsrnaolap = params.cdsrnaolap ? "--cdsrnaolap" : ""
    norrna = params.norrna ? "--norrna" : ""
    notrna = params.notrna ? "--notrna" : ""
    rnammer = params.rnammer ? "--rnammer" : ""
    rfam = params.rnammer ? "--rfam" : ""
    template(task.ext.template)

    stub:
    """
    mkdir annotation
    mkdir ${task.process}
    touch annotation/${sample}*
    touch annotation/${sample}.ffn
    touch annotation/${sample}.ffn.gz
    touch annotation/${sample}.faa
    touch annotation/${sample}.faa.gz
    touch "${task.process}/*"
    """
}