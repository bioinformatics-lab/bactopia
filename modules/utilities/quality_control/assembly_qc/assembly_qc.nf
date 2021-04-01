process assembly_qc {
    /* Assess the quality of the assembly using QUAST and CheckM */
    tag "${sample} - ${method}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/assembly", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${method}/*"

    input:
    set val(sample), file(fasta), file(genome_size) from ASSEMBLY_QC
    each method from Channel.fromList(['checkm', 'quast'])

    output:
    file "${method}/*"
    file "${task.process}/*" optional true

    shell:
    //CheckM Related
    full_tree = params.full_tree ? '' : '--reduced_tree'
    checkm_ali = params.checkm_ali ? '--ali' : ''
    checkm_nt = params.checkm_nt ? '--nt' : ''
    force_domain = params.force_domain ? '--force_domain' : ''
    no_refinement = params.no_refinement ? '--no_refinement' : ''
    individual_markers = params.individual_markers ? '--individual_markers' : ''
    skip_adj_correction = params.skip_adj_correction ? '--skip_adj_correction' : ''
    skip_pseudogene_correction = params.skip_pseudogene_correction ? '--skip_pseudogene_correction' : ''
    ignore_thresholds = params.ignore_thresholds ? '--ignore_thresholds' : ''
    template(task.ext.template)

    stub:
    """
    mkdir ${method}
    mkdir ${task.process}
    touch ${method}/*
    touch ${task.process}/*
    """
}