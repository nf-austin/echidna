process RUN_ECHIDNA {
    tag { sample_id }
    publishDir { "${params.outdir}/${sample_id}" }, mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container params.echidna_container

    input:
    // The script is staged as a path input rather than referenced via
    // ${moduleDir}: moduleDir is not bind-mounted into the container, so a
    // moduleDir reference is a file-not-found under -profile docker/singularity.
    tuple val(sample_id), path(h5ad), path(wgs_csv), val(inverse_gamma)
    path run_script

    output:
    tuple val(sample_id), path("${sample_id}_echidna.h5ad"),     emit: h5ad
    tuple val(sample_id), path("${sample_id}_echidna_cnv.csv"),  emit: cnv
    tuple val(sample_id), path("${sample_id}_gmm_neutrals.csv"), emit: neutrals
    tuple val(sample_id), path("${sample_id}_gene_dosage.pt"),   emit: dosage

    script:
    def wgs_arg       = wgs_csv          ? "--wgs_csv ${wgs_csv}"                    : ""
    def patience_arg  = params.patience  != null ? "--patience ${params.patience}"   : ""
    def num_genes_arg = params.num_genes != null ? "--num_genes ${params.num_genes}" : ""
    """
    python3 ${run_script} \\
        --h5ad ${h5ad} \\
        --sample_id ${sample_id} \\
        --timepoint_label ${params.timepoint_label} \\
        --counts_layer ${params.counts_layer} \\
        --clusters ${params.clusters} \\
        --n_steps ${params.n_steps} \\
        --learning_rate ${params.learning_rate} \\
        --val_split ${params.val_split} \\
        --seed ${params.seed} \\
        --inverse_gamma ${inverse_gamma} \\
        --n_comps ${params.n_comps} \\
        --phenograph_k ${params.phenograph_k} \\
        --n_neighbors ${params.n_neighbors} \\
        --n_hmm_components ${params.n_hmm_components} \\
        --n_gmm_components ${params.n_gmm_components} \\
        --gaussian_smoothing ${params.gaussian_smoothing} \\
        --filter_quantile ${params.filter_quantile} \\
        --smoother_sigma ${params.smoother_sigma} \\
        --smoother_radius ${params.smoother_radius} \\
        --neut_method ${params.neut_method} \\
        --threads ${task.cpus} \\
        ${wgs_arg} \\
        ${patience_arg} \\
        ${num_genes_arg}
    """

    // Lets `nextflow run ... -stub-run` exercise channel wiring, the WGS/no-WGS
    // branch and the optional tumor-calling step without running inference.
    // Filenames must stay in sync with the output: block above.
    stub:
    """
    touch ${sample_id}_echidna.h5ad ${sample_id}_gene_dosage.pt
    echo 'gene,clone,cn'    > ${sample_id}_echidna_cnv.csv
    echo 'gene,neutral_val' > ${sample_id}_gmm_neutrals.csv
    """
}
