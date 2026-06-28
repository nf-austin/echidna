process RUN_ECHIDNA {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nf-austin/echidna:1.0.3'

    input:
    tuple val(sample_id), path(h5ad), path(wgs_csv)

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
    python3 ${moduleDir}/run_echidna.py \\
        --h5ad ${h5ad} \\
        --sample_id ${sample_id} \\
        --timepoint_label ${params.timepoint_label} \\
        --counts_layer ${params.counts_layer} \\
        --clusters ${params.clusters} \\
        --n_steps ${params.n_steps} \\
        --learning_rate ${params.learning_rate} \\
        --val_split ${params.val_split} \\
        --seed ${params.seed} \\
        --inverse_gamma ${params.inverse_gamma} \\
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
}
