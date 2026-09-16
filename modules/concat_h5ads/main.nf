process CONCAT_H5ADS {
    // Spans every sample, so no `tag`.
    publishDir "${params.outdir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container params.echidna_container

    input:
    path h5ads
    path run_script

    output:
    path "combined_annotated.h5ad", emit: combined_h5ad

    script:
    """
    python3 ${run_script} \\
        --inputs ${h5ads} \\
        --out_h5ad combined_annotated.h5ad
    """

    stub:
    """
    touch combined_annotated.h5ad
    """
}
