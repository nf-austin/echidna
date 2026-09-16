process SEG_TO_GENE_CN {
    tag { sample_id }

    conda "${moduleDir}/environment.yml"
    container params.echidna_container

    input:
    // The script is staged as a path input rather than referenced via
    // ${moduleDir}: moduleDir is not bind-mounted into the container, so a
    // moduleDir reference is a file-not-found under -profile docker/singularity.
    tuple val(sample_id), path(seg_txt)
    path gene_bed
    path run_script

    output:
    tuple val(sample_id), path("${sample_id}_W.csv"), emit: wgs_csv

    script:
    """
    python3 ${run_script} \\
        --seg_txt ${seg_txt} \\
        --gene_bed ${gene_bed} \\
        --sample_id ${sample_id} \\
        --out_csv ${sample_id}_W.csv
    """

    stub:
    """
    echo 'gene,cn' > ${sample_id}_W.csv
    """
}
