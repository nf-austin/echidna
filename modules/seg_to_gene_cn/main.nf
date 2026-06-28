process SEG_TO_GENE_CN {
    tag { sample_id }

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(seg_txt)
    path gene_bed

    output:
    tuple val(sample_id), path("${sample_id}_W.csv"), emit: wgs_csv

    script:
    """
    python3 ${moduleDir}/seg_to_gene_cn.py \\
        --seg_txt ${seg_txt} \\
        --gene_bed ${gene_bed} \\
        --sample_id ${sample_id} \\
        --out_csv ${sample_id}_W.csv
    """
}
