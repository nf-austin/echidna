process CALL_TUMOR_CELLS {
    tag { sample_id }
    publishDir { "${params.outdir}/${sample_id}" }, mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container params.echidna_container

    input:
    // Staged as a path input, not referenced via ${moduleDir}, which is not
    // bind-mounted into the container.
    tuple val(sample_id), path(h5ad), path(cnv_csv)
    path run_script

    output:
    tuple val(sample_id), path("${sample_id}_echidna_clone_calls.csv"), emit: clone_calls
    tuple val(sample_id), path("${sample_id}_echidna_cell_calls.csv"),  emit: cell_calls
    tuple val(sample_id), path("${sample_id}_echidna_annotated.h5ad"),  emit: h5ad

    script:
    """
    python3 ${run_script} \\
        --h5ad ${h5ad} \\
        --cnv_csv ${cnv_csv} \\
        --sample_id ${sample_id} \\
        --frac_altered_threshold ${params.tumor_frac_altered_threshold}
    """

    stub:
    """
    echo 'clone,frac_altered,call' > ${sample_id}_echidna_clone_calls.csv
    echo 'cell,clone,call'         > ${sample_id}_echidna_cell_calls.csv
    touch ${sample_id}_echidna_annotated.h5ad
    """
}
