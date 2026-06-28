#!/usr/bin/env nextflow

include { SEG_TO_GENE_CN } from './modules/seg_to_gene_cn/main'
include { RUN_ECHIDNA }    from './modules/run_echidna/main'

workflow {
    // ── Input discovery ───────────────────────────────────────────────────────
    // Mode 1: auto-discover from nf-austin/scrnaseq and nf-austin/wgs-cna output dirs
    // Mode 2: explicit samplesheet CSV (sample,h5ad,seg_txt) for multi-timepoint or custom inputs
    if (params.scrna_dir) {
        // scrnaseq layout: {scrna_dir}/{sample_id}/{sample_id}_annotated.h5ad
        ch_h5ad = Channel.fromPath("${params.scrna_dir}/*/*_annotated.h5ad")
            | map { f -> tuple(f.parent.name, f) }

        if (params.wgs_dir) {
            // wgs-cna layout: {wgs_dir}/{sample_id}.seg.txt
            ch_seg = Channel.fromPath("${params.wgs_dir}/*.seg.txt")
                | map { f -> tuple(f.name.replaceFirst(/\.seg\.txt$/, ''), f) }
            // Left-join: samples without a matching seg.txt get null (no-WGS mode)
            ch_input = ch_h5ad
                .join(ch_seg, remainder: true)
                .filter { _id, h5ad, _seg -> h5ad != null }
        } else {
            ch_input = ch_h5ad.map { id, h5ad -> tuple(id, h5ad, null) }
        }
    } else {
        // Explicit samplesheet — required for multi-timepoint (pre-concatenated h5ads)
        ch_input = Channel.fromPath(params.input)
            | splitCsv(header: true)
            | map { row -> tuple(row.sample, file(row.h5ad), row.seg_txt ?: null) }
    }

    // ── Branch on WGS availability ────────────────────────────────────────────
    ch_input.branch {
        with_wgs:    it[2] != null
        without_wgs: true
    }.set { ch_branched }

    ch_gene_bed = Channel.value(file(params.gene_bed ?: 'NO_GENE_BED'))
    SEG_TO_GENE_CN(
        ch_branched.with_wgs.map { id, _h5ad, seg -> tuple(id, file(seg)) },
        ch_gene_bed
    )

    ch_with_w = ch_branched.with_wgs
        .map    { id, h5ad, _seg -> tuple(id, h5ad) }
        .join   (SEG_TO_GENE_CN.out.wgs_csv)
        .map    { id, h5ad, wcsv -> tuple(id, h5ad, wcsv) }

    ch_without_w = ch_branched.without_wgs
        .map { id, h5ad, _null -> tuple(id, h5ad, []) }

    ch_with_w.mix(ch_without_w) | RUN_ECHIDNA
}
