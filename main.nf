#!/usr/bin/env nextflow

include { SEG_TO_GENE_CN } from './modules/seg_to_gene_cn/main'
include { RUN_ECHIDNA }    from './modules/run_echidna/main'

process DOWNLOAD_GENE_BED {
    storeDir "${params.outdir}/reference"
    container 'ubuntu:22.04'

    input:
    val genome

    output:
    path "${genome}_refGene.bed", emit: bed

    script:
    """
    wget -qO- "https://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/refGene.txt.gz" \\
        | gunzip -c \\
        | awk 'BEGIN{OFS="\\t"} {print \$3, \$5, \$6, \$13}' \\
        | sort -k1,1 -k2,2n \\
        > ${genome}_refGene.bed
    """
}

workflow {
    // ── Input discovery ───────────────────────────────────────────────────────
    // Mode 1: auto-discover from nf-austin/scrnaseq and nf-austin/wgs-cna output dirs
    //         --sample_map CSV bridges mismatched naming conventions
    // Mode 2: explicit samplesheet CSV (sample,h5ad,seg_txt) for multi-timepoint or custom inputs
    if (params.scrna_dir) {
        ch_h5ad = Channel.fromPath("${params.scrna_dir}/*/*_annotated.h5ad")
            | map { f -> tuple(f.parent.name, f) }

        if (params.wgs_dir) {
            ch_seg = Channel.fromPath("${params.wgs_dir}/*.seg.txt")
                | map { f -> tuple(f.name.replaceFirst(/\.seg\.txt$/, ''), f) }

            if (params.sample_map) {
                // Mapping CSV columns: scrna_sample, wgs_sample
                // Blank wgs_sample → run that scRNA sample in no-WGS mode
                // scRNA samples absent from map → run in no-WGS mode (not skipped)
                // WGS samples absent from map → ignored
                ch_map = Channel.fromPath(params.sample_map)
                    | splitCsv(header: true)
                    | map { row -> tuple(row.scrna_sample, row.wgs_sample ?: null) }

                // Left-join from h5ad: every scRNA sample is processed;
                // map provides an optional WGS sample name. Map-only rows (no h5ad) are dropped.
                ch_mapped = ch_h5ad.join(ch_map, remainder: true)
                    // [scrna_sample, h5ad, wgs_sample_or_null]
                    .filter { _id, h5ad, _wgs -> h5ad != null }

                ch_mapped.branch {
                    has_wgs_name: it[2] != null
                    no_wgs_name:  true
                }.set { ch_map_branched }

                // Rekey by wgs_sample to look up seg files; fall back to no-WGS if not found
                ch_with_seg = ch_map_branched.has_wgs_name
                    .map    { scrna, h5ad, wgs -> tuple(wgs, scrna, h5ad) }
                    .join   (ch_seg, remainder: true)
                    .filter { _wgs, scrna, _h5ad, _seg -> scrna != null }
                    .map    { _wgs, scrna, h5ad, seg -> tuple(scrna, h5ad, seg) }
                    // seg is null when mapped wgs_sample has no matching .seg.txt → no-WGS

                ch_no_seg = ch_map_branched.no_wgs_name
                    .map { scrna, h5ad, _wgs -> tuple(scrna, h5ad, null) }

                ch_input = ch_with_seg.mix(ch_no_seg)
            } else {
                // No mapping — match scRNA and WGS samples by name
                ch_input = ch_h5ad
                    .join(ch_seg, remainder: true)
                    .filter { _id, h5ad, _seg -> h5ad != null }
            }
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

    // ── Gene BED — use provided file or auto-download from UCSC ──────────────
    if (params.gene_bed) {
        ch_gene_bed = Channel.value(file(params.gene_bed))
    } else {
        ch_gene_bed = DOWNLOAD_GENE_BED(Channel.value(params.genome)).bed
    }

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
