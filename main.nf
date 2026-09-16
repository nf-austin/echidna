#!/usr/bin/env nextflow

include { DOWNLOAD_GENE_BED } from './modules/download_gene_bed/main'
include { SEG_TO_GENE_CN }    from './modules/seg_to_gene_cn/main'
include { RUN_ECHIDNA }       from './modules/run_echidna/main'
include { CALL_TUMOR_CELLS }  from './modules/call_tumor_cells/main'
include { CONCAT_H5ADS }      from './modules/concat_h5ads/main'

def helpMessage() {
    log.info """
    nf-austin/echidna -- Bayesian gene dosage inference from scRNA-seq (+ optional bulk WGS)

    Usage, samplesheet (recommended; this is what Seqera Platform launches with):
      nextflow run main.nf -profile docker --input samplesheet.csv --outdir results

      samplesheet.csv columns: sample, h5ad, seg_txt (seg_txt optional -- leave
      blank to run that sample in no-WGS mode)

    Usage, auto-discovery from sibling pipeline output directories:
      nextflow run main.nf -profile docker \\
          --scrna_dir /path/to/scrnaseq/results \\
          --wgs_dir   /path/to/wgs-cna/results

    Required (one of):
      --input         Samplesheet CSV.
      --scrna_dir     nf-austin/scrnaseq output dir ({dir}/{sample}/{sample}_annotated.h5ad).

    Common options:
      --wgs_dir       nf-austin/wgs-cna output dir ({dir}/{sample}.seg.txt). Optional.
      --sample_map    CSV (scrna_sample,wgs_sample) bridging mismatched naming.
      --gene_bed      Existing refGene BED. Prefer this on HPC -- compute nodes
                      often cannot reach UCSC.
      --genome        UCSC genome for the BED auto-download (default: ${params.genome}).
      --call_tumor_cells  Derive aneuploid/diploid calls and a CIN diversity index
                      (default: ${params.call_tumor_cells}).
      --outdir        Output directory (default: ${params.outdir}).

    On a cluster, add `-profile slurm,singularity --slurm_queue <partition>` and use
    absolute paths. See the README for the HPC notes.
    """.stripIndent()
}

/**
 * Resolve one samplesheet entry to a file.
 *
 * A relative entry is resolved against the samplesheet's OWN directory first,
 * which is what someone editing that sheet expects. Nextflow's default is the
 * launch directory, and on Seqera Platform the launch directory is the work
 * directory -- so a relative path there silently resolves somewhere unrelated.
 * Falls back to launch-dir resolution, and only then reports the entry missing.
 */
def resolveInput(path, sheet_dir, row_num, column) {
    // Absolute POSIX path, or a remote URI (s3://, gs://, az://): take as-is.
    if (path.startsWith('/') || path ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/) {
        return file(path, checkIfExists: true)
    }
    def beside_sheet = sheet_dir.resolve(path)
    if (beside_sheet.exists()) {
        return beside_sheet
    }
    def from_launch = file(path)
    if (from_launch.exists()) {
        return from_launch
    }
    error "Samplesheet row ${row_num}: '${column}' not found as '${beside_sheet}' (relative to the samplesheet) nor as '${from_launch}' (relative to the launch directory). Use an absolute path."
}

/**
 * Turn samplesheet rows into (sample_id, h5ad, seg_txt|null) tuples.
 *
 * Validated eagerly over the fully-read row list rather than inside a channel
 * closure: errors raised in a closure are lazy -- they never fire under
 * -preview, and in a real run they surface only once the channel is consumed.
 * A bad samplesheet must fail at launch, before Platform provisions compute.
 */
def buildSamples(rows, sheet_dir) {
    if (!rows) {
        error "Samplesheet is empty: ${params.input}"
    }
    def required = ['sample', 'h5ad']
    def missing = required.findAll { c -> !rows[0].containsKey(c) }
    if (missing) {
        error "Samplesheet is missing column(s): ${missing.join(', ')}. Found: ${rows[0].keySet().join(', ')}"
    }

    def seen = [] as Set
    return rows.withIndex().collect { row, idx ->
        def sample_id = row.sample?.trim()
        if (!sample_id) {
            error "Samplesheet row ${idx + 1} has an empty 'sample' value"
        }
        if (!seen.add(sample_id)) {
            error "Samplesheet has a duplicate sample id: '${sample_id}'. Sample ids become output paths and must be unique."
        }
        def h5ad = row.h5ad?.trim()
        if (!h5ad) {
            error "Samplesheet row ${idx + 1} ('${sample_id}') has an empty 'h5ad' value"
        }
        // seg_txt is optional: a blank cell, or no column at all, means no-WGS mode.
        def seg = row.seg_txt?.trim()
        tuple(
            sample_id,
            resolveInput(h5ad, sheet_dir, idx + 1, 'h5ad'),
            seg ? resolveInput(seg, sheet_dir, idx + 1, 'seg_txt') : null
        )
    }
}

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    // Validate eagerly, before any channel is built. Previously, giving neither
    // --scrna_dir nor --input reached channel.fromPath(null) and failed with an
    // unhelpful NullPointerException once the channel was consumed.
    if (params.scrna_dir && params.input) {
        error "Use either --scrna_dir (auto-discovery) or --input (samplesheet), not both."
    }
    if (!params.scrna_dir && !params.input) {
        error "No input given. Provide --input samplesheet.csv or --scrna_dir /path/to/scrnaseq/results. Run with --help for details."
    }
    if (params.wgs_dir && params.input) {
        error "--wgs_dir applies to --scrna_dir auto-discovery only. With --input, put the seg.txt path in the samplesheet's seg_txt column."
    }
    if (params.sample_map && !params.scrna_dir) {
        error "--sample_map applies to --scrna_dir auto-discovery only."
    }

    // ── Input discovery ───────────────────────────────────────────────────────
    // Mode 1: auto-discover from nf-austin/scrnaseq and nf-austin/wgs-cna output dirs
    //         --sample_map CSV bridges mismatched naming conventions
    // Mode 2: explicit samplesheet CSV (sample,h5ad,seg_txt) for multi-timepoint or custom inputs
    if (params.scrna_dir) {
        ch_h5ad = channel.fromPath("${params.scrna_dir}/*/*_annotated.h5ad", checkIfExists: true)
            | map { f -> tuple(f.parent.name, f) }

        if (params.wgs_dir) {
            ch_seg = channel.fromPath("${params.wgs_dir}/*.seg.txt", checkIfExists: true)
                | map { f -> tuple(f.name.replaceFirst(/\.seg\.txt$/, ''), f) }

            if (params.sample_map) {
                // Mapping CSV columns: scrna_sample, wgs_sample
                // Blank wgs_sample → run that scRNA sample in no-WGS mode
                // scRNA samples absent from map → run in no-WGS mode (not skipped)
                // WGS samples absent from map → ignored
                def map_file = file(params.sample_map, checkIfExists: true)
                def map_rows = map_file.splitCsv(header: true, strip: true)
                if (map_rows && !map_rows[0].containsKey('scrna_sample')) {
                    error "--sample_map needs a 'scrna_sample' column. Found: ${map_rows[0].keySet().join(', ')}"
                }
                ch_map = channel.fromList(
                    map_rows.collect { row -> tuple(row.scrna_sample, row.wgs_sample ?: null) }
                )

                // Left-join from h5ad: every scRNA sample is processed;
                // map provides an optional WGS sample name. Map-only rows (no h5ad) are dropped.
                ch_mapped = ch_h5ad.join(ch_map, remainder: true)
                    // [scrna_sample, h5ad, wgs_sample_or_null]
                    .filter { _id, h5ad, _wgs -> h5ad != null }

                ch_mapped.branch { row ->
                    has_wgs_name: row[2] != null
                    no_wgs_name:  true
                }.set { ch_map_branched }

                // Rekey by wgs_sample to look up seg files; fall back to no-WGS if not found
                ch_with_seg = ch_map_branched.has_wgs_name
                    .map    { scrna, h5ad, wgs -> tuple(wgs, scrna, h5ad) }
                    .join   (ch_seg, remainder: true)
                    .filter { row -> row[1] != null }
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
        // Explicit samplesheet — required for multi-timepoint (pre-concatenated
        // h5ads). Read and validated synchronously; see buildSamples above.
        def sheet = file(params.input, checkIfExists: true)
        def rows = sheet.splitCsv(header: true, strip: true)
        ch_input = channel.fromList(buildSamples(rows, sheet.parent))
    }

    log.info """
    P I P E L I N E   nf-austin/echidna
    ===================================
    input        : ${params.input ?: params.scrna_dir}
    wgs          : ${params.wgs_dir ?: (params.input ? 'from samplesheet seg_txt' : 'none')}
    gene bed     : ${params.gene_bed ?: "auto-download (${params.genome})"}
    tumor calls  : ${params.call_tumor_cells}
    outdir       : ${params.outdir}
    """.stripIndent()

    // ── Branch on WGS availability ────────────────────────────────────────────
    ch_input.branch { row ->
        with_wgs:    row[2] != null
        without_wgs: true
    }.set { ch_branched }

    // ── Gene BED — use provided file or auto-download from UCSC ──────────────
    if (params.gene_bed) {
        ch_gene_bed = channel.value(file(params.gene_bed, checkIfExists: true))
    } else {
        ch_gene_bed = DOWNLOAD_GENE_BED(channel.value(params.genome)).bed
    }

    SEG_TO_GENE_CN(
        ch_branched.with_wgs.map { id, _h5ad, seg -> tuple(id, file(seg)) },
        ch_gene_bed,
        channel.value(file("${projectDir}/modules/seg_to_gene_cn/seg_to_gene_cn.py"))
    )

    ch_with_w = ch_branched.with_wgs
        .map    { id, h5ad, _seg -> tuple(id, h5ad) }
        .join   (SEG_TO_GENE_CN.out.wgs_csv)
        .map    { id, h5ad, wcsv -> tuple(id, h5ad, wcsv, params.inverse_gamma) }

    ch_without_w = ch_branched.without_wgs
        .map { id, h5ad, _null -> tuple(id, h5ad, [], true) }

    RUN_ECHIDNA(
        ch_with_w.mix(ch_without_w),
        channel.value(file("${projectDir}/modules/run_echidna/run_echidna.py"))
    )

    // ── Optional: aneuploid/diploid calling + CIN diversity index ────────────
    if (params.call_tumor_cells) {
        CALL_TUMOR_CELLS(
            RUN_ECHIDNA.out.h5ad.join(RUN_ECHIDNA.out.cnv),
            channel.value(file("${projectDir}/modules/call_tumor_cells/call_tumor_cells.py"))
        )
        ch_final_h5ad = CALL_TUMOR_CELLS.out.h5ad
    } else {
        ch_final_h5ad = RUN_ECHIDNA.out.h5ad
    }

    // ── Concatenate every sample's h5ad into one combined AnnData for scanpy ──
    CONCAT_H5ADS(
        ch_final_h5ad.map { _id, h5ad -> h5ad }.collect(),
        channel.value(file("${projectDir}/modules/concat_h5ads/concat_h5ads.py"))
    )
}
