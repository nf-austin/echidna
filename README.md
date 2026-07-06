# nf-austin/echidna

A Nextflow DSL2 pipeline wrapping [Echidna](https://github.com/azizilab/echidna) — a Bayesian framework for integrative inference of copy number alterations (CNAs) and gene dosage effects from scRNA-seq and bulk WGS data.

## Pipeline steps

1. **SEG_TO_GENE_CN** (`seg_to_gene_cn`) — converts ichorCNA segment-level copy numbers (`*.seg.txt` from [wgs-cna](https://github.com/nf-austin/wgs-cna)) to gene-level W matrix via overlap-weighted averaging against a gene annotation BED. Skipped for samples without WGS.
2. **RUN_ECHIDNA** (`run_echidna`) — runs pre-processing, SVI training, CNV inference (HMM/GMM), and gene dosage effect scoring. Produces per-sample outputs.
3. **CALL_TUMOR_CELLS** (`call_tumor_cells`, optional, `--call_tumor_cells true`) — derives a per-clone aneuploid/diploid call and CIN diversity index from RUN_ECHIDNA's CNV states, and propagates both onto every cell. See [Tumor cell calling & chromosomal instability](#tumor-cell-calling--chromosomal-instability).
4. **CONCAT_H5ADS** (`concat_h5ads`) — merges every sample's h5ad (post-`CALL_TUMOR_CELLS` if enabled, otherwise straight from `RUN_ECHIDNA`) into a single `combined_annotated.h5ad` for easy loading in scanpy, matching [nf-austin/copykat](https://github.com/nf-austin/copykat)'s final concatenation step.

## Requirements

- Nextflow >= 24.04.0
- Docker, Singularity, or Conda

## Usage

### Automatic integration with nf-austin/scrnaseq and nf-austin/wgs-cna

Point `--scrna_dir` and `--wgs_dir` at the output directories of the companion pipelines — no samplesheet or custom scripts needed. Samples are matched by name automatically; any sample missing a `seg.txt` runs in [no-WGS mode](#no-wgs-mode).

```bash
# scrnaseq + wgs-cna (full integration — gene BED downloaded automatically)
nextflow run nf-austin/echidna \
    -profile docker \
    --scrna_dir scrna_results/qc \
    --wgs_dir   wgs_results/cna/ichorcna

# scrnaseq only (no WGS — neutral diploid W for all samples)
nextflow run nf-austin/echidna \
    -profile docker \
    --scrna_dir     scrna_results/qc \
    --inverse_gamma true
```

The pipeline expects the standard output layouts produced by each pipeline:

| Pipeline | Expected layout |
| --- | --- |
| nf-austin/scrnaseq | `{scrna_dir}/{sample_id}/{sample_id}_annotated.h5ad` |
| nf-austin/wgs-cna | `{wgs_dir}/{sample_id}.seg.txt` |

`obs["passing_qc"]` cells from scrnaseq are filtered automatically. Raw counts in `X` are handled by `pre_process` without any `--counts_layer` override.

### Mismatched sample names (`--sample_map`)

If the scrnaseq and wgs-cna runs used different naming conventions, supply a two-column CSV that maps between them:

```bash
nextflow run nf-austin/echidna \
    -profile docker \
    --scrna_dir  scrna_results/qc \
    --wgs_dir    wgs_results/cna/ichorcna \
    --sample_map sample_map.csv
```

`sample_map.csv` format (`scrna_sample,wgs_sample`):

```csv
scrna_sample,wgs_sample
tumor1_scRNA,PATIENT1-WGS
tumor2_scRNA,PATIENT2-WGS
tumor3_scRNA,
```

Behaviour:
- **Both columns set** — scRNA sample is run with the matched WGS `.seg.txt`; falls back to no-WGS mode if the named `.seg.txt` is absent on disk
- **Blank `wgs_sample`** — scRNA sample runs in no-WGS mode (no WGS available)
- **scRNA sample not in map** — run in no-WGS mode; h5ad is still processed
- **WGS sample not in map** — ignored

### Multi-timepoint run (samplesheet mode)

For longitudinal analyses where multiple scrnaseq runs from the same patient must be modelled jointly, pre-concatenate the per-timepoint h5ads and supply them via a samplesheet:

```bash
nextflow run nf-austin/echidna \
    -profile docker \
    --input samplesheet.csv \
    --timepoint_label timepoint
```

Samplesheet format (`sample,h5ad,seg_txt`; `seg_txt` is optional):

```csv
sample,h5ad,seg_txt
patient1,patient1_combined.h5ad,wgs_results/cna/ichorcna/patient1.seg.txt
patient2,patient2_combined.h5ad,
```

To produce the combined h5ad from scrnaseq output:

```python
import anndata as ad
import scanpy as sc

pre  = sc.read_h5ad("scrna_results/qc/patient1_pre/patient1_pre_annotated.h5ad")
post = sc.read_h5ad("scrna_results/qc/patient1_post/patient1_post_annotated.h5ad")

# Filter before concat — run_echidna.py handles this automatically for single h5ads,
# but must be done here to avoid barcode collisions across timepoints
pre  = pre[pre.obs["passing_qc"]].copy()
post = post[post.obs["passing_qc"]].copy()

pre.obs["timepoint"]  = "pre"
post.obs["timepoint"] = "post"

ad.concat([pre, post], index_unique="-").write_h5ad("patient1_combined.h5ad")
```

## Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--scrna_dir` | `null` | scrnaseq output `qc/` directory for auto-discovery |
| `--wgs_dir` | `null` | wgs-cna output `cna/ichorcna/` directory for auto-discovery |
| `--sample_map` | `null` | CSV (`scrna_sample,wgs_sample`) bridging mismatched naming conventions |
| `--input` | `null` | Samplesheet CSV (alternative to `--scrna_dir`; required for multi-timepoint) |
| `--outdir` | `results` | Output directory |
| `--genome` | `hg38` | UCSC genome name; used to auto-download refGene BED when `--gene_bed` is not set |
| `--gene_bed` | `null` | Path to existing gene annotation BED; overrides auto-download |
| `--num_genes` | `null` | Highly variable genes to retain (null = keep all) |
| `--n_comps` | `15` | PCA components |
| `--phenograph_k` | `60` | k for PhenoGraph clustering |
| `--n_neighbors` | `15` | Neighbours for UMAP |
| `--timepoint_label` | `timepoint` | `adata.obs` column for timepoint |
| `--counts_layer` | `counts` | `adata.layers` key for raw counts |
| `--clusters` | `pheno_louvain` | `adata.obs` column for cluster assignments |
| `--n_steps` | `10000` | Max SVI iterations |
| `--learning_rate` | `0.1` | Adam learning rate |
| `--val_split` | `0.1` | Fraction held out for validation |
| `--patience` | `null` | Early stopping patience (null = disabled) |
| `--seed` | `42` | Random seed |
| `--inverse_gamma` | `false` | Inverse-Gamma prior on eta variance; automatically `true` for any sample without WGS |
| `--n_hmm_components` | `5` | HMM states for CNV inference |
| `--n_gmm_components` | `5` | GMM components for neutral CNA estimation |
| `--gaussian_smoothing` | `true` | Gaussian smoothing before HMM |
| `--filter_quantile` | `0.7` | Gene-level variance filter quantile |
| `--smoother_sigma` | `6` | Gaussian kernel sigma |
| `--smoother_radius` | `8` | Gaussian kernel radius |
| `--neut_method` | `peak` | Neutral GMM component method (`peak` or `mode`) |
| `--call_tumor_cells` | `false` | Run `CALL_TUMOR_CELLS` — aneuploid/diploid calling + CIN diversity index (see below) |
| `--tumor_frac_altered_threshold` | `0.05` | Clone-level fraction of non-neutral genes above which a clone is called `aneuploid` |
| `--max_memory` | `128.GB` | Resource cap |
| `--max_cpus` | `32` | Resource cap |
| `--max_time` | `72.h` | Resource cap |

## No-WGS mode

When a sample has no matching `seg.txt` (either `--wgs_dir` is unset or no file matches the sample name), a neutral diploid W matrix (all genes = 2.0) is used. The WGS likelihood term then anchors the **cluster-proportion-weighted average** of gene dosage to ≈2.0 per gene, while individual cluster-level dosages are still inferred from scRNA-seq correlations. Clone reconstruction and relative CNA inference still work; what is lost is absolute copy number anchoring. Set `--inverse_gamma true` when running without WGS.

Check `ichorCNA_summary.tsv` (from wgs-cna) before running — samples with `qc_status = FAIL` (MAD > 0.30) have unreliable copy number calls and should be treated as no-WGS.

## Tumor cell calling & chromosomal instability

With `--call_tumor_cells true`, an additional `CALL_TUMOR_CELLS` step reads RUN_ECHIDNA's
`{sample}_echidna_cnv.csv` and, per clone, computes:

- **`frac_altered`** — the fraction of genes with a non-neutral (`amp`/`del`) HMM call
- **`cnv_diversity_index`** — Shannon entropy (log2) of the clone's neut/amp/del call distribution
  across genes; a simple proxy for chromosomal instability (CIN). An all-neutral clone has entropy 0
  (stable); a clone with a broad mix of amplifications and deletions has higher entropy (unstable).
  A clone that's uniformly amplified or deleted genome-wide also scores near 0 despite being
  maximally altered — `cnv_diversity_index` measures *heterogeneity* of CNA calls, not *burden*, so
  it complements `frac_altered` rather than replacing it.
- **`echidna_prediction`** — `aneuploid` if `frac_altered` exceeds `--tumor_frac_altered_threshold`
  (default `0.05`), else `diploid`

Both values are then mapped back onto every cell via the `adata.obs` cluster column actually used
for training (recorded in `adata.uns['echidna']['config']['clusters']`, which may differ from
`--clusters` if `run_echidna.py` fell back to an auto-detected column). This mirrors the
`copykat_prediction`/`cnv_diversity_index` columns from the sibling
[nf-austin/copykat](https://github.com/nf-austin/copykat) pipeline (CopyKAT-based tumor/normal
calling) so results from both pipelines use consistent naming if ever compared.

`--tumor_frac_altered_threshold` is a manually-tuned heuristic, not a model-based classification —
there's no universal cutoff, so inspect `{sample}_echidna_clone_calls.csv` per cohort before
trusting the call, and adjust the threshold to the background noise level implied by your
`--n_hmm_components`/`--n_gmm_components`/`--filter_quantile` settings and gene BED resolution.

## Output structure

```text
results/
├── {sample}/
│   ├── {sample}_echidna.h5ad             # updated AnnData with .uns['echidna'] model results, adata.obs['sample'] set
│   ├── {sample}_echidna_cnv.csv          # per-gene CNV states per clone
│   ├── {sample}_gmm_neutrals.csv         # neutral state statistics per clone
│   ├── {sample}_gene_dosage.pt           # GDX variance ratios [genes × timepoints × clones]
│   ├── {sample}_echidna_clone_calls.csv  # [--call_tumor_cells] per-clone frac_altered, cnv_diversity_index, echidna_prediction
│   ├── {sample}_echidna_cell_calls.csv   # [--call_tumor_cells] per-cell echidna_prediction, cnv_diversity_index
│   └── {sample}_echidna_annotated.h5ad   # [--call_tumor_cells] echidna.h5ad + the two obs columns above
└── combined_annotated.h5ad               # every sample's h5ad concatenated (barcode collisions resolved via '-<index>' suffix)
```
