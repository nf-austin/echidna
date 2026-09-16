# nf-austin/echidna

A Nextflow DSL2 pipeline wrapping [Echidna](https://github.com/azizilab/echidna) — a Bayesian framework for integrative inference of copy number alterations (CNAs) and gene dosage effects from scRNA-seq and bulk WGS data.

## Pipeline steps

1. **SEG_TO_GENE_CN** (`seg_to_gene_cn`) — converts ichorCNA segment-level copy numbers (`*.seg.txt` from [wgs-cna](https://github.com/nf-austin/wgs-cna)) to gene-level W matrix via overlap-weighted averaging against a gene annotation BED. Skipped for samples without WGS.
2. **RUN_ECHIDNA** (`run_echidna`) — runs pre-processing, SVI training, CNV inference (HMM/GMM), and gene dosage effect scoring. Produces per-sample outputs.
3. **CALL_TUMOR_CELLS** (`call_tumor_cells`, optional, `--call_tumor_cells true`) — derives a per-clone aneuploid/diploid call and CIN diversity index from RUN_ECHIDNA's CNV states, and propagates both onto every cell. See [Tumor cell calling & chromosomal instability](#tumor-cell-calling--chromosomal-instability).
4. **CONCAT_H5ADS** (`concat_h5ads`) — merges every sample's h5ad (post-`CALL_TUMOR_CELLS` if enabled, otherwise straight from `RUN_ECHIDNA`) into a single `combined_annotated.h5ad` for easy loading in scanpy, matching [nf-austin/copykat](https://github.com/nf-austin/copykat)'s final concatenation step.

## Requirements

- Nextflow >= 24.04.0
- Docker (local) or Singularity/Apptainer (HPC); Conda works as a fallback
- Optional: a SLURM cluster — see [HPC / SLURM](#hpc--slurm)

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

### Samplesheet mode (recommended on Seqera Platform; required for multi-timepoint)

A samplesheet is the explicit alternative to directory discovery. It is validated at launch, and it
is what Seqera Platform launches with — see [Seqera Platform](#seqera-platform-nextflow-tower).
It is also required for longitudinal analyses, where multiple scrnaseq runs from the same patient
must be modelled jointly: pre-concatenate the per-timepoint h5ads and supply them via the sheet.

```bash
nextflow run nf-austin/echidna \
    -profile docker \
    --input samplesheet.csv \
    --timepoint_label timepoint
```

Samplesheet format (`sample,h5ad,seg_txt`; `seg_txt` is optional -- leave it blank to run that
sample in [no-WGS mode](#no-wgs-mode)). `sample` must be unique. Relative paths are resolved
against the samplesheet's own directory first, then the launch directory; prefer absolute paths,
especially on Seqera Platform, where the launch directory is the work directory. A copy-pasteable
example lives in `assets/samplesheet_example.csv`.

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
| `--clusters` | `pheno_leiden` | `adata.obs` column for cluster assignments |
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
| `--echidna_container` | `ghcr.io/nf-austin/echidna:0.1.0` | Image used by every process |
| `--max_memory` | `128.GB` | Resource cap |
| `--max_cpus` | `32` | Resource cap |
| `--max_time` | `72.h` | Resource cap |
| `--max_forks_echidna` | `1` | Concurrent `RUN_ECHIDNA` jobs; raise if the GPU has headroom or you are on CPU |
| `--slurm_queue` | *(cluster default)* | SLURM partition (`sbatch --partition`). Used by `-profile slurm` |
| `--slurm_account` | *(none)* | SLURM account to charge (`sbatch --account`) |
| `--cluster_options` | *(none)* | Raw sbatch options added to every job, e.g. `--qos=long` or `--gres=gpu:1` |
| `--singularity_cache_dir` | `$NXF_SINGULARITY_CACHEDIR` | Shared directory for pulled images |
| `--conda_cache_dir` | `$NXF_CONDA_CACHEDIR` | Shared directory for conda environments |
| `--singularity_bind` | *(none)* | Extra bind mounts, comma-separated, e.g. `/mnt/gpfs,/scratch` |

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
├── combined_annotated.h5ad               # every sample's h5ad concatenated (barcode collisions resolved via '-<index>' suffix)
├── reference/                            # auto-downloaded refGene BED (skipped when --gene_bed is given)
└── pipeline_info/                        # Nextflow execution report, timeline, trace and DAG
```

## Seqera Platform (Nextflow Tower)

The repo ships everything Platform needs:

- **`nextflow_schema.json`** — renders the launch form. `--input` appears as a file picker wired to
  Data Explorer, options are grouped by pipeline stage, and tuning knobs are marked hidden so the
  default form stays short.
- **`assets/schema_input.json`** — the samplesheet contract (`sample`, `h5ad`, optional `seg_txt`),
  so a malformed sheet is caught before compute is provisioned.
- **`tower.yml`** — puts the per-sample CNV tables, clone/cell calls and the Nextflow execution
  report in the run's **Reports** tab.

To add it: **Pipelines → Add pipeline**, point at this repository, and pick a compute environment.

**Prefer `--input` over `--scrna_dir` on Platform.** A samplesheet gets a file browser and is
validated at launch; directory-glob discovery is convenient on the command line but gives Platform
nothing to check. Use **absolute paths** for `--input`, the files it references, and `--outdir`.

Samplesheet problems — a missing column, a duplicate sample id, an unreadable `h5ad` — and
mutually exclusive flag combinations are raised at launch, before Platform provisions any compute.

## HPC / SLURM

The `slurm` profile sets only the executor and queue, so it composes with an engine profile in
either order:

```bash
nextflow run nf-austin/echidna \
    -profile slurm,singularity \
    --slurm_queue normal \
    --input /mnt/gpfs/project/sheet.csv \
    --gene_bed /mnt/gpfs/refs/hg38_refGene.bed \
    --outdir /mnt/gpfs/project/results \
    --singularity_cache_dir /mnt/gpfs/shared/singularity
```

Points that matter on a cluster:

- **Pass `--gene_bed`.** `DOWNLOAD_GENE_BED` fetches the refGene table from UCSC, and compute nodes
  on most clusters have no outbound network. Download it once on a login node and point at it; the
  process fails with that instruction if the download cannot complete.
- **Use absolute paths** for `--input`, the files it lists, and `--outdir`. The data is expected to
  live on the shared filesystem; nothing here assumes object storage.
- **Put `--singularity_cache_dir` on shared storage.** `$HOME` is usually quota-limited and is not
  always mounted on compute nodes. `NXF_SINGULARITY_CACHEDIR` is honored if you would rather set it
  site-wide.
- **`--singularity_bind` is the escape hatch for symlinked filesystems.** `autoMounts` binds only
  the paths Nextflow resolved itself; if `/data` is a symlink to `/mnt/gpfs/...`, the container sees
  a dangling link and reports a missing file even though the host path is fine. Bind the real
  parent: `--singularity_bind /mnt/gpfs`.
- **Requesting a GPU for `RUN_ECHIDNA`.** The SLURM executor ignores Nextflow's `accelerator`
  directive, so ask for the GPU through sbatch options instead:
  `--cluster_options '--gres=gpu:1'`. `--max_forks_echidna` stays at 1 by default so concurrent
  jobs do not contend for one device; raise it if the GPU has headroom or you are running on CPU.
- **Seqera Platform already sets the executor** when you launch against a SLURM compute
  environment, so `-profile slurm` is mainly for launching by hand from a login node.

## Container image

Built from `modules/run_echidna/Dockerfile` and published to GHCR by `.github/workflows/docker.yml`
as `ghcr.io/nf-austin/echidna:<ver>` (amd64 + arm64). One image serves every process: the echidna
stack is a superset of what the other steps need, and it also carries the `wget` that
`DOWNLOAD_GENE_BED` uses.

The GHCR package must be **public** for `nextflow run` to pull it without credentials.

Every module also ships an `environment.yml`, so `-profile conda` remains a working fallback.

To test a change to the image before it is published:

```bash
docker build -t echidna-nf:test modules/run_echidna
nextflow run . -profile docker --input samplesheet.csv --echidna_container echidna-nf:test
```

## Notes

- `nextflow run . -stub-run --input samplesheet.csv` exercises the real channel wiring, the
  WGS/no-WGS branch and publishing with no containers and no data — useful on a laptop.
- `nextflow lint main.nf nextflow.config modules/*/main.nf` catches config errors that `-preview`
  accepts.
