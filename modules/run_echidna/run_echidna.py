#!/usr/bin/env python3
"""
Nextflow wrapper for the echidna Bayesian CNA inference pipeline.

Runs (in sequence, sharing adata.uns state):
  1. ec.tl.pre_process  — normalization, PCA, clustering
  2. ec.tl.echidna_train — SVI training
  3. ec.tl.echi_cnv      — CNV inference (HMM/GMM)
  4. ec.tl.gene_dosage_effect — GDX scores

If --wgs_csv is omitted, a neutral diploid W matrix (all genes = 2.0) is
used. The WGS constraint then anchors the cluster-proportion-weighted average
of gene dosage to 2.0, allowing relative inter-cluster CNA differences to be
learned from scRNA-seq alone. In this mode --inverse_gamma true is recommended.
"""

import argparse
import shutil
from pathlib import Path

import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser()
    # I/O
    p.add_argument("--h5ad",       required=True)
    p.add_argument("--sample_id",  required=True)
    p.add_argument("--wgs_csv",    default=None)
    # Pre-processing
    p.add_argument("--num_genes",    type=int,   default=None)
    p.add_argument("--n_comps",      type=int,   default=15)
    p.add_argument("--phenograph_k", type=int,   default=60)
    p.add_argument("--n_neighbors",  type=int,   default=15)
    # Training
    p.add_argument("--timepoint_label", default="timepoint")
    p.add_argument("--counts_layer",    default="counts")
    p.add_argument("--clusters",        default="pheno_louvain")
    p.add_argument("--n_steps",     type=int,   default=10000)
    p.add_argument("--learning_rate", type=float, default=0.1)
    p.add_argument("--val_split",   type=float, default=0.1)
    p.add_argument("--patience",    type=int,   default=None)
    p.add_argument("--seed",        type=int,   default=42)
    p.add_argument("--inverse_gamma", type=lambda x: x.lower() == "true", default=False)
    p.add_argument("--threads",     type=int,   default=1)
    # CNV inference
    p.add_argument("--n_hmm_components",   type=int,   default=5)
    p.add_argument("--n_gmm_components",   type=int,   default=5)
    p.add_argument("--gaussian_smoothing", type=lambda x: x.lower() == "true", default=True)
    p.add_argument("--filter_quantile",    type=float, default=0.7)
    p.add_argument("--smoother_sigma",     type=float, default=6.0)
    p.add_argument("--smoother_radius",    type=float, default=8.0)
    p.add_argument("--neut_method",        default="peak")
    return p.parse_args()


def main():
    args = parse_args()
    sid  = args.sample_id

    import echidna as ec
    from echidna.tools import EchidnaConfig

    # ── Load scRNA-seq ────────────────────────────────────────────────────────
    adata = sc.read_h5ad(args.h5ad)

    # Filter QC-failing cells if the scrnaseq pipeline's passing_qc flag is present
    if "passing_qc" in adata.obs.columns:
        n_before = adata.n_obs
        adata = adata[adata.obs["passing_qc"]].copy()
        print(f"Filtered to {adata.n_obs}/{n_before} passing-QC cells", flush=True)

    # ── Load or synthesise W matrix ───────────────────────────────────────────
    if args.wgs_csv:
        wdf = pd.read_csv(args.wgs_csv, index_col=0)
        # Intersect and deduplicate genes
        shared = adata.var.index.intersection(wdf.index)
        wdf    = wdf.loc[shared]
        wdf    = wdf.loc[~wdf.index.duplicated(keep=False)]
        adata  = adata[:, wdf.index].copy()
    else:
        # Neutral diploid W: cluster-proportion-weighted average anchored to 2.0;
        # individual cluster dosages (eta) still inferred from scRNA-seq.
        wdf = pd.DataFrame({"counts": 2.0}, index=adata.var.index)

    # echidna's pre_process calls sc.pp.calculate_qc_metrics(layer="counts") internally;
    # scrnaseq stores raw counts in X rather than a named layer, so backfill if absent.
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Ensure a timepoint column exists (single-timepoint mode)
    if args.timepoint_label not in adata.obs.columns:
        adata.obs[args.timepoint_label] = "single_tp"

    # ── Pre-process ──────────────────────────────────────────────────────────
    adata = ec.tl.pre_process(
        adata,
        num_genes=args.num_genes,
        n_comps=args.n_comps,
        phenograph_k=args.phenograph_k,
        n_neighbors=args.n_neighbors,
    )

    # ── Train ────────────────────────────────────────────────────────────────
    config = EchidnaConfig(
        timepoint_label=args.timepoint_label,
        counts_layer=args.counts_layer,
        clusters=args.clusters,
        n_steps=args.n_steps,
        learning_rate=args.learning_rate,
        val_split=args.val_split,
        patience=args.patience,
        seed=args.seed,
        inverse_gamma=args.inverse_gamma,
    )
    ec.tl.echidna_train(adata, wdf, config=config)

    # ── CNV inference ─────────────────────────────────────────────────────────
    ec.tl.echi_cnv(
        adata,
        n_hmm_components=args.n_hmm_components,
        n_gmm_components=args.n_gmm_components,
        gaussian_smoothing=args.gaussian_smoothing,
        filter_genes=True,
        filter_quantile=args.filter_quantile,
        smoother_sigma=args.smoother_sigma,
        smoother_radius=args.smoother_radius,
        neut_method=args.neut_method,
    )

    # ── Gene dosage effect ───────────────────────────────────────────────────
    ec.tl.gene_dosage_effect(adata)

    # ── Rename outputs to sample-scoped names ────────────────────────────────
    save_data = adata.uns.get("echidna", {}).get("save_data", {})

    cnv_src = save_data.get("echi_cnv")
    if cnv_src and Path(cnv_src).exists():
        shutil.copy(cnv_src, f"{sid}_echidna_cnv.csv")

    neut_src = save_data.get("gmm_neutrals")
    if neut_src and Path(neut_src).exists():
        shutil.copy(neut_src, f"{sid}_gmm_neutrals.csv")

    dosage_src = save_data.get("gene_dosage")
    if dosage_src and Path(dosage_src).exists():
        shutil.copy(dosage_src, f"{sid}_gene_dosage.pt")

    # ── Save updated AnnData ──────────────────────────────────────────────────
    adata.write_h5ad(f"{sid}_echidna.h5ad")
    print(f"echidna complete for {sid}", flush=True)


if __name__ == "__main__":
    main()
