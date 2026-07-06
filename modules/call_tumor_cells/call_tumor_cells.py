#!/usr/bin/env python3
"""
Derive per-clone aneuploid/diploid calls and a CIN (Shannon diversity) index
from echidna's echi_cnv output, then propagate both onto every cell.

Input:  {sample}_echidna_cnv.csv  (from RUN_ECHIDNA / ec.tl.echi_cnv) — one row
                                   per gene, with 'states_echidna_clone_i' HMM
                                   calls (neut/amp/del) per clone i
        {sample}_echidna.h5ad     (from RUN_ECHIDNA) — adata.uns['echidna']
                                   ['config']['clusters'] names the obs column
                                   actually used for training (may differ from
                                   the --clusters CLI flag if run_echidna.py
                                   fell back to an auto-detected column)
Output: {sample}_echidna_clone_calls.csv  — per clone: frac_altered,
                                             cnv_diversity_index, echidna_prediction
        {sample}_echidna_cell_calls.csv   — per cell: echidna_prediction,
                                             cnv_diversity_index
        {sample}_echidna_annotated.h5ad   — input h5ad with both columns added
                                             to adata.obs
"""

import argparse

import anndata as ad
import numpy as np
import pandas as pd

STATE_PREFIX = "states_echidna_clone_"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cnv_csv", required=True)
    p.add_argument("--sample_id", required=True)
    p.add_argument("--frac_altered_threshold", type=float, required=True)
    return p.parse_args()


def shannon_diversity(states: pd.Series) -> float:
    """Shannon entropy (log2, matching nf-austin/copykat) of the neut/amp/del calls."""
    p = states.value_counts(normalize=True)
    return float(-(p * np.log2(p)).sum())


def main():
    args = parse_args()
    sid = args.sample_id

    cnv = pd.read_csv(args.cnv_csv)
    state_cols = [c for c in cnv.columns if c.startswith(STATE_PREFIX)]
    if not state_cols:
        raise ValueError(
            f"No '{STATE_PREFIX}*' columns found in {args.cnv_csv}. "
            f"Available columns: {list(cnv.columns)}. The installed sc-echidna "
            f"version's echi_cnv output format may not match what this script expects."
        )

    clone_calls = pd.DataFrame({
        "frac_altered":        {c: (cnv[c] != "neut").mean() for c in state_cols},
        "cnv_diversity_index": {c: shannon_diversity(cnv[c]) for c in state_cols},
    })
    clone_calls.index = (
        clone_calls.index.str.replace(STATE_PREFIX, "", regex=False).astype(int)
    )
    clone_calls.index.name = "clone"
    clone_calls["echidna_prediction"] = np.where(
        clone_calls["frac_altered"] > args.frac_altered_threshold, "aneuploid", "diploid"
    )
    clone_calls = clone_calls.sort_index()
    clone_calls.to_csv(f"{sid}_echidna_clone_calls.csv")

    adata = ad.read_h5ad(args.h5ad)

    try:
        clusters_col = adata.uns["echidna"]["config"]["clusters"]
    except (KeyError, TypeError) as e:
        raise ValueError(
            "Could not resolve adata.uns['echidna']['config']['clusters'] — this h5ad "
            "does not look like RUN_ECHIDNA output (echi_cnv/echidna_train must have run "
            "on it first)."
        ) from e

    # Clone index i in echi_cnv corresponds to the i-th element of
    # sorted(adata.obs[clusters_col].unique()) — echidna's create_z_pi() builds the
    # per-clone tensors via np.unique(), which sorts. Do not assume int(cluster_label) == i;
    # cluster labels may be strings, and lexicographic order != numeric order (e.g. "10" < "2").
    sorted_labels = sorted(adata.obs[clusters_col].unique())
    label_to_clone = {label: i for i, label in enumerate(sorted_labels)}

    cell_clone = adata.obs[clusters_col].map(label_to_clone)
    adata.obs["echidna_prediction"] = cell_clone.map(clone_calls["echidna_prediction"])
    adata.obs["cnv_diversity_index"] = cell_clone.map(clone_calls["cnv_diversity_index"])

    adata.obs[["echidna_prediction", "cnv_diversity_index"]].to_csv(
        f"{sid}_echidna_cell_calls.csv"
    )
    adata.write_h5ad(f"{sid}_echidna_annotated.h5ad")
    print(f"call_tumor_cells complete for {sid}", flush=True)


if __name__ == "__main__":
    main()
