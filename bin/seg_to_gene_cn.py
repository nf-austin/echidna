#!/usr/bin/env python3
"""
Convert ichorCNA segment-level copy numbers to gene-level W CSV for echidna.

Input:  ichorCNA {sample}.seg.txt  (tab-separated; columns include chr, start, end,
                                     Corrected_Copy_Number)
        Gene annotation BED file   (≥4 columns: chr, start, end, gene_name)
Output: gene-indexed CSV with columns: gene, counts
        where counts = overlap-length-weighted mean Corrected_Copy_Number across
        all segments overlapping that gene.
"""

import argparse
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--seg_txt",   required=True)
    p.add_argument("--gene_bed",  required=True)
    p.add_argument("--sample_id", required=True)
    p.add_argument("--out_csv",   required=True)
    return p.parse_args()


def load_segments(path: str) -> pd.DataFrame:
    seg = pd.read_csv(path, sep="\t")
    # ichorCNA seg.txt may have a leading sample-name column; detect by header
    if seg.columns[0].lower() not in ("chr", "chrom", "chromosome"):
        # First column is sample ID — drop it
        seg = seg.iloc[:, 1:]
    # Normalise column names
    col_map = {}
    for c in seg.columns:
        lc = c.lower()
        if lc in ("chr", "chrom", "chromosome"):
            col_map[c] = "chr"
        elif lc == "start":
            col_map[c] = "start"
        elif lc == "end":
            col_map[c] = "end"
        elif "corrected_copy_number" in lc or lc == "corrected.copy.number":
            col_map[c] = "cn"
    seg = seg.rename(columns=col_map)[["chr", "start", "end", "cn"]]
    seg["start"] = seg["start"].astype(int)
    seg["end"]   = seg["end"].astype(int)
    seg["cn"]    = seg["cn"].astype(float)
    # Normalise chromosome names: ensure "chr" prefix
    seg["chr"] = seg["chr"].astype(str).str.strip()
    seg["chr"] = seg["chr"].apply(lambda c: c if c.startswith("chr") else f"chr{c}")
    return seg


def load_genes(path: str) -> pd.DataFrame:
    genes = pd.read_csv(path, sep="\t", header=None,
                        usecols=[0, 1, 2, 3],
                        names=["chr", "start", "end", "gene"])
    genes["start"] = genes["start"].astype(int)
    genes["end"]   = genes["end"].astype(int)
    genes["chr"]   = genes["chr"].astype(str).str.strip()
    genes["chr"]   = genes["chr"].apply(lambda c: c if c.startswith("chr") else f"chr{c}")
    # Deduplicate genes by taking the union span per gene name per chrom
    genes = (genes.groupby(["chr", "gene"], sort=False)
             .agg(start=("start", "min"), end=("end", "max"))
             .reset_index())
    return genes


def weighted_cn(gene_chr: str, gene_start: int, gene_end: int,
                segs: pd.DataFrame, diploid: float = 2.0) -> float:
    """Overlap-weighted mean segment copy number for one gene."""
    chr_segs = segs[segs["chr"] == gene_chr]
    if chr_segs.empty:
        return diploid
    # Compute overlap lengths
    ov_start = np.maximum(chr_segs["start"].values, gene_start)
    ov_end   = np.minimum(chr_segs["end"].values, gene_end)
    lengths   = np.maximum(0, ov_end - ov_start)
    total = lengths.sum()
    if total == 0:
        return diploid
    return float(np.dot(lengths, chr_segs["cn"].values) / total)


def main():
    args = parse_args()

    segs  = load_segments(args.seg_txt)
    genes = load_genes(args.gene_bed)

    counts = [
        weighted_cn(row.chr, row.start, row.end, segs)
        for row in genes.itertuples()
    ]

    out = pd.DataFrame({"gene": genes["gene"].values, "counts": counts})
    out = out.drop_duplicates(subset="gene").set_index("gene")
    out.to_csv(args.out_csv)
    print(f"Wrote {len(out)} gene-level copy numbers to {args.out_csv}", flush=True)


if __name__ == "__main__":
    main()
