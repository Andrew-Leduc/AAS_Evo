#!/usr/bin/env python3
"""
Shard an ESM-MSR input CSV by structure and concatenate the per-structure score
CSVs back together. ESM-MSR's inference.py processes ONE structure per execution
("Found multiple in CSV" assertion), so a multi-protein input must be split into
one CSV per `code` (UniProt accession / structure), scored as a SLURM array, then
merged.

Usage:
    python3 esm_msr_shard.py split  <input.csv> <shard_dir>
    python3 esm_msr_shard.py concat <score_shard_dir> <combined_out.csv>
"""
import argparse
import glob
import os
import pandas as pd


def split(inp, outdir):
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(inp)
    codes = sorted(df["code"].astype(str).unique())
    with open(os.path.join(outdir, "list.txt"), "w") as f:
        for c in codes:
            df[df["code"].astype(str) == c].to_csv(
                os.path.join(outdir, f"{c}.csv"), index=False)
            f.write(c + "\n")
    print(f"{len(codes)} shards ({len(df):,} rows) -> {outdir}")


def concat(scoredir, out):
    files = sorted(glob.glob(os.path.join(scoredir, "*.csv")))
    if not files:
        raise SystemExit(f"no shard outputs in {scoredir}")
    df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    df.to_csv(out, index=False)
    print(f"{len(files)} shards -> {len(df):,} rows -> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split");  s.add_argument("input"); s.add_argument("outdir")
    c = sub.add_parser("concat"); c.add_argument("scoredir"); c.add_argument("out")
    a = ap.parse_args()
    (split if a.cmd == "split" else concat)(
        *( (a.input, a.outdir) if a.cmd == "split" else (a.scoredir, a.out) ))
