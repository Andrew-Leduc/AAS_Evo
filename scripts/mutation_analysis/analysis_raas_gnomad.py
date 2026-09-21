#!/usr/bin/env python3
"""
#6 (part 2)  Average precursor RAAS by gnomAD allele frequency, with a
"not in gnomAD" bin.

Reads the VEP-annotated candidate-SNV VCF (from build_candidate_snvs.py + VEP),
pulls gnomADe_AF per SNV, takes the MAX AF over all SNVs that produce each swap
(acc|pos|wt|alt encoded in the VCF ID). A swap with no gnomAD-observed SNV -> the
"not in gnomAD" bin. Joins per-swap mean RAAS and compares.

Run VEP on the candidate VCF first (in the venv shell, on a compute node):

  DATA=/scratch/leduc.an/AAS_Evo
  apptainer exec \
    -B $DATA/SEQ_FILES/vep_cache:/cache \
    -B $DATA/SEQ_FILES:/ref \
    -B $DATA/ANALYSIS/contact_saap:/io \
    -B $DATA/SEQ_FILES:/am \
    /scratch/leduc.an/tools/vep/ensembl-vep.sif \
    vep -i /io/candidate_snvs.vcf -o /io/candidate_snvs.vep.vcf.gz \
        --vcf --compress_output bgzip --cache --dir_cache /cache \
        --assembly GRCh38 --offline --fasta /ref/hg38.fa \
        --symbol --canonical --af_gnomade --pick \
        --fields "Consequence,SYMBOL,gnomADe_AF" --fork 8

Then: python3 analysis_raas_gnomad.py
"""
import argparse
import gzip
import os
import re
import numpy as np
import pandas as pd
from scipy import stats

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"
LOG10_2 = 0.3010299957


def parse_vep_vcf(path):
    """Yield (swap_id, gnomADe_AF or nan) per record."""
    op = gzip.open if path.endswith(".gz") else open
    csq_fields = None
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                if "ID=CSQ" in line and "Format:" in line:
                    fmt = line.split("Format:")[1].strip().rstrip('">')
                    csq_fields = fmt.split("|")
                continue
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            swap_id, info = c[2], c[7]
            m = re.search(r"CSQ=([^;]+)", info)
            af = np.nan
            if m and csq_fields:
                # take max AF across annotations for this SNV
                for ann in m.group(1).split(","):
                    parts = ann.split("|")
                    d = dict(zip(csq_fields, parts))
                    v = d.get("gnomADe_AF", "")
                    if v not in ("", ".", None):
                        try:
                            af = np.nanmax([af, float(v)])
                        except ValueError:
                            pass
            yield swap_id, af


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=f"{B}/precursor_pairs.tsv")
    ap.add_argument("--vep", default=f"{B}/candidate_snvs.vep.vcf.gz")
    ap.add_argument("--out", default=B)
    a = ap.parse_args()

    if not os.path.exists(a.vep):
        raise SystemExit(f"VEP output not found: {a.vep}\n"
                         "run VEP on candidate_snvs.vcf first (see header).")

    # max gnomAD AF per swap id
    af_by_swap = {}
    for sid, af in parse_vep_vcf(a.vep):
        prev = af_by_swap.get(sid, np.nan)
        af_by_swap[sid] = np.nanmax([prev, af]) if not (np.isnan(prev) and np.isnan(af)) \
            else np.nan
    sk = pd.DataFrame(
        [(s, v) for s, v in af_by_swap.items()], columns=["swap_id", "gnomad_af"])
    parts = sk["swap_id"].str.split("|", expand=True)
    sk["acc"], sk["pos"], sk["wt"], sk["alt"] = parts[0], parts[1].astype(int), parts[2], parts[3]
    sk["has_snv"] = True
    print(f"swaps with >=1 candidate SNV: {len(sk):,} | "
          f"in gnomAD: {sk['gnomad_af'].notna().sum():,}")

    pairs = pd.read_csv(a.pairs, sep="\t")
    pairs = pairs[np.isfinite(pairs["raas_precursor"])]
    sw = pairs.groupby(["acc", "pos", "wt", "alt"], as_index=False).agg(
        raas_precursor=("raas_precursor", "mean"))
    sw["raas10"] = sw["raas_precursor"] * LOG10_2
    df = sw.merge(sk[["acc", "pos", "wt", "alt", "gnomad_af", "has_snv"]],
                  on=["acc", "pos", "wt", "alt"], how="left")
    df["has_snv"] = df["has_snv"].fillna(False)
    print(f"unique SAAPs: {len(df):,} | with a 1-nt SNV: {int(df['has_snv'].sum()):,} "
          f"| in gnomAD: {int(df['gnomad_af'].notna().sum()):,}")

    # bins: no 1-nt SNV (>=2 nt) | 1-nt SNV not in gnomAD | AF tertiles among observed
    obs = df[df["gnomad_af"].notna()].copy()
    df["af_bin"] = np.where(df["has_snv"], "not_in_gnomAD", "no_1nt_snv")
    if len(obs) >= 6:
        obs["af_bin"] = pd.qcut(obs["gnomad_af"], 3, duplicates="drop",
                                labels=["rare", "mid", "common"]).astype(str)
        df.loc[obs.index, "af_bin"] = obs["af_bin"].values
    order = ["no_1nt_snv", "not_in_gnomAD", "rare", "mid", "common"]
    order = [o for o in order if o in set(df["af_bin"])]

    g = df.groupby("af_bin")["raas10"].agg(["mean", "median", "size"]).reindex(order)
    print(g.to_string())
    g.reset_index().to_csv(f"{a.out}/raas_by_gnomad.tsv", sep="\t", index=False)
    if obs["gnomad_af"].notna().sum() >= 10:
        rho, p = stats.spearmanr(obs["gnomad_af"], obs["raas10"])
        print(f"Spearman(gnomAD AF, log10 RAAS | observed) = {rho:+.3f} (p={p:.2g})")
    ni = df.loc[df["af_bin"] == "not_in_gnomAD", "raas10"].dropna()
    inn = df.loc[df["af_bin"] != "not_in_gnomAD", "raas10"].dropna()
    if len(ni) and len(inn):
        u, pu = stats.mannwhitneyu(ni, inn)
        print(f"not-in-gnomAD vs in-gnomAD RAAS: MWU p={pu:.2g} "
              f"(medians {ni.median():.2f} vs {inn.median():.2f})")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        grp = [df.loc[df["af_bin"] == o, "raas10"].dropna().values for o in order]
        fig, ax = plt.subplots(figsize=(8, 5.5))
        bp = ax.boxplot(grp, tick_labels=order, showfliers=False,
                        patch_artist=True, medianprops=dict(color="black"))
        for patch in bp["boxes"]:
            patch.set_facecolor("#8a5a2a"); patch.set_alpha(.7)
        for i, x in enumerate(grp):
            ax.text(i + 1, ax.get_ylim()[1], f"n={len(x)}", ha="center",
                    va="top", fontsize=8, color="#555")
        ax.axhline(df["raas10"].median(), color="#c44", ls="--", lw=.8)
        ax.set_ylabel("log10 RAAS (SAAP/BP)")
        ax.set_title("#6  RAAS by gnomAD allele frequency")
        plt.tight_layout()
        plt.savefig(f"{a.out}/raas_by_gnomad.png", dpi=200, bbox_inches="tight")
        print(f"wrote {a.out}/raas_by_gnomad.png + raas_by_gnomad.tsv")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
