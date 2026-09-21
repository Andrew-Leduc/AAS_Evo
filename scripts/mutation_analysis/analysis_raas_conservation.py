#!/usr/bin/env python3
"""
#5  Average precursor RAAS by evolutionary conservation (phyloP) of the position.

Uses the codon's genomic coordinates (codon_map.tsv, from build_codon_map.py) to
look up per-base phyloP from the UCSC 100-way bigwig, averaged over the 3 codon
bases = residue-level conservation. Joins the same-run precursor RAAS, collapses
to one mean-RAAS value per unique SAAP, and bins RAAS by conservation.

Hypothesis: if detected SAAPs are real, more conserved (high phyloP, less
substitution-tolerant) positions might behave differently from unconstrained ones.

Setup (in the venv):
    pip install pyBigWig
    # phyloP bigwig (~9.2 GB), one-time:
    cd /scratch/leduc.an/AAS_Evo/SEQ_FILES
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"
DEFAULT_BW = "/scratch/leduc.an/AAS_Evo/SEQ_FILES/hg38.phyloP100way.bw"
LOG10_2 = 0.3010299957


def ucsc_chrom(c):
    c = str(c)
    if c == "MT":
        return "chrM"
    return c if c.startswith("chr") else "chr" + c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=f"{B}/precursor_pairs.tsv")
    ap.add_argument("--codons", default=f"{B}/codon_map.tsv")
    ap.add_argument("--bigwig", default=DEFAULT_BW)
    ap.add_argument("--out", default=B)
    a = ap.parse_args()

    if not os.path.exists(a.bigwig):
        raise SystemExit(f"phyloP bigwig not found: {a.bigwig}\n"
                         "download it (see header) or pass --bigwig.")
    import pyBigWig
    bw = pyBigWig.open(a.bigwig)
    have = set(bw.chroms().keys())

    cm = pd.read_csv(a.codons, sep="\t")
    cm = cm.dropna(subset=["g1", "g2", "g3"])
    for c in ("g1", "g2", "g3"):
        cm[c] = cm[c].astype(int)

    def site_phylop(chrom, g1, g2, g3):
        ch = ucsc_chrom(chrom)
        if ch not in have:
            return np.nan
        vals = []
        for g in (g1, g2, g3):
            try:
                v = bw.values(ch, g - 1, g)[0]      # 1-based -> 0-based
            except Exception:
                v = None
            if v is not None and not np.isnan(v):
                vals.append(v)
        return float(np.mean(vals)) if vals else np.nan

    cm["phylop"] = [site_phylop(c, a1, a2, a3) for c, a1, a2, a3
                    in zip(cm["chrom"], cm["g1"], cm["g2"], cm["g3"])]
    bw.close()
    site = cm.dropna(subset=["phylop"])[["acc", "pos", "wt", "phylop"]]
    print(f"sites with phyloP: {len(site):,} / {len(cm):,}")

    pairs = pd.read_csv(a.pairs, sep="\t")
    pairs = pairs.dropna(subset=["raas_precursor"])
    pairs = pairs[np.isfinite(pairs["raas_precursor"])]
    df = pairs.merge(site, on=["acc", "pos", "wt"], how="inner")
    # collapse to one mean-RAAS per unique SAAP
    sw = df.groupby(["acc", "pos", "wt", "alt"], as_index=False).agg(
        raas_precursor=("raas_precursor", "mean"), phylop=("phylop", "first"))
    sw["raas10"] = sw["raas_precursor"] * LOG10_2
    print(f"unique SAAPs with phyloP + RAAS: {len(sw):,}")

    rho, p = stats.spearmanr(sw["phylop"], sw["raas10"])
    print(f"Spearman(phyloP, log10 RAAS) = {rho:+.3f} (p={p:.2g})")

    sw["phylop_bin"] = pd.qcut(sw["phylop"], 5, duplicates="drop")
    g = sw.groupby("phylop_bin", observed=True)["raas10"].agg(["mean", "median", "size"])
    print(g.to_string())
    g.reset_index().astype({"phylop_bin": str}).to_csv(
        f"{a.out}/raas_by_phylop.tsv", sep="\t", index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cats = list(sw["phylop_bin"].cat.categories)
        grp = [sw.loc[sw["phylop_bin"] == c, "raas10"].dropna().values for c in cats]
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        bp = ax.boxplot(grp, tick_labels=[f"Q{i+1}" for i in range(len(cats))],
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black"))
        for patch in bp["boxes"]:
            patch.set_facecolor("#5a3a89"); patch.set_alpha(.7)
        for i, x in enumerate(grp):
            ax.text(i + 1, ax.get_ylim()[1], f"n={len(x)}", ha="center",
                    va="top", fontsize=8, color="#555")
        ax.axhline(sw["raas10"].median(), color="#c44", ls="--", lw=.8)
        ax.set_xlabel("phyloP conservation (Q1=low → Q5=high)")
        ax.set_ylabel("log10 RAAS (SAAP/BP)")
        ax.set_title(f"#5  RAAS by conservation (ρ={rho:+.3f}, p={p:.2g})")
        plt.tight_layout()
        plt.savefig(f"{a.out}/raas_by_phylop.png", dpi=200, bbox_inches="tight")
        print(f"wrote {a.out}/raas_by_phylop.png + raas_by_phylop.tsv")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
