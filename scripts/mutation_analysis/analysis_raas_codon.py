#!/usr/bin/env python3
"""
#3  Average precursor RAAS by # nucleotide mismatches for the translation error
#4  Average precursor RAAS by WT-codon usage frequency

Joins the same-run precursor RAAS (contact_saap/precursor_pairs.tsv) with the
reference codon at each site (codon_map.tsv from build_codon_map.py).

#3: min # nt changes from the ACTUAL wt codon to the nearest codon encoding the
    observed alt aa. Near-cognate mistranslation is dominated by 1-mismatch, so
    if these are real translation errors 1-mismatch swaps should carry the signal
    (higher/more RAAS) and >=2-mismatch swaps should be depleted.
#4: usage frequency (per-1000, human) of the wt codon being decoded. Rare codons
    are decoded more slowly / error-prone -> hypothesis: higher RAAS.

RAAS = log2(SAAP precursor intensity / BP precursor intensity), per same-run pair.
"""
import argparse
import os
import numpy as np
import pandas as pd

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
for _i, _a in enumerate(_AAS):
    CODON_TABLE[_BASES[_i >> 4] + _BASES[(_i >> 2) & 3] + _BASES[_i & 3]] = _a
AA2CODONS = {}
for _c, _a in CODON_TABLE.items():
    AA2CODONS.setdefault(_a, []).append(_c)

# human codon usage, per-1000 (Kazusa Homo sapiens)
CODON_FREQ = {
    'TTT':17.6,'TTC':20.3,'TTA':7.7,'TTG':12.9,'CTT':13.2,'CTC':19.6,'CTA':7.2,'CTG':39.6,
    'ATT':16.0,'ATC':20.8,'ATA':7.5,'ATG':22.0,'GTT':11.0,'GTC':14.5,'GTA':7.1,'GTG':28.1,
    'TCT':15.2,'TCC':17.7,'TCA':12.2,'TCG':4.4,'CCT':17.5,'CCC':19.8,'CCA':16.9,'CCG':6.9,
    'ACT':13.1,'ACC':18.9,'ACA':15.1,'ACG':6.1,'GCT':18.4,'GCC':27.7,'GCA':15.8,'GCG':7.4,
    'TAT':12.2,'TAC':15.3,'TAA':1.0,'TAG':0.8,'CAT':10.9,'CAC':15.1,'CAA':12.3,'CAG':34.2,
    'AAT':17.0,'AAC':19.1,'AAA':24.4,'AAG':31.9,'GAT':21.8,'GAC':25.1,'GAA':29.0,'GAG':39.6,
    'TGT':10.6,'TGC':12.6,'TGA':1.6,'TGG':13.2,'CGT':4.5,'CGC':10.4,'CGA':6.2,'CGG':11.4,
    'AGT':12.1,'AGC':19.5,'AGA':12.2,'AGG':12.0,'GGT':10.8,'GGC':22.2,'GGA':16.5,'GGG':16.5,
}


def min_mismatch(wt_codon, alt_aa):
    cands = AA2CODONS.get(alt_aa, [])
    if not cands:
        return np.nan
    return min(sum(x != y for x, y in zip(wt_codon, c)) for c in cands)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=f"{B}/precursor_pairs.tsv")
    ap.add_argument("--codons", default=f"{B}/codon_map.tsv")
    ap.add_argument("--out", default=B)
    a = ap.parse_args()

    pairs = pd.read_csv(a.pairs, sep="\t")
    cm = pd.read_csv(a.codons, sep="\t", usecols=["acc", "pos", "wt", "codon"])
    df = pairs.merge(cm, on=["acc", "pos", "wt"], how="inner")
    df = df.dropna(subset=["raas_precursor", "codon"])
    df = df[np.isfinite(df["raas_precursor"])]          # drop +/-inf from zero-intensity pairs
    df = df[df["codon"].str.len() == 3]
    print(f"pairs with codon + RAAS: {len(df):,} "
          f"({df.groupby(['acc','pos','wt','alt']).ngroups:,} unique SAAPs)")

    df["n_mismatch"] = [min_mismatch(c, al) for c, al in zip(df["codon"], df["alt"])]
    df["wt_codon_freq"] = df["codon"].map(CODON_FREQ)

    # collapse to ONE value per unique SAAP (mean RAAS across its runs) so swaps
    # seen in many runs are not pseudo-replicated.
    df = df.groupby(["acc", "pos", "wt", "alt"], as_index=False).agg(
        raas_precursor=("raas_precursor", "mean"),
        n_mismatch=("n_mismatch", "first"),
        wt_codon_freq=("wt_codon_freq", "first"),
        codon=("codon", "first"))
    print(f"unique SAAPs (mean RAAS per swap): {len(df):,}")

    # RAAS is stored as log2(SAAP/BP); convert to log10 for the plots/tables.
    LOG10_2 = 0.3010299957
    df["raas10"] = df["raas_precursor"] * LOG10_2

    # ── #3 RAAS by # nt mismatch ──
    print("\n#3  log10 RAAS by # nucleotide mismatches")
    g3 = df.groupby("n_mismatch")["raas10"].agg(["mean", "median", "size"])
    print(g3.to_string())
    g3.reset_index().to_csv(f"{a.out}/raas_by_nt_mismatch.tsv", sep="\t", index=False)
    from scipy import stats
    grp3 = [df.loc[df["n_mismatch"] == k, "raas10"].dropna().values for k in (1, 2, 3)]
    if all(len(x) for x in grp3):
        h, ph = stats.kruskal(*grp3)
        print(f"Kruskal-Wallis across mismatch classes: H={h:.2f} (p={ph:.2g})")

    # ── #4 RAAS by wt codon frequency (quintiles) ──
    print("\n#4  log10 RAAS by WT codon usage frequency (quintiles)")
    d4 = df.dropna(subset=["wt_codon_freq"]).copy()
    d4["freq_bin"] = pd.qcut(d4["wt_codon_freq"], 5, duplicates="drop")
    g4 = d4.groupby("freq_bin", observed=True)["raas10"].agg(["mean", "median", "size"])
    print(g4.to_string())
    g4.reset_index().astype({"freq_bin": str}).to_csv(
        f"{a.out}/raas_by_codon_freq.tsv", sep="\t", index=False)
    rho, p = stats.spearmanr(d4["wt_codon_freq"], d4["raas10"])
    print(f"Spearman(wt_codon_freq, log10 RAAS) = {rho:+.3f} (p={p:.1e})")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        overall = df["raas10"].median()
        fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))

        # #3 boxplot by nt-mismatch
        bp = ax[0].boxplot(grp3, tick_labels=["1", "2", "3"], showfliers=False,
                           patch_artist=True, medianprops=dict(color="black"))
        for patch in bp["boxes"]:
            patch.set_facecolor("#3a8489"); patch.set_alpha(.7)
        for i, x in enumerate(grp3):
            ax[0].text(i + 1, ax[0].get_ylim()[1], f"n={len(x)}", ha="center",
                       va="top", fontsize=8, color="#555")
        ax[0].axhline(overall, color="#c44", ls="--", lw=.8)
        ax[0].set_xlabel("# nt mismatches (wt codon → alt)")
        ax[0].set_ylabel("log10 RAAS (SAAP/BP)")
        ph_txt = f"KW p={ph:.2g}" if all(len(x) for x in grp3) else ""
        ax[0].set_title(f"#3  RAAS by codon mismatch   {ph_txt}")

        # #4 boxplot by codon-freq quintile
        cats = list(d4["freq_bin"].cat.categories)
        grp4 = [d4.loc[d4["freq_bin"] == c, "raas10"].dropna().values for c in cats]
        bp4 = ax[1].boxplot(grp4, tick_labels=[f"Q{i+1}" for i in range(len(cats))],
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black"))
        for patch in bp4["boxes"]:
            patch.set_facecolor("#33627f"); patch.set_alpha(.7)
        for i, x in enumerate(grp4):
            ax[1].text(i + 1, ax[1].get_ylim()[1], f"n={len(x)}", ha="center",
                       va="top", fontsize=8, color="#555")
        ax[1].axhline(overall, color="#c44", ls="--", lw=.8)
        ax[1].set_xlabel("WT codon usage (Q1=rare → Q5=common)")
        ax[1].set_ylabel("log10 RAAS (SAAP/BP)")
        ax[1].set_title(f"#4  RAAS by codon freq (ρ={rho:+.3f}, p={p:.2g})")
        plt.tight_layout()
        plt.savefig(f"{a.out}/raas_by_codon.png", dpi=200, bbox_inches="tight")
        print(f"\nwrote {a.out}/raas_by_codon.png + raas_by_nt_mismatch.tsv + raas_by_codon_freq.tsv")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
