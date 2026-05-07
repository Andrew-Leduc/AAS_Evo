#!/usr/bin/env python3
"""
Violin plot of RAAS difference (with missense - without missense) for
structurally contacted SAAP-missense pairs, split by AM/SPURS concordance group.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path

REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

AM_THRESHOLD    = 0.564
SPURS_THRESHOLD = 0.5
AM_BENIGN_MAX   = 0.1
CONTACT_ANG     = 8.0

GDC_META    = str(REPO_DIR / "metadata" / "GDC_meta" / "gdc_meta_matched.tsv")
COOC        = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_contact_annotated.tsv"
OUT_DEFAULT = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence"

GROUP_ORDER  = ["Both AM+SPURS", "SPURS only", "AM only", "Benign"]
GROUP_COLORS = {"Both AM+SPURS": "#d62728", "SPURS only": "#9467bd",
                "AM only": "#ff7f0e", "Benign": "#2ca02c"}


def classify(am, spurs):
    am_destab    = pd.notna(am)    and am    >= AM_THRESHOLD
    spurs_destab = pd.notna(spurs) and spurs >= SPURS_THRESHOLD
    am_benign    = pd.notna(am)    and am    <= AM_BENIGN_MAX
    if am_destab and spurs_destab:
        return "Both AM+SPURS"
    if spurs_destab and not am_destab:
        return "SPURS only"
    if am_destab and not spurs_destab:
        return "AM only"
    if am_benign:
        return "Benign"
    return None  # ambiguous — exclude


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooc",     default=COOC)
    ap.add_argument("--gdc-meta", default=GDC_META)
    ap.add_argument("-o",         default=OUT_DEFAULT)
    args = ap.parse_args()

    # ── Load data ─────────────────────────────────────────────────────────────
    print("Loading tables...", flush=True)
    p2pat = pd.read_csv(TSOUR_DIR / "peptide_to_patient.csv",
                        usecols=["SAAP","Sample name","RAAS"])
    p2pat["RAAS"] = pd.to_numeric(p2pat["RAAS"], errors="coerce")
    p2pat = p2pat[np.isfinite(p2pat["RAAS"])]
    p2pat["case_id"] = p2pat["Sample name"].str.rsplit("_", n=1).str[0]

    gdc = pd.read_csv(args.gdc_meta, sep="\t",
                      usecols=["gdc_file_id","case_submitter_id"])
    uuid_to_case = gdc.set_index("gdc_file_id")["case_submitter_id"].to_dict()

    df = pd.read_csv(args.cooc, sep="\t")
    df["missense_am"]        = pd.to_numeric(df["missense_am"],        errors="coerce")
    df["missense_spurs_ddg"] = pd.to_numeric(df["missense_spurs_ddg"], errors="coerce")

    # Filter: in contact, not same-site
    contacts = df[(df["in_structural_contact"] == True) &
                  (df["structural_distance_ang"] > 0.0)].copy()
    contacts["group"] = contacts.apply(
        lambda r: classify(r["missense_am"], r["missense_spurs_ddg"]), axis=1)
    contacts = contacts[contacts["group"].notna()]
    contacts["missense_case_id"] = contacts["missense_sample_id"].map(uuid_to_case)

    print(f"Contacted pairs (excl. same-site): {len(contacts):,}", flush=True)
    for g in GROUP_ORDER:
        print(f"  {g}: {(contacts['group']==g).sum()}", flush=True)

    # ── Compute RAAS diff per (SAAP, missense_swap, group) ────────────────────
    group_diffs = {g: [] for g in GROUP_ORDER}

    for (saap, miss_swap, group), grp in contacts.groupby(["SAAP","missense_swap","group"]):
        missense_cases = set(grp["missense_case_id"].dropna())
        saap_pat = p2pat[p2pat["SAAP"] == saap].copy()
        saap_pat["has_missense"] = saap_pat["case_id"].isin(missense_cases)
        with_m    = saap_pat[saap_pat["has_missense"]]["RAAS"].dropna()
        without_m = saap_pat[~saap_pat["has_missense"]]["RAAS"].dropna()
        if len(with_m) == 0 or len(without_m) == 0:
            continue
        diff = with_m.mean() - without_m.mean()
        group_diffs[group].append(diff)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Panel A: violin
    ax = axes[0]
    data_for_violin = [group_diffs[g] for g in GROUP_ORDER if len(group_diffs[g]) >= 2]
    labels_violin   = [g for g in GROUP_ORDER if len(group_diffs[g]) >= 2]

    if data_for_violin:
        parts = ax.violinplot(data_for_violin, positions=range(len(data_for_violin)),
                              showmedians=True, showextrema=False)
        for pc, g in zip(parts["bodies"], labels_violin):
            pc.set_facecolor(GROUP_COLORS[g]); pc.set_alpha(0.7)
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(2)

        # overlay individual points
        for i, (data, g) in enumerate(zip(data_for_violin, labels_violin)):
            ax.scatter([i] * len(data), data, color=GROUP_COLORS[g],
                       alpha=0.5, s=25, zorder=3)

    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xticks(range(len(labels_violin)))
    ax.set_xticklabels(labels_violin, fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("RAAS diff (with missense − without missense)", fontsize=9)
    ax.set_title("RAAS enrichment by AM/SPURS concordance\n(structurally contacted pairs, excl. same-site)", fontsize=9)

    # Panel B: mean ± SEM bar chart + stats
    ax2 = axes[1]
    xs, means, sems, ns = [], [], [], []
    for i, g in enumerate(GROUP_ORDER):
        d = group_diffs[g]
        if len(d) == 0:
            continue
        xs.append(i); means.append(np.mean(d))
        sems.append(stats.sem(d) if len(d) > 1 else 0)
        ns.append(len(d))

    bars = ax2.bar(xs, means, yerr=sems, capsize=5,
                   color=[GROUP_COLORS[GROUP_ORDER[x]] for x in xs],
                   alpha=0.8, error_kw={"linewidth": 1.5})
    for x, m, n in zip(xs, means, ns):
        t_val, p_val = stats.ttest_1samp(group_diffs[GROUP_ORDER[x]], popmean=0) \
                       if n >= 3 else (np.nan, np.nan)
        label = f"n={n}\np={p_val:.3f}" if not np.isnan(p_val) else f"n={n}"
        ax2.text(x, (m or 0) + (sems[xs.index(x)] or 0) + 0.05,
                 label, ha="center", va="bottom", fontsize=7)

    ax2.axhline(0, color="k", lw=0.8, ls="--")
    ax2.set_xticks(xs)
    ax2.set_xticklabels([GROUP_ORDER[x] for x in xs], fontsize=8,
                        rotation=15, ha="right")
    ax2.set_ylabel("mean RAAS diff", fontsize=9)
    ax2.set_title("Mean RAAS enrichment ± SEM\n(t-test vs 0)", fontsize=9)

    plt.suptitle("SAAP RAAS enrichment in patients with structurally proximal missense\n"
                 "(positive = higher mistranslation in missense carriers)", fontsize=10)
    plt.tight_layout()

    out = Path(args.o) / "saap_contact_raas_violin.pdf"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved {out}", flush=True)

    # Print summary
    print("\nSummary:")
    for g in GROUP_ORDER:
        d = group_diffs[g]
        if not d: continue
        t, p = stats.ttest_1samp(d, popmean=0) if len(d) >= 3 else (np.nan, np.nan)
        print(f"  {g:<18}: n={len(d):>3}  mean={np.mean(d):>7.3f}  "
              f"median={np.median(d):>7.3f}  p={p:.3f}" if not np.isnan(p)
              else f"  {g:<18}: n={len(d):>3}  mean={np.mean(d):>7.3f}")


if __name__ == "__main__":
    main()
