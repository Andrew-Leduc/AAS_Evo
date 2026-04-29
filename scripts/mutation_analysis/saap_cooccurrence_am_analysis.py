#!/usr/bin/env python3
"""
AM score distribution: co-occurring missense vs background
==========================================================
For each TMT set, we have:
  - Co-occurring: missense mutations in genes that ALSO have a SAAP in that set
  - Background:   missense mutations in genes that do NOT have any SAAP in that set

Both groups are drawn from the same patients / same TMT sets, so the only
difference is whether the gene is SAAP-affected or not.

Outputs
-------
am_distribution_plot.pdf   — histogram + CDF comparison
am_stats.tsv               — summary statistics
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooccurrence",
                    default="/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_missense_cooccurrence.tsv")
    ap.add_argument("--summary",
                    default="/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_set_summary.tsv",
                    help="saap_set_summary.tsv — needed to know which (gene, TMT set) were tested")
    ap.add_argument("--pep-to-patient",
                    default=str(TSOUR_DIR / "peptide_to_patient.csv"),
                    help="Tsour peptide_to_patient.csv — for TMT set membership")
    ap.add_argument("--pep-to-protein",
                    default=str(TSOUR_DIR / "pep_to_protein.csv"),
                    help="Tsour pep_to_protein.csv — for SAAP gene mapping")
    ap.add_argument("--gdc-meta",
                    default=str(REPO_DIR / "metadata/GDC_meta/gdc_meta_matched.tsv"))
    ap.add_argument("--missense",
                    required=True,
                    help="all_missense_with_spurs.tsv")
    ap.add_argument("-o", "--out-dir",
                    default="/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence")
    return ap.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load co-occurrence results ────────────────────────────────────────────
    print("Loading co-occurrence table...", flush=True)
    cooc = pd.read_csv(args.cooccurrence, sep="\t")
    cooc["missense_am"] = pd.to_numeric(cooc["missense_am"], errors="coerce")
    print(f"  Co-occurrence records: {len(cooc):,}", flush=True)

    # Co-occurring AM scores — missense in same gene as a SAAP in same TMT set
    cooc_am = cooc["missense_am"].dropna()
    print(f"  Co-occurring AM scores: {len(cooc_am):,}", flush=True)

    # ── Build background: missense in same TMT sets but different genes ───────
    # Step 1: which (Dataset, TMT set) groups appear in any tested SAAP
    print("Loading summary table...", flush=True)
    summ = pd.read_csv(args.summary, sep="\t")
    tested_sets = set(zip(summ["Dataset"], summ["TMT set"]))
    # Gene × set combinations that have a SAAP — these are the co-occurring genes per set
    saap_gene_sets = set(zip(summ["gene"], summ["Dataset"], summ["TMT set"]))
    print(f"  Tested (Dataset, TMT set) groups: {len(tested_sets)}", flush=True)

    # Step 2: get all sample_ids in each TMT set via Tsour patient table
    print("Loading Tsour patient table...", flush=True)
    p2pat  = pd.read_csv(args.pep_to_patient, index_col=0)
    p2prot = pd.read_csv(args.pep_to_protein, index_col=0)
    p2prot = p2prot[["SAAP","name"]].drop_duplicates("SAAP").rename(columns={"name":"gene_symbol"})
    p2pat  = p2pat.merge(p2prot[["SAAP","gene_symbol"]], on="SAAP", how="left")
    p2pat["case_submitter_id"] = p2pat["Sample name"].str.rsplit("_", n=1).str[0]

    gdc = pd.read_csv(args.gdc_meta, sep="\t",
                      usecols=["case_submitter_id","gdc_file_id","sample_type"])
    gdc_map = (pd.concat([
        gdc[gdc["sample_type"]=="Primary Tumor"].drop_duplicates("case_submitter_id"),
        gdc.drop_duplicates("case_submitter_id")
    ]).drop_duplicates("case_submitter_id")
      .set_index("case_submitter_id")["gdc_file_id"].to_dict())
    p2pat["sample_id"] = p2pat["case_submitter_id"].map(gdc_map)

    # sample_ids per (Dataset, TMT set)
    set_to_uuids = (p2pat[p2pat["sample_id"].notna()]
                    .groupby(["Dataset","TMT set"])["sample_id"]
                    .apply(frozenset).to_dict())

    # All sample_ids that appear in any tested TMT set
    tested_uuids = set()
    for k, v in set_to_uuids.items():
        if k in tested_sets:
            tested_uuids |= set(v)
    print(f"  Unique sample UUIDs across tested sets: {len(tested_uuids):,}", flush=True)

    # Step 3: load missense for those samples only
    print("Loading missense table...", flush=True)
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["sample_id","SYMBOL","am_pathogenicity"])
    miss["am_pathogenicity"] = pd.to_numeric(miss["am_pathogenicity"], errors="coerce")
    miss = miss[miss["sample_id"].isin(tested_uuids) & miss["am_pathogenicity"].notna()]
    print(f"  Missense rows in tested samples with AM score: {len(miss):,}", flush=True)

    # Step 4: for each row, determine which (Dataset, TMT set) it belongs to
    # and whether its gene has a SAAP in that set
    # Build sample_id -> list of (Dataset, TMT set) it appears in
    uuid_to_sets = {}
    for (dataset, tmt_set), uuids in set_to_uuids.items():
        if (dataset, tmt_set) not in tested_sets:
            continue
        for uid in uuids:
            uuid_to_sets.setdefault(uid, []).append((dataset, tmt_set))

    print("Classifying missense as co-occurring vs background...", flush=True)
    is_cooc = []
    for row in miss[["sample_id","SYMBOL"]].itertuples(index=False):
        sets_for_sample = uuid_to_sets.get(row.sample_id, [])
        # co-occurring = this gene has a SAAP in at least one set this sample belongs to
        hit = any((row.SYMBOL, ds, ts) in saap_gene_sets for ds, ts in sets_for_sample)
        is_cooc.append(hit)

    miss["is_cooc"] = is_cooc
    background_am = miss[~miss["is_cooc"]]["am_pathogenicity"]
    cooc_am_recomputed = miss[miss["is_cooc"]]["am_pathogenicity"]
    print(f"  Co-occurring (same-set, SAAP-gene): {len(cooc_am_recomputed):,}", flush=True)
    print(f"  Background (same-set, non-SAAP gene): {len(background_am):,}", flush=True)

    # ── Statistics ────────────────────────────────────────────────────────────
    u, p_mw = stats.mannwhitneyu(cooc_am_recomputed, background_am, alternative="two-sided")
    ks, p_ks = stats.ks_2samp(cooc_am_recomputed, background_am)

    print(f"\n── AM SCORE STATISTICS ──────────────────────────────────────────────")
    print(f"Co-occurring:  n={len(cooc_am_recomputed):,}  "
          f"mean={cooc_am_recomputed.mean():.3f}  median={cooc_am_recomputed.median():.3f}")
    print(f"Background:    n={len(background_am):,}  "
          f"mean={background_am.mean():.3f}  median={background_am.median():.3f}")
    print(f"Mann-Whitney U (two-sided): p={p_mw:.2e}")
    print(f"KS test:                   D={ks:.3f}  p={p_ks:.2e}")

    bins_labels = [(0.0, 0.34, "likely_benign"),
                   (0.34, 0.564, "ambiguous"),
                   (0.564, 1.01, "likely_pathogenic")]
    print(f"\n{'Category':<22} {'Co-occurring':>14} {'Background':>14}")
    for lo, hi, label in bins_labels:
        n_c = ((cooc_am_recomputed >= lo) & (cooc_am_recomputed < hi)).sum()
        n_b = ((background_am >= lo) & (background_am < hi)).sum()
        pct_c = 100 * n_c / len(cooc_am_recomputed)
        pct_b = 100 * n_b / len(background_am)
        print(f"{label:<22} {n_c:>6} ({pct_c:4.1f}%)   {n_b:>7} ({pct_b:4.1f}%)")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    bins = np.linspace(0, 1, 41)
    ax = axes[0]
    ax.hist(background_am, bins=bins, density=True, alpha=0.5, color="#888888",
            label=f"Background\n(same patients, non-SAAP genes)\nn={len(background_am):,}  "
                  f"median={background_am.median():.3f}")
    ax.hist(cooc_am_recomputed, bins=bins, density=True, alpha=0.7, color="#d62728",
            label=f"Co-occurring\n(same patients, SAAP-affected gene)\nn={len(cooc_am_recomputed):,}  "
                  f"median={cooc_am_recomputed.median():.3f}")
    for threshold, ls in [(0.34, "--"), (0.564, "-.")]:
        ax.axvline(threshold, color="k", lw=0.8, ls=ls, alpha=0.5)
    ax.set_xlabel("AlphaMissense pathogenicity score", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(f"AM score distribution\nMann-Whitney p={p_mw:.2e} | KS D={ks:.3f} p={p_ks:.2e}",
                 fontsize=9)
    ax.legend(fontsize=8)

    ax2 = axes[1]
    for data, color, label in [
        (background_am,       "#888888", f"Background (n={len(background_am):,})"),
        (cooc_am_recomputed,  "#d62728", f"Co-occurring (n={len(cooc_am_recomputed):,})"),
    ]:
        xs = np.sort(data.values)
        ys = np.arange(1, len(xs)+1) / len(xs)
        ax2.plot(xs, ys, color=color, lw=1.5, label=label)
    for threshold, ls in [(0.34, "--"), (0.564, "-.")]:
        ax2.axvline(threshold, color="k", lw=0.8, ls=ls, alpha=0.5)
    ax2.set_xlabel("AlphaMissense pathogenicity score", fontsize=10)
    ax2.set_ylabel("Cumulative fraction", fontsize=10)
    ax2.set_title("Empirical CDF", fontsize=9)
    ax2.legend(fontsize=8)

    plt.suptitle("AlphaMissense scores: co-occurring missense vs background\n"
                 "(same TMT-set patients; co-occurring = SAAP-affected gene, "
                 "background = all other genes)",
                 fontsize=9)
    plt.tight_layout()
    plt.savefig(out_dir / "am_distribution_plot.pdf", dpi=150)
    plt.savefig(out_dir / "am_distribution_plot.png", dpi=150)
    print(f"\nPlot saved to {out_dir}/am_distribution_plot.pdf", flush=True)

    stat_rows = []
    for lo, hi, label in bins_labels:
        n_c = ((cooc_am_recomputed >= lo) & (cooc_am_recomputed < hi)).sum()
        n_b = ((background_am >= lo) & (background_am < hi)).sum()
        stat_rows.append({
            "am_class": label,
            "cooc_n": n_c, "cooc_pct": round(100*n_c/len(cooc_am_recomputed), 2),
            "bg_n": n_b,   "bg_pct":   round(100*n_b/len(background_am), 2),
        })
    pd.DataFrame(stat_rows).to_csv(out_dir / "am_stats.tsv", sep="\t", index=False)
    print(f"Stats saved to {out_dir}/am_stats.tsv", flush=True)


if __name__ == "__main__":
    main()
