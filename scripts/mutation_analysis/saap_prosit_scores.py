#!/usr/bin/env python3
"""
Before/after-Prosit (MSBooster) scores for the filtered SAAP list.

FragPipe's MSBooster writes an `_edited.pin` per run holding BOTH:
  - the original search score  (BEFORE Prosit):  hyperscore  (higher = better)
  - the deep-learning features (the Prosit evidence, AFTER):
        delta_RT_loess          predicted-RT error (|.| smaller = better)
        spectral similarity     (e.g. brayCurtis / spectral_entropy / cosine;
                                 higher = predicted spectrum matches observed)

For each SAAP PSM (protein accession = the {ACC}-{SWAP}-{HASH} mock entry) we pull
these and compare against a baseline of ordinary target peptides (canonical
proteins). If the swaps are real, their spectral/RT agreement should look like the
baseline; if they're misidentifications, the predicted spectrum/RT won't match and
those features will be worse than baseline even when the raw hyperscore passed.

Streams the pins (regex-prefilter) so it's fast. Auto-detects column names and
prints what it used. Applies the same PTM-mass clean filter as the rest.

Run (in the venv):  python3 saap_prosit_scores.py
"""
import glob
import os
import random
import re
import numpy as np
import pandas as pd
from scipy import stats

RESULTS_BASE = "/scratch/leduc.an/AAS_Evo/MS_SEARCH/results_contact"
OUT_DIR = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"
SWAP_RE = re.compile(r'\b([A-Z0-9]+)-([A-Z]\d+[A-Z])-[0-9A-F]{4}\b')
BASELINE_KEEP = 0.02          # subsample fraction for baseline target PSMs
random.seed(0)

_AA_MASS = {'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,'F':147.06841,
            'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,
            'M':131.04049,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,
            'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06333}
_PTM = [14.01565,15.99491,2.01565,0.98402,42.01057,79.96633,18.01056,
        14.99967,28.03130,16.97893,31.98983,1.96804]
_SUSP = {(w, a) for w in _AA_MASS for a in _AA_MASS
         if w != a and any(abs(abs(_AA_MASS[a]-_AA_MASS[w])-p) < 0.05 for p in _PTM)}

# candidate column names (lowercased)
COL_BEFORE = ["hyperscore"]
COL_RT = ["delta_rt_loess", "abs_rt_diff", "rt_diff", "braycurtis_rt"]
COL_SPEC = ["braycurtis", "spectral_entropy_similarity", "unweighted_spectral_entropy",
            "cosine_similarity", "cosine", "pearson_corr", "spectral_angle",
            "dot_product", "entropy"]


def find_edited_pins():
    pins = glob.glob(os.path.join(RESULTS_BASE, "*", "*_1", "*.pin")) \
        or glob.glob(os.path.join(RESULTS_BASE, "*", "*.pin"))
    return sorted(pins)


def pick_cols(header):
    low = {c.lower(): c for c in header}
    def first(cands):
        for c in cands:
            if c in low:
                return low[c]
        return None
    return (first(COL_BEFORE), first(COL_RT),
            [low[c] for c in low if any(s in c for s in COL_SPEC)][:1])


def main():
    pins = find_edited_pins()
    if not pins:
        raise SystemExit(f"no .pin files under {RESULTS_BASE} — MSBooster pins may "
                         "have been cleaned; check FragPipe output retention.")
    print(f"{len(pins)} pin files")

    # detect columns from the first pin
    with open(pins[0]) as f:
        header = f.readline().rstrip("\n").split("\t")
    c_before, c_rt, c_spec = pick_cols(header)
    c_spec = c_spec[0] if c_spec else None
    idx = {c: i for i, c in enumerate(header)}
    print(f"detected columns:\n  before (hyperscore): {c_before}\n"
          f"  RT feature         : {c_rt}\n  spectral feature   : {c_spec}")
    if not any([c_before, c_rt, c_spec]):
        print("  (none of the expected feature columns found — dumping header:)")
        print("  " + " | ".join(header))
        raise SystemExit("no usable feature columns; paste the header above.")

    pep_i = idx.get("Peptide", len(header) - 2)
    want = {"before": idx.get(c_before), "rt": idx.get(c_rt), "spec": idx.get(c_spec)}
    rows = {"saap": {k: [] for k in want}, "base": {k: [] for k in want}}

    for pi, pf in enumerate(pins):
        if pi % 20 == 0:
            print(f"  {pi}/{len(pins)}", flush=True)
        with open(pf) as f:
            next(f, None)
            for line in f:
                is_swap = bool(SWAP_RE.search(line))
                if not is_swap and random.random() > BASELINE_KEEP:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= pep_i:
                    continue
                grp = None
                if is_swap:
                    m = SWAP_RE.search(line)
                    sw = m.group(2)
                    mm = re.match(r'^([A-Z])\d+([A-Z])$', sw)
                    if mm and (mm.group(1), mm.group(2)) in _SUSP:
                        continue                       # drop PTM-mass confounds
                    grp = "saap"
                else:
                    grp = "base"
                for k, ci in want.items():
                    if ci is not None and ci < len(parts):
                        try:
                            rows[grp][k].append(float(parts[ci]))
                        except ValueError:
                            pass

    os.makedirs(OUT_DIR, exist_ok=True)
    summ = []
    print("\nfeature            group   n        median      mean")
    for k, name in [("before", f"{c_before} (before)"),
                    ("rt", f"{c_rt} (RT err)"),
                    ("spec", f"{c_spec} (spectral sim)")]:
        if want[k] is None:
            continue
        s = np.array(rows["saap"][k]); b = np.array(rows["base"][k])
        s = s[np.isfinite(s)]; b = b[np.isfinite(b)]
        if len(s) and len(b):
            u, p = stats.mannwhitneyu(s, b)
            print(f"{name:26s} SAAP  {len(s):7,d}  {np.median(s):9.3f}  {np.mean(s):9.3f}")
            print(f"{'':26s} base  {len(b):7,d}  {np.median(b):9.3f}  {np.mean(b):9.3f}   MWU p={p:.1e}")
            summ.append(dict(feature=name, saap_n=len(s), saap_median=np.median(s),
                             base_n=len(b), base_median=np.median(b), mwu_p=p))
    pd.DataFrame(summ).to_csv(f"{OUT_DIR}/saap_prosit_scores.tsv", sep="\t", index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        feats = [(k, nm) for k, nm in
                 [("before", c_before), ("rt", c_rt), ("spec", c_spec)] if want[k] is not None]
        fig, ax = plt.subplots(1, len(feats), figsize=(5.2 * len(feats), 5))
        if len(feats) == 1:
            ax = [ax]
        for j, (k, nm) in enumerate(feats):
            s = np.array(rows["saap"][k]); b = np.array(rows["base"][k])
            s = s[np.isfinite(s)]; b = b[np.isfinite(b)]
            bp = ax[j].boxplot([b, s], tick_labels=["baseline", "SAAP"],
                               showfliers=False, patch_artist=True,
                               medianprops=dict(color="black"))
            bp["boxes"][0].set_facecolor("#7aa"); bp["boxes"][1].set_facecolor("#c86")
            ax[j].set_title(nm); ax[j].set_ylabel(nm)
        fig.suptitle("SAAP vs baseline: before-search & Prosit/MSBooster features",
                     fontweight="bold")
        plt.tight_layout()
        plt.savefig(f"{OUT_DIR}/saap_prosit_scores.png", dpi=200, bbox_inches="tight")
        print(f"\nwrote {OUT_DIR}/saap_prosit_scores.png + saap_prosit_scores.tsv")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
