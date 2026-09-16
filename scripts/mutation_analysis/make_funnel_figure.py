#!/usr/bin/env python3
"""
Processing-funnel figure: genomics -> contact-proximal AAS search space.

Edit GENOMICS and STAGES below (rough numbers are pre-filled), then run:
    python3 make_funnel_figure.py
-> writes aas_processing_funnel.pdf in the current directory.
Self-contained (numbers are hard-coded), so it runs anywhere.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "aas_processing_funnel.pdf"

# ── EDIT: genomics input context (header band) ───────────────────────────────
GENOMICS = ("~1,570 tumor/normal WXS genomes   ·   ~1,057-1,192 patients"
            "   ·   ~1,672 tissue BAMs")

# ── EDIT: funnel stages — (label, count, note). Keep one unit (missense counts) ─
STAGES = [
    ("Total missense calls",      15_000_000, "bcftools + VEP, all samples   [CONFIRM]"),
    ("Unique missense variants",  305_437,    "distinct gene+pos+substitution (VAF>=0.3)"),
    ("Eligible missense",         61_417,     ">=5 patients & <75% of cohort"),
    ("With a long-range contact", 20_761,     ">=21 aa apart & <5 A (AlphaFold)"),
    ("Seeded the AAS search",     19_012,     "-> swap peptides, 7,636 genes"),
]

labels = [s[0] for s in STAGES]
counts = [s[1] for s in STAGES]
notes = [s[2] for s in STAGES]

# widths scaled by log10(count) so 3 orders of magnitude stay visible
lg = np.array([np.log10(max(c, 1)) for c in counts])
w = 0.28 + 0.72 * (lg - lg.min()) / (lg.max() - lg.min() + 1e-9)

fig, ax = plt.subplots(figsize=(10, 6.8))
colors = plt.cm.viridis(np.linspace(0.12, 0.78, len(STAGES)))
bh, gap = 0.65, 0.42
for i in range(len(STAGES)):
    y = -i * (bh + gap)
    if i < len(STAGES) - 1:                       # connector trapezoid to next stage
        yn = -(i + 1) * (bh + gap)
        ax.fill([-w[i] / 2, w[i] / 2, w[i + 1] / 2, -w[i + 1] / 2],
                [y - bh / 2, y - bh / 2, yn + bh / 2, yn + bh / 2],
                color="#e8e8e8", zorder=1)
    ax.add_patch(plt.Rectangle((-w[i] / 2, y - bh / 2), w[i], bh,
                               color=colors[i], ec="white", lw=1.5, zorder=2))
    ax.text(0, y + 0.09, f"{counts[i]:,}", ha="center", va="center",
            color="white", fontweight="bold", fontsize=13, zorder=3)
    ax.text(0, y - 0.17, notes[i], ha="center", va="center",
            color="white", fontsize=7.5, zorder=3)
    ax.text(-0.60, y, labels[i], ha="right", va="center", fontsize=10.5, fontweight="bold")
    if i > 0 and counts[i - 1]:
        ax.text(0.60, y, f"{100 * counts[i] / counts[i - 1]:.1f}% of prev",
                ha="left", va="center", fontsize=8, color="#666")

ax.text(0, bh / 2 + 0.62, "GENOMICS PROCESSED", ha="center", fontsize=9.5,
        fontweight="bold", color="#333")
ax.text(0, bh / 2 + 0.34, GENOMICS, ha="center", fontsize=9, color="#333")

ax.set_xlim(-1.55, 1.25)
ax.set_ylim(-len(STAGES) * (bh + gap) + 0.15, bh / 2 + 1.0)
ax.axis("off")
ax.set_title("AAS_Evo: genomics → contact-proximal AAS search funnel",
             fontsize=13, fontweight="bold", pad=12)
plt.tight_layout()
plt.savefig(OUT, bbox_inches="tight", dpi=200)
print(f"wrote {OUT}")
