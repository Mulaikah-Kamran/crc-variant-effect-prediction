"""
Generate the four core figures from the final, audited scored dataset.
Run after src/05_primary_analysis.py and src/06_hotspot_concordance.py.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

TIER_A = ["ALDOB", "FBXL16", "IL1RN", "METTL27", "MSLN", "RAC3", "SLC38A11", "WNT11"]
TIER_B = ["TP53", "KRAS", "PIK3CA", "FBXW7", "BRAF", "SMAD4"]
TIER_A_COLOR = "#4C72B0"
TIER_B_COLOR = "#C44E52"


def fig1_pipeline_overview():
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.axis("off")

    stages = [
        "TCGA COADREAD\n(colon subset)\nsomatic mutations",
        "Filter to\nmissense,\nQC + WT check",
        "UniProt\ncanonical\nsequences",
        "ESM2\nmasked-marginal\nscoring",
        "Per-gene\nscore\ndistributions",
    ]
    n = len(stages)
    box_w, box_h = 1.7, 1.1
    xs = np.linspace(0.3, 0.3 + (n - 1) * 2.1, n)
    y = 0.5

    for x, label in zip(xs, stages):
        box = FancyBboxPatch((x, y - box_h / 2), box_w, box_h,
                              boxstyle="round,pad=0.02,rounding_size=0.08",
                              linewidth=1.2, edgecolor="#333333", facecolor="#F2F2F2")
        ax.add_patch(box)
        ax.text(x + box_w / 2, y, label, ha="center", va="center", fontsize=9)

    for i in range(n - 1):
        x_start = xs[i] + box_w
        x_end = xs[i + 1]
        arrow = FancyArrowPatch((x_start, y), (x_end, y),
                                 arrowstyle="-|>", mutation_scale=14,
                                 color="#555555", linewidth=1.2)
        ax.add_patch(arrow)

    ax.set_xlim(0, xs[-1] + box_w + 0.3)
    ax.set_ylim(-0.3, 1.3)
    ax.text(xs[0], 1.15, "799 raw records \u2192 673 missense \u2192 259 unique variants",
            fontsize=8.5, color="#555555", ha="left")
    plt.tight_layout()
    plt.savefig("results/figures/01_pipeline_overview.png", dpi=200)
    plt.close()


def fig2_per_gene_distributions():
    with open("data/scored/scored_variants.json") as f:
        scored = json.load(f)

    from collections import defaultdict
    by_gene = defaultdict(list)
    for r in scored:
        by_gene[r["gene"]].append(r["esm2_score"])

    genes_ordered = TIER_A + TIER_B
    data = [by_gene.get(g, []) for g in genes_ordered]
    colors = [TIER_A_COLOR] * len(TIER_A) + [TIER_B_COLOR] * len(TIER_B)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    positions = range(1, len(genes_ordered) + 1)
    bp = ax.boxplot(data, positions=positions, widths=0.55, showfliers=False,
                     patch_artist=True, medianprops=dict(color="black", linewidth=1.5))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor(color)

    # overlay individual points
    rng = np.random.default_rng(0)
    for pos, vals, color in zip(positions, data, colors):
        if not vals:
            continue
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter([pos + j for j in jitter], vals, s=14, color=color, alpha=0.7, zorder=3)

    ax.axhline(0, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_xticks(positions)
    ax.set_xticklabels(genes_ordered, rotation=40, ha="right")
    ax.set_ylabel("ESM2 masked-marginal score\n(log P(mutant) \u2212 log P(WT))")
    ax.set_title("Per-gene distribution of scored missense variants (n=259 unique variants)", fontsize=11)

    handles = [mpatches.Patch(color=TIER_A_COLOR, alpha=0.5, label="Tier A \u2014 published EOCRC signature"),
               mpatches.Patch(color=TIER_B_COLOR, alpha=0.5, label="Tier B \u2014 established CRC drivers")]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=9)
    plt.tight_layout()
    plt.savefig("results/figures/02_per_gene_distributions.png", dpi=200)
    plt.close()


def fig3_tier_b_heterogeneity():
    with open("results/tables/per_gene_summary.json") as f:
        summary = json.load(f)

    tier_b_summary = [s for s in summary if s["tier"] == "B"]
    tier_b_summary.sort(key=lambda s: s["median_score"])

    genes = [s["gene"] for s in tier_b_summary]
    medians = [s["median_score"] for s in tier_b_summary]
    q1s = [s["q1"] for s in tier_b_summary]
    q3s = [s["q3"] for s in tier_b_summary]
    err_low = [med - q1 for med, q1 in zip(medians, q1s)]
    err_high = [q3 - med for q3, med in zip(q3s, medians)]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    y_pos = range(len(genes))
    ax.errorbar(medians, y_pos, xerr=[err_low, err_high], fmt="o", color=TIER_B_COLOR,
                markersize=8, capsize=4, elinewidth=1.5, ecolor="#888888")
    ax.axvline(0, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes)
    ax.set_xlabel("Median ESM2 score (bars: Tukey's hinges IQR)")
    ax.set_title("Tier B is heterogeneous, not uniformly negative", fontsize=11)
    ax.set_ylim(-1.1, len(genes) - 0.3)
    ax.text(0.02, 0.04,
            "TP53 / KRAS / FBXW7 / BRAF cluster strongly negative; PIK3CA / SMAD4 sit close to zero.",
            transform=ax.transAxes, fontsize=8.5, color="#555555", va="bottom")
    plt.tight_layout()
    plt.savefig("results/figures/03_tier_b_heterogeneity.png", dpi=200)
    plt.close()


def fig4_hotspot_concordance():
    with open("results/tables/hotspot_concordance.json") as f:
        hotspot_data = json.load(f)

    genes = [h["gene"] for h in hotspot_data]
    hs_medians = [h["hotspot_median_score"] for h in hotspot_data]
    nonhs_medians = [h["nonhotspot_median_score"] for h in hotspot_data]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = np.arange(len(genes))
    width = 0.32
    ax.bar(x - width / 2, hs_medians, width, label="Known hotspot positions", color="#55A868")
    ax.bar(x + width / 2, nonhs_medians, width, label="Other observed positions", color="#8172B2")
    ax.axhline(0, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(genes)
    ax.set_ylabel("Median ESM2 score")
    ax.set_title("Hotspot concordance is mixed across genes", fontsize=11)
    ax.legend(frameon=False, fontsize=9)

    notes = ["hotspots NOT\nmore negative", "directionally\nconsistent", "no separation"]
    for xi, gene, note in zip(x, genes, notes):
        top = max(hs_medians[list(genes).index(gene)], nonhs_medians[list(genes).index(gene)])
        ax.text(xi, top + 0.4, note, ha="center", fontsize=7.5, color="#555555", style="italic", va="bottom")

    ax.set_ylim(min(min(hs_medians), min(nonhs_medians)) - 0.6, 1.8)
    plt.tight_layout()
    plt.savefig("results/figures/04_hotspot_concordance.png", dpi=200)
    plt.close()


if __name__ == "__main__":
    fig1_pipeline_overview()
    print("Saved results/figures/01_pipeline_overview.png")
    fig2_per_gene_distributions()
    print("Saved results/figures/02_per_gene_distributions.png")
    fig3_tier_b_heterogeneity()
    print("Saved results/figures/03_tier_b_heterogeneity.png")
    fig4_hotspot_concordance()
    print("Saved results/figures/04_hotspot_concordance.png")
