"""
Classical baseline comparison: BLOSUM62 vs ESM2 zero-shot score.

Scientific question: does the ESM2 zero-shot mutation score behave similarly
to, or differently from, a simple classical amino-acid substitution-matrix
score, for the same 259 already-scored, frozen missense variants?

This is a characterization/baseline comparison, not a claim that either
method is biologically or clinically superior. It does NOT touch the frozen
ESM2 scores, the gene panel, or the mutation filtering - it reads the
existing scored dataset and adds one new derived column plus summary
statistics.

BLOSUM62 interpretation:
  - higher values  -> more evolutionarily tolerated/common substitution
  - lower/negative -> less favored substitution
  - NOT a probability, NOT a pathogenicity score, NOT a functional assay
  - a classical pairwise substitution baseline, not a cancer-specific predictor
"""
import json
import csv
from collections import defaultdict
from scipy.stats import spearmanr

from pipeline_utils import blosum62_score

MIN_N_FOR_PER_GENE_CORRELATION = 5  # below this, a correlation coefficient
                                     # is not reported as meaningful - see
                                     # docs/methodology.md


def main():
    with open("data/scored/scored_variants.json") as f:
        scored = json.load(f)

    assert len(scored) == 259, f"Expected 259 frozen variants, found {len(scored)} - stopping."

    # ---- compute BLOSUM62 score for every variant, alongside (not replacing) ESM2 ----
    enriched = []
    for r in scored:
        b62 = blosum62_score(r["wt_residue"], r["mutant_residue"])
        enriched.append({**r, "blosum62_score": b62})

    assert len(enriched) == 259
    assert all("blosum62_score" in r and r["blosum62_score"] is not None for r in enriched)
    assert all(r["esm2_score"] == orig["esm2_score"] for r, orig in zip(enriched, scored)), \
        "ESM2 scores must be unchanged - this script only adds a derived column."

    # duplicate-key check, same convention as the original QC (script 03)
    keys = [(r["gene"], r["protein_position"], r["wt_residue"], r["mutant_residue"]) for r in enriched]
    assert len(set(keys)) == len(keys), "Duplicate unique-variant key detected - stopping."

    # ---- overall Spearman correlation, n=259 ----
    esm2_scores = [r["esm2_score"] for r in enriched]
    b62_scores = [r["blosum62_score"] for r in enriched]
    overall_rho, overall_p = spearmanr(esm2_scores, b62_scores)

    print("=" * 80)
    print("BLOSUM62 vs ESM2 — classical substitution-matrix baseline comparison")
    print("=" * 80)
    print(f"\nOverall (n={len(enriched)}): Spearman rho = {overall_rho:.3f}, p = {overall_p:.3e}")

    overall_result = {"n": len(enriched), "spearman_rho": overall_rho, "p_value": overall_p}

    # ---- per-gene Spearman correlation, with small-n handling ----
    by_gene = defaultdict(list)
    for r in enriched:
        by_gene[r["gene"]].append(r)

    per_gene_results = []
    print(f"\n{'Gene':<10}{'Tier':<6}{'n':<6}{'Spearman rho':<16}{'p-value':<12}{'Note'}")
    for gene, rows in sorted(by_gene.items()):
        n = len(rows)
        tier = rows[0]["tier"]
        gene_esm2 = [r["esm2_score"] for r in rows]
        gene_b62 = [r["blosum62_score"] for r in rows]

        if n < MIN_N_FOR_PER_GENE_CORRELATION:
            note = f"n<{MIN_N_FOR_PER_GENE_CORRELATION}: correlation not reported as meaningful"
            print(f"{gene:<10}{tier:<6}{n:<6}{'--':<16}{'--':<12}{note}")
            per_gene_results.append({
                "gene": gene, "tier": tier, "n": n,
                "spearman_rho": None, "p_value": None,
                "note": note,
            })
            continue

        # a correlation requires variance in both variables; guard against
        # a degenerate all-identical-score edge case rather than letting
        # scipy emit a silent NaN
        if len(set(gene_esm2)) < 2 or len(set(gene_b62)) < 2:
            note = "insufficient variation in scores for a meaningful correlation"
            print(f"{gene:<10}{tier:<6}{n:<6}{'--':<16}{'--':<12}{note}")
            per_gene_results.append({
                "gene": gene, "tier": tier, "n": n,
                "spearman_rho": None, "p_value": None,
                "note": note,
            })
            continue

        rho, p = spearmanr(gene_esm2, gene_b62)
        note = ""
        print(f"{gene:<10}{tier:<6}{n:<6}{rho:<16.3f}{p:<12.3e}{note}")
        per_gene_results.append({
            "gene": gene, "tier": tier, "n": n,
            "spearman_rho": rho, "p_value": p, "note": note,
        })

    # ---- save outputs ----
    with open("results/tables/blosum62_comparison_overall.json", "w") as f:
        json.dump(overall_result, f, indent=2)
    with open("results/tables/blosum62_comparison_per_gene.json", "w") as f:
        json.dump(per_gene_results, f, indent=2)

    fieldnames = ["gene", "tier", "variant_notation", "protein_position", "wt_residue",
                  "mutant_residue", "esm2_score", "blosum62_score", "recurrence_count",
                  "context_window_applied"]
    with open("results/tables/esm2_blosum62_comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in enriched:
            w.writerow({k: r.get(k) for k in fieldnames})

    with open("results/tables/blosum62_comparison_per_gene.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "tier", "n", "spearman_rho", "p_value", "note"])
        w.writeheader()
        w.writerows(per_gene_results)

    print(f"\nSaved: results/tables/esm2_blosum62_comparison.csv ({len(enriched)} rows)")
    print(f"Saved: results/tables/blosum62_comparison_overall.json")
    print(f"Saved: results/tables/blosum62_comparison_per_gene.json + .csv")

    print("\nInterpretation note (record verbatim in any writeup):")
    print("A strong positive correlation would mean ESM2's ranking of these mutations")
    print("broadly resembles the classical substitution-matrix baseline. A weak")
    print("correlation would mean ESM2 is capturing sequence-context information that")
    print("differs from simple pairwise substitution preferences. Correlation alone")
    print("does not establish that either method is biologically or clinically superior.")

    generate_figure(enriched, overall_rho, len(enriched))


def generate_figure(enriched, overall_rho, n):
    """One scatter figure: BLOSUM62 (x) vs ESM2 (y), colored by tier (not by
    all 14 individual genes, to avoid an oversized legend on the primary
    figure) - consistent with the Tier A / Tier B color convention already
    used in results/figures/02_per_gene_distributions.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })
    TIER_A_COLOR = "#4C72B0"
    TIER_B_COLOR = "#C44E52"

    fig, ax = plt.subplots(figsize=(7, 6))
    for tier, color, label in [("A", TIER_A_COLOR, "Tier A — published EOCRC signature"),
                                ("B", TIER_B_COLOR, "Tier B — established CRC drivers")]:
        xs = [r["blosum62_score"] for r in enriched if r["tier"] == tier]
        ys = [r["esm2_score"] for r in enriched if r["tier"] == tier]
        ax.scatter(xs, ys, s=22, color=color, alpha=0.6, label=label, edgecolors="none")

    ax.set_xlabel("BLOSUM62 substitution score (classical baseline)")
    ax.set_ylabel("ESM2 masked-marginal score")
    ax.set_title("ESM2 vs. BLOSUM62: classical substitution-matrix baseline comparison", fontsize=11)
    ax.axhline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=0)
    ax.axvline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=0)
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
    ax.text(0.03, 0.97, f"Spearman ρ = {overall_rho:.2f}\nn = {n}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc"))

    plt.tight_layout()
    plt.savefig("results/figures/05_blosum62_vs_esm2.png", dpi=200)
    plt.close()
    print(f"\nSaved: results/figures/05_blosum62_vs_esm2.png")


if __name__ == "__main__":
    main()
