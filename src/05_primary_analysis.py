"""
Step 6: Primary biological/statistical analysis.

Per the locked methodology: per-gene distributions are the primary unit of
analysis (not pooled scores). Median + IQR, not mean. No percentile-threshold
"strongly disfavored" category (removed per the last design decision) - any
variant callouts use a fixed rank ("top 5 lowest-scoring") with ranking
language only, never "predicted pathogenic."
"""
import json
import statistics as stats
from collections import defaultdict
from pipeline_utils import tukey_hinges_iqr as iqr

TIER_A = ["ALDOB", "FBXL16", "IL1RN", "METTL27", "MSLN", "RAC3", "SLC38A11", "WNT11"]
TIER_B = ["TP53", "KRAS", "PIK3CA", "FBXW7", "BRAF", "SMAD4"]


def main():
    with open("data/scored/scored_variants.json") as f:
        scored = json.load(f)

    by_gene = defaultdict(list)
    for r in scored:
        by_gene[r["gene"]].append(r)

    print("=" * 100)
    print("PER-GENE ESM2 MASKED-MARGINAL SCORE DISTRIBUTIONS (primary analysis)")
    print("Score = logP(mutant) - logP(WT). More negative = more model-disfavored.")
    print("=" * 100)

    gene_summary = []
    for tier, genes in [("A", TIER_A), ("B", TIER_B)]:
        print(f"\n--- TIER {tier} ---")
        for gene in genes:
            rows = by_gene.get(gene, [])
            if not rows:
                print(f"  {gene}: NO SCORED VARIANTS")
                continue
            scores = [r["esm2_score"] for r in rows]
            n = len(scores)
            med = stats.median(scores)
            q1q3 = iqr(scores)
            total_recurrence = sum(r["recurrence_count"] for r in rows)
            windowed = sum(1 for r in rows if r.get("context_window_applied"))

            print(f"  {gene:<10} n_unique={n:<4} n_patient_instances={total_recurrence:<5} "
                  f"median={med:>7.2f}  IQR=[{q1q3[0]:.2f}, {q1q3[1]:.2f}]"
                  f"{'  (windowed)' if windowed else ''}")

            gene_summary.append({
                "gene": gene, "tier": tier, "n_unique_variants": n,
                "n_patient_instances": total_recurrence,
                "median_score": med, "q1": q1q3[0], "q3": q1q3[1],
                "min_score": min(scores), "max_score": max(scores),
                "context_window_applied": windowed > 0,
            })

    # ---- fixed-rank callouts: top 5 most model-disfavored PER GENE, ranking language only ----
    print("\n" + "=" * 100)
    print("TOP 5 MOST MODEL-DISFAVORED VARIANTS PER GENE (rank-based, NOT a pathogenicity claim)")
    print("=" * 100)
    top5_rows = []
    for tier, genes in [("A", TIER_A), ("B", TIER_B)]:
        for gene in genes:
            rows = sorted(by_gene.get(gene, []), key=lambda r: r["esm2_score"])[:5]
            if not rows:
                continue
            print(f"\n  {gene} (Tier {tier}):")
            for rank, r in enumerate(rows, 1):
                print(f"    #{rank}  {r['variant_notation']:<10} score={r['esm2_score']:>7.2f}  "
                      f"seen in {r['recurrence_count']} patient sample(s)")
                top5_rows.append({**r, "rank_within_gene": rank})

    with open("results/tables/per_gene_summary.json", "w") as f:
        json.dump(gene_summary, f, indent=2)
    with open("results/tables/top5_per_gene.json", "w") as f:
        json.dump(top5_rows, f, indent=2)

    # ---- simple CSV exports for easy viewing ----
    import csv
    with open("results/tables/per_gene_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(gene_summary[0].keys()))
        w.writeheader()
        w.writerows(gene_summary)

    with open("results/tables/scored_variants_full.csv", "w", newline="") as f:
        fieldnames = ["gene", "tier", "variant_notation", "uniprot_accession",
                      "protein_length", "context_window_applied", "esm2_score",
                      "recurrence_count", "sample_ids"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in scored:
            row = {k: r.get(k) for k in fieldnames}
            row["sample_ids"] = ";".join(row["sample_ids"])
            w.writerow(row)

    print(f"\nSaved: results/tables/per_gene_summary.json + .csv")
    print(f"Saved: results/tables/top5_per_gene.json")
    print(f"Saved: results/tables/scored_variants_full.csv")


if __name__ == "__main__":
    main()
