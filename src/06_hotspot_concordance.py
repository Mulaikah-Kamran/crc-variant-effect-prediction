"""
Step 7: Hotspot concordance check (SECONDARY analysis, sanity-check role only).

Hotspot positions below are taken from established oncology literature,
decided BEFORE looking at this gene's score distribution:
  - KRAS: codons 12, 13, 61 (Muzny et al. 2012, TCGA colorectal landmark paper)
  - BRAF: codon 600 (V600E, the classic activating mutation)
  - PIK3CA: codons 545, 1047 (E545K/H1047R, the classic activating hotspots)

This is explicitly a concordance/positive-control check, NOT a discovery.
ESM2 flagging these as disfavored confirms the pipeline is working correctly;
it is not new biological information, since these mutations are already
extremely well characterized experimentally.
"""
import json
import statistics as stats

HOTSPOT_POSITIONS = {
    "KRAS": {12, 13, 61},
    "BRAF": {600},
    "PIK3CA": {545, 1047},
}


def main():
    with open("data/scored/scored_variants.json") as f:
        scored = json.load(f)

    print("=" * 90)
    print("HOTSPOT CONCORDANCE CHECK (secondary analysis — sanity check, not a discovery)")
    print("Hotspot positions defined a priori from oncology literature, independent of ESM2 results.")
    print("=" * 90)

    summary = []
    for gene, hotspot_pos in HOTSPOT_POSITIONS.items():
        rows = [r for r in scored if r["gene"] == gene]
        hotspot_rows = [r for r in rows if r["protein_position"] in hotspot_pos]
        other_rows = [r for r in rows if r["protein_position"] not in hotspot_pos]

        print(f"\n{gene} (hotspot codons: {sorted(hotspot_pos)}):")
        if hotspot_rows:
            print(f"  Hotspot variants observed in this cohort ({len(hotspot_rows)} unique):")
            for r in sorted(hotspot_rows, key=lambda x: x["esm2_score"]):
                print(f"    {r['variant_notation']:<10} score={r['esm2_score']:>7.2f}  "
                      f"(n={r['recurrence_count']} patients)")
            hs_scores = [r["esm2_score"] for r in hotspot_rows]
            print(f"  Hotspot median: {stats.median(hs_scores):.2f}")
        else:
            print("  No hotspot-position variants observed in this cohort.")

        if other_rows:
            other_scores = [r["esm2_score"] for r in other_rows]
            print(f"  Non-hotspot median (n={len(other_rows)} unique variants): "
                  f"{stats.median(other_scores):.2f}")

        summary.append({
            "gene": gene,
            "hotspot_positions": sorted(hotspot_pos),
            "n_hotspot_variants_observed": len(hotspot_rows),
            "hotspot_median_score": stats.median([r["esm2_score"] for r in hotspot_rows]) if hotspot_rows else None,
            "n_nonhotspot_variants": len(other_rows),
            "nonhotspot_median_score": stats.median([r["esm2_score"] for r in other_rows]) if other_rows else None,
        })

    print("\n" + "=" * 90)
    print("INTERPRETATION NOTE (record this verbatim in the writeup):")
    print("ESM2 assigning strongly negative scores to known hotspots (e.g. BRAF V600E,")
    print("KRAS G12/G13/Q61 substitutions) is expected and serves only as a pipeline")
    print("sanity check — it does not validate ESM2 for the genuinely open questions")
    print("elsewhere in this panel (Tier A, or non-hotspot Tier B positions), since")
    print("these specific mutations are already among the most functionally")
    print("characterized substitutions in all of cancer biology.")
    print("=" * 90)

    with open("results/tables/hotspot_concordance.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved: results/tables/hotspot_concordance.json")


if __name__ == "__main__":
    main()
