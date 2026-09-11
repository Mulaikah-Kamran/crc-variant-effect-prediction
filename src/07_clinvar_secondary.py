"""
Step 8: SECONDARY ClinVar cross-check (coverage-gated, exact-match only).

Critical distinction preserved throughout: ClinVar classifications are
predominantly GERMLINE clinical variant interpretations, not CRC somatic
tumor biology. This checks whether ESM2's sequence-level score externally
correlates with an independent clinical labeling system - it is NOT treated
as a gold-standard label for what these mutations do in a colorectal tumor.

Eligibility: exact match on (gene, protein_position, wt_residue, mt_residue)
to one of our already-scored variants. No fuzzy matching.

Gating rule (locked): a gene is excluded from the comparison if it has fewer
than 10 Pathogenic/Likely_Pathogenic OR fewer than 10 Benign/Likely_Benign
*matching* variants - reported as insufficient coverage, not forced into a
comparison with a handful of points per side.
"""
import json
import re
import time
import requests
import statistics as stats

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MIN_PER_CATEGORY = 10

PROTEIN_CHANGE_RE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

GENES_TO_CHECK = ["TP53", "KRAS", "BRAF", "PIK3CA", "FBXW7", "SMAD4"]  # Tier B only - Tier A too sparse/obscure for meaningful ClinVar coverage


def fetch_clinvar_variants(gene, classification_group):
    """classification_group: 'pathogenic' or 'benign' search term."""
    term = f"{gene}[gene] AND {classification_group}[Clinical_Significance] AND single_gene[Type]"
    r = requests.get(f"{EUTILS}/esearch.fcgi", params={
        "db": "clinvar", "term": term, "retmax": 500, "retmode": "json"
    }, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    results = []
    # esummary in batches of 100
    for i in range(0, len(ids), 100):
        batch = ids[i:i+100]
        r2 = requests.get(f"{EUTILS}/esummary.fcgi", params={
            "db": "clinvar", "id": ",".join(batch), "retmode": "json"
        }, timeout=30)
        r2.raise_for_status()
        data = r2.json().get("result", {})
        for uid in batch:
            entry = data.get(uid)
            if not entry:
                continue
            title = entry.get("title", "")
            m = PROTEIN_CHANGE_RE.search(title)
            if not m:
                continue
            wt3, pos, mt3 = m.group(1), int(m.group(2)), m.group(3)
            wt1 = AA3_TO_1.get(wt3)
            mt1 = AA3_TO_1.get(mt3)
            if not wt1 or not mt1:
                continue
            results.append({
                "gene": gene, "protein_position": pos,
                "wt_residue": wt1, "mutant_residue": mt1,
                "classification": classification_group,
                "clinvar_title": title,
            })
        time.sleep(0.35)  # be polite to NCBI eutils rate limits
    return results


def main():
    with open("data/scored/scored_variants.json") as f:
        scored = json.load(f)
    scored_lookup = {
        (r["gene"], r["protein_position"], r["wt_residue"], r["mutant_residue"]): r
        for r in scored
    }

    print("=" * 90)
    print("SECONDARY CLINVAR CROSS-CHECK (germline clinical labels vs. ESM2 sequence score)")
    print("NOT a gold-standard label for CRC somatic biology - see methodology Section 12.")
    print("=" * 90)

    all_matches = []
    gene_reports = []

    for gene in GENES_TO_CHECK:
        print(f"\n{gene}:")
        path_variants = fetch_clinvar_variants(gene, "pathogenic")
        benign_variants = fetch_clinvar_variants(gene, "benign")

        matched_path = [v for v in path_variants
                         if (v["gene"], v["protein_position"], v["wt_residue"], v["mutant_residue"]) in scored_lookup]
        matched_benign = [v for v in benign_variants
                            if (v["gene"], v["protein_position"], v["wt_residue"], v["mutant_residue"]) in scored_lookup]

        # dedupe (same variant can appear multiple times in ClinVar under different submissions)
        def dedupe(vs):
            seen = {}
            for v in vs:
                key = (v["protein_position"], v["wt_residue"], v["mutant_residue"])
                seen[key] = v
            return list(seen.values())

        matched_path = dedupe(matched_path)
        matched_benign = dedupe(matched_benign)

        n_path, n_benign = len(matched_path), len(matched_benign)
        print(f"  ClinVar P/LP variants matching our scored set: {n_path}")
        print(f"  ClinVar B/LB variants matching our scored set: {n_benign}")

        if n_path < MIN_PER_CATEGORY or n_benign < MIN_PER_CATEGORY:
            print(f"  -> INSUFFICIENT COVERAGE (need >={MIN_PER_CATEGORY} per category). Excluded from comparison.")
            gene_reports.append({
                "gene": gene, "n_pathogenic_matched": n_path, "n_benign_matched": n_benign,
                "sufficient_coverage": False,
            })
            continue

        path_scores = [scored_lookup[(gene, v["protein_position"], v["wt_residue"], v["mutant_residue"])]["esm2_score"]
                       for v in matched_path]
        benign_scores = [scored_lookup[(gene, v["protein_position"], v["wt_residue"], v["mutant_residue"])]["esm2_score"]
                         for v in matched_benign]
        print(f"  Pathogenic/LP median ESM2 score: {stats.median(path_scores):.2f}")
        print(f"  Benign/LB median ESM2 score:     {stats.median(benign_scores):.2f}")

        gene_reports.append({
            "gene": gene, "n_pathogenic_matched": n_path, "n_benign_matched": n_benign,
            "sufficient_coverage": True,
            "pathogenic_median_score": stats.median(path_scores),
            "benign_median_score": stats.median(benign_scores),
        })
        all_matches.extend(matched_path)
        all_matches.extend(matched_benign)

    with open("results/tables/clinvar_secondary_check.json", "w") as f:
        json.dump(gene_reports, f, indent=2)
    print(f"\nSaved: results/tables/clinvar_secondary_check.json")
    print("\nReminder for writeup: any gene shown here reflects GERMLINE clinical")
    print("interpretation, cross-referenced against our SOMATIC CRC-observed variant")
    print("set purely by coincidence of exact protein position/change overlap.")


if __name__ == "__main__":
    main()
