"""
Step 2: Retrieve canonical UniProt sequences for the 14-gene panel.

Safeguard (added after the FBXW7 / FBXW7-AS1 near-miss during panel design):
we do NOT trust the top search hit by relevance. We require that the returned
entry's own gene_names field contains the EXACT queried symbol, and we flag
(not silently pick) any query returning more than one distinct matching entry.
"""
import requests
import json
import time
from datetime import datetime, timezone

TIER_A = ["ALDOB", "FBXL16", "IL1RN", "METTL27", "MSLN", "RAC3", "SLC38A11", "WNT11"]
TIER_B = ["TP53", "KRAS", "PIK3CA", "FBXW7", "BRAF", "SMAD4"]
PANEL = TIER_A + TIER_B

ESM2_CONTEXT_LIMIT = 1022  # confirmed empirically: max_position_embeddings=1026, minus BOS/EOS


def fetch_candidates(gene_symbol):
    r = requests.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={
            "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
            "fields": "accession,gene_names,length,sequence,protein_name",
            "format": "json",
            "size": 10,
        },
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("results", [])


def verify_and_select(gene_symbol, candidates):
    """Return the verified canonical entry, or raise/flag if ambiguous."""
    exact_matches = []
    for e in candidates:
        gene_names_block = e.get("genes", [{}])[0]
        primary_name = gene_names_block.get("geneName", {}).get("value", "")
        synonyms = [s.get("value", "") for s in gene_names_block.get("synonyms", [])]
        all_names = [primary_name] + synonyms
        if gene_symbol in all_names and primary_name == gene_symbol:
            exact_matches.append(e)

    if len(exact_matches) == 0:
        return None, "NO_EXACT_SYMBOL_MATCH"
    if len(exact_matches) > 1:
        # multiple reviewed entries with this exact primary gene name - genuinely
        # ambiguous, needs manual review, do not auto-pick
        return None, f"AMBIGUOUS: {len(exact_matches)} exact matches"
    return exact_matches[0], "OK"


def main():
    results = {}
    problems = []

    for gene in PANEL:
        candidates = fetch_candidates(gene)
        entry, status = verify_and_select(gene, candidates)

        if entry is None:
            problems.append((gene, status, [c["primaryAccession"] for c in candidates]))
            print(f"  {gene}: *** {status} *** candidates seen: {[c['primaryAccession'] for c in candidates]}")
            time.sleep(0.15)
            continue

        acc = entry["primaryAccession"]
        seq = entry["sequence"]["value"]
        length = entry["sequence"]["length"]
        within_limit = length <= ESM2_CONTEXT_LIMIT

        results[gene] = {
            "gene": gene,
            "uniprot_accession": acc,
            "sequence": seq,
            "protein_length": length,
            "within_esm2_limit": within_limit,
            "context_window_required": not within_limit,
            "protein_name": entry.get("proteinDescription", {})
                                  .get("recommendedName", {})
                                  .get("fullName", {})
                                  .get("value", ""),
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        flag = "OK" if within_limit else "*** EXCEEDS LIMIT - windowing required ***"
        print(f"  {gene}: {acc} | {length} aa | {flag}")
        time.sleep(0.15)

    if problems:
        print("\nFATAL: unresolved gene-symbol ambiguities, cannot proceed silently:")
        for p in problems:
            print(" ", p)
        raise SystemExit(1)

    assert len(results) == 14, f"Expected 14 resolved sequences, got {len(results)}"

    with open("data/raw/sequences.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nAll 14 sequences verified and saved to data/raw/sequences.json")
    long_proteins = [g for g, r in results.items() if r["context_window_required"]]
    print(f"Proteins requiring sliding-window handling (>{ESM2_CONTEXT_LIMIT} aa): {long_proteins}")


if __name__ == "__main__":
    main()
