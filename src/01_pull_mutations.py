"""
Step 1: Pull raw mutation data for the locked 14-gene panel.

Study: coadread_tcga_pan_can_atlas_2018 (TCGA COADREAD Pan-Cancer Atlas)
Sample list: coadread_tcga_pan_can_atlas_2018_coad (colon-only subset, 439 patients)
Molecular profile: coadread_tcga_pan_can_atlas_2018_mutations

Pulls ALL mutation types (not just missense) for full transparency / audit trail,
per the locked methodology (Section 2). Filtering to missense-only happens in
script 03, not here.
"""
import requests
import json
import time
from datetime import datetime, timezone

BASE = "https://www.cbioportal.org/api"
STUDY_ID = "coadread_tcga_pan_can_atlas_2018"
PROFILE_ID = f"{STUDY_ID}_mutations"
SAMPLE_LIST_ID = f"{STUDY_ID}_coad"

TIER_A = ["ALDOB", "FBXL16", "IL1RN", "METTL27", "MSLN", "RAC3", "SLC38A11", "WNT11"]
TIER_B = ["TP53", "KRAS", "PIK3CA", "FBXW7", "BRAF", "SMAD4"]
PANEL = TIER_A + TIER_B
TIER_MAP = {g: "A" for g in TIER_A} | {g: "B" for g in TIER_B}

assert len(PANEL) == 14, f"Expected 14 genes, got {len(PANEL)}"
assert len(set(PANEL)) == 14, "Duplicate gene in panel!"


def resolve_gene_ids(genes):
    """Resolve gene symbols to Entrez IDs, with strict verification that the
    returned gene's symbol matches what we asked for (safeguard against
    symbol collisions, e.g. FBXW7 vs FBXW7-AS1)."""
    resolved = {}
    failures = []
    for g in genes:
        r = requests.get(f"{BASE}/genes/{g}", headers={"Accept": "application/json"}, timeout=30)
        if r.status_code != 200:
            failures.append((g, f"HTTP {r.status_code}"))
            continue
        data = r.json()
        if data.get("hugoGeneSymbol") != g:
            failures.append((g, f"symbol mismatch: got {data.get('hugoGeneSymbol')}"))
            continue
        resolved[g] = data["entrezGeneId"]
        time.sleep(0.1)
    return resolved, failures


def fetch_all_mutations(entrez_ids):
    payload = {"entrezGeneIds": entrez_ids, "sampleListId": SAMPLE_LIST_ID}
    r = requests.post(
        f"{BASE}/molecular-profiles/{PROFILE_ID}/mutations/fetch",
        params={"projection": "DETAILED"},
        json=payload,
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        timeout=120,
    )
    r.raise_for_status()
    return r.json()


def main():
    print(f"Resolving {len(PANEL)} gene symbols to Entrez IDs...")
    gene_ids, failures = resolve_gene_ids(PANEL)
    if failures:
        print("FATAL: gene resolution failures:", failures)
        raise SystemExit(1)
    print(f"  Resolved all {len(gene_ids)} genes cleanly (symbol-verified).")
    for g, eid in gene_ids.items():
        print(f"    {g}: entrez {eid} (Tier {TIER_MAP[g]})")

    print(f"\nFetching ALL mutation types from {PROFILE_ID}, sample list {SAMPLE_LIST_ID}...")
    muts = fetch_all_mutations(list(gene_ids.values()))
    print(f"  Retrieved {len(muts)} raw mutation records.")

    # annotate each record with our panel/tier metadata for downstream traceability
    entrez_to_gene = {v: k for k, v in gene_ids.items()}
    for m in muts:
        entrez = m.get("entrezGeneId") or m.get("gene", {}).get("entrezGeneId")
        gene_symbol = entrez_to_gene.get(entrez, m.get("gene", {}).get("hugoGeneSymbol"))
        m["_panel_gene"] = gene_symbol
        m["_panel_tier"] = TIER_MAP.get(gene_symbol, "UNKNOWN")

    out = {
        "metadata": {
            "study_id": STUDY_ID,
            "sample_list_id": SAMPLE_LIST_ID,
            "molecular_profile_id": PROFILE_ID,
            "panel_genes": PANEL,
            "tier_map": TIER_MAP,
            "n_records": len(muts),
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "source": "cBioPortal public API, https://www.cbioportal.org/api",
        },
        "mutations": muts,
    }

    with open("data/raw/mutations_raw.json", "w") as f:
        json.dump(out, f)
    print(f"\nSaved to data/raw/mutations_raw.json")

    # quick per-gene, per-type tally for the transparency table
    from collections import Counter
    tally = Counter()
    for m in muts:
        tally[(m["_panel_gene"], m.get("mutationType", "UNKNOWN"))] += 1

    print(f"\n{'Gene':<10}{'Tier':<6}{'Type':<20}{'Count'}")
    for g in PANEL:
        rows = sorted([(t, c) for (gg, t), c in tally.items() if gg == g], key=lambda x: -x[1])
        if not rows:
            print(f"{g:<10}{TIER_MAP[g]:<6}{'(no mutations found)':<20}")
            continue
        for i, (t, c) in enumerate(rows):
            label = g if i == 0 else ""
            tier = TIER_MAP[g] if i == 0 else ""
            print(f"{label:<10}{tier:<6}{t:<20}{c}")


if __name__ == "__main__":
    main()
