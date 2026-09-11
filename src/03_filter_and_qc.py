"""
Step 3: Filter, map, and QC every mutation record.

For every record: parse protein change -> validate WT residue against the
canonical UniProt sequence -> check position bounds -> check valid amino acids
-> assign qc_status. Missense-only records that pass QC are deduplicated to
unique (gene, position, wt, mt) tuples for ESM2 scoring; full instance-level
detail (including non-missense types) is retained in the audit table.
"""
import json
from collections import defaultdict
from pipeline_utils import parse_protein_change, validate_variant


def main():
    with open("data/raw/mutations_raw.json") as f:
        raw = json.load(f)
    with open("data/raw/sequences.json") as f:
        sequences = json.load(f)

    mutations = raw["mutations"]
    tier_map = raw["metadata"]["tier_map"]

    audit_rows = []
    exclusion_counts = defaultdict(int)

    for m in mutations:
        gene = m["_panel_gene"]
        tier = m["_panel_tier"]
        mtype = m.get("mutationType", "UNKNOWN")
        pc = m.get("proteinChange", "")
        sample_id = m.get("sampleId") or m.get("patientId") or "UNKNOWN"

        row = {
            "gene": gene,
            "tier": tier,
            "mutation_type": mtype,
            "protein_change_raw": pc,
            "sample_id": sample_id,
            "qc_status": None,
            "exclusion_reason": None,
            "protein_position": None,
            "wt_residue": None,
            "mutant_residue": None,
        }

        # Rule 1: primary experiment is missense-only. Non-missense retained
        # for the audit trail but explicitly marked, not scored.
        if mtype != "Missense_Mutation":
            row["qc_status"] = "excluded_non_missense"
            row["exclusion_reason"] = f"mutation_type={mtype}, not eligible for masked-marginal scoring"
            exclusion_counts[row["qc_status"]] += 1
            audit_rows.append(row)
            continue

        parsed = parse_protein_change(pc)
        if parsed is None:
            row["qc_status"] = "excluded_unparseable"
            row["exclusion_reason"] = f"could not parse protein_change='{pc}' as simple WT{{pos}}MT"
            exclusion_counts[row["qc_status"]] += 1
            audit_rows.append(row)
            continue

        wt, pos, mt = parsed
        row["protein_position"] = pos
        row["wt_residue"] = wt
        row["mutant_residue"] = mt

        seq_info = sequences.get(gene)
        if seq_info is None:
            row["qc_status"] = "excluded_no_sequence"
            row["exclusion_reason"] = "no UniProt sequence resolved for this gene"
            exclusion_counts[row["qc_status"]] += 1
            audit_rows.append(row)
            continue

        seq = seq_info["sequence"]
        qc_status, exclusion_reason = validate_variant(seq, pos, wt, mt)

        if qc_status != "pass":
            row["qc_status"] = qc_status
            row["exclusion_reason"] = exclusion_reason.replace(
                "sequence has", f"UniProt {seq_info['uniprot_accession']} has"
            )
            exclusion_counts[row["qc_status"]] += 1
            audit_rows.append(row)
            continue

        # passed all checks
        row["qc_status"] = "pass"
        row["uniprot_accession"] = seq_info["uniprot_accession"]
        row["protein_length"] = seq_info["protein_length"]
        row["context_window_applied"] = seq_info["context_window_required"]
        audit_rows.append(row)

    # ---- systematic-mismatch check: flag any gene where wt_mismatch failures
    # cluster, since that suggests a transcript/isoform problem, not noise ----
    mismatch_by_gene = defaultdict(int)
    missense_total_by_gene = defaultdict(int)
    for r in audit_rows:
        if r["mutation_type"] == "Missense_Mutation":
            missense_total_by_gene[r["gene"]] += 1
            if r["qc_status"] == "excluded_wt_mismatch":
                mismatch_by_gene[r["gene"]] += 1

    print("=== WT-residue validation summary per gene ===")
    systematic_flags = []
    for gene in sorted(missense_total_by_gene):
        total = missense_total_by_gene[gene]
        mismatches = mismatch_by_gene.get(gene, 0)
        rate = mismatches / total if total else 0
        flag = ""
        if rate > 0.2 and total >= 3:
            flag = "  *** SYSTEMATIC MISMATCH RATE - possible transcript issue, needs manual review ***"
            systematic_flags.append(gene)
        print(f"  {gene:<10} {total-mismatches}/{total} passed WT check{flag}")

    if systematic_flags:
        print(f"\nFATAL: systematic WT-mismatch detected in {systematic_flags} - stopping for manual review.")
        raise SystemExit(1)

    # ---- deduplicate passing missense variants to unique (gene,pos,wt,mt) for scoring ----
    unique_variants = {}
    for r in audit_rows:
        if r["qc_status"] != "pass":
            continue
        key = (r["gene"], r["protein_position"], r["wt_residue"], r["mutant_residue"])
        if key not in unique_variants:
            unique_variants[key] = {
                "gene": r["gene"],
                "tier": r["tier"],
                "uniprot_accession": r["uniprot_accession"],
                "protein_length": r["protein_length"],
                "context_window_applied": r["context_window_applied"],
                "protein_position": r["protein_position"],
                "wt_residue": r["wt_residue"],
                "mutant_residue": r["mutant_residue"],
                "variant_notation": f"{r['wt_residue']}{r['protein_position']}{r['mutant_residue']}",
                "sample_ids": [],
            }
        unique_variants[key]["sample_ids"].append(r["sample_id"])

    for v in unique_variants.values():
        v["recurrence_count"] = len(v["sample_ids"])

    with open("data/processed/audit_all_records.json", "w") as f:
        json.dump(audit_rows, f, indent=2)
    with open("data/processed/unique_variants_to_score.json", "w") as f:
        json.dump(list(unique_variants.values()), f, indent=2)

    print(f"\n=== Overall summary ===")
    print(f"Total raw records:              {len(mutations)}")
    print(f"Excluded (non-missense):        {exclusion_counts['excluded_non_missense']}")
    print(f"Excluded (unparseable):         {exclusion_counts['excluded_unparseable']}")
    print(f"Excluded (invalid residue):     {exclusion_counts['excluded_invalid_residue']}")
    print(f"Excluded (out of bounds):       {exclusion_counts['excluded_out_of_bounds']}")
    print(f"Excluded (WT mismatch):         {exclusion_counts['excluded_wt_mismatch']}")
    passed = sum(1 for r in audit_rows if r["qc_status"] == "pass")
    print(f"Passed QC (missense instances): {passed}")
    print(f"Unique variants to score:       {len(unique_variants)}")
    print(f"\nSaved: data/processed/audit_all_records.json")
    print(f"Saved: data/processed/unique_variants_to_score.json")


if __name__ == "__main__":
    main()
