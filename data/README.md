# Data manifest

## `raw/`
Unmodified API responses. Never edited after retrieval.

- `mutations_raw.json` — all mutation types (not pre-filtered) for the 14-gene panel, from cBioPortal study `coadread_tcga_pan_can_atlas_2018`, sample list `coadread_tcga_pan_can_atlas_2018_coad` (439 patients). 799 records. Includes retrieval timestamp and the gene/tier mapping used.
- `sequences.json` — canonical reviewed UniProt sequences for all 14 genes, with accession, length, and retrieval timestamp.

## `processed/`
- `audit_all_records.json` — every one of the 799 raw records, each with a `qc_status` (`pass`, `excluded_non_missense`, `excluded_unparseable`, `excluded_invalid_residue`, `excluded_out_of_bounds`, or `excluded_wt_mismatch`) and an explicit reason. Nothing is silently dropped — this file is the full audit trail.
- `unique_variants_to_score.json` — the 259 unique (gene, position, WT, mutant) variants that passed QC, deduplicated for scoring, each carrying a `recurrence_count` and the full list of patient sample IDs it was observed in.

## `scored/`
- `scored_variants.json` — the 259 unique variants with their ESM2 masked-marginal scores, including `context_window_applied` for the 51 PIK3CA variants scored via local windowing.

## Regenerating

Run `src/01` through `src/04` in order; each stage reads only the previous stage's output. `src/04_score_esm2.py` checkpoints its progress to `data/scored/_checkpoint.json` (not tracked in git) and can be safely re-run to resume an interrupted scoring run.
