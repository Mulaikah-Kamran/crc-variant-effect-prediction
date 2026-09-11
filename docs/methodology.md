# Methodology

## 1. Mutation data

**Source**: cBioPortal public REST API (`https://www.cbioportal.org/api`), mirroring TCGA's MC3 mutation calls.

**Study**: `coadread_tcga_pan_can_atlas_2018` — TCGA colorectal adenocarcinoma, Pan-Cancer Atlas. This groups colon and rectal adenocarcinoma together, following TCGA's own project structure; there is no standalone "TCGA-COAD" study on the public cBioPortal instance.

**Cohort used**: the study's built-in colon-only sample list, `coadread_tcga_pan_can_atlas_2018_coad` (439 patients), not the full 594-patient colon+rectal set. This keeps the cohort colon-specific without misrepresenting its provenance as a separate TCGA project.

**Retrieved**: all mutation types for the 14 panel genes in one pull (`data/raw/mutations_raw.json`), not pre-filtered to missense — this keeps a complete, auditable record of what was excluded and why. 799 records total, spanning 8 mutation types.

## 2. Mutation filtering

Primary experiment is missense-only. Every record is classified into one of:

`pass` · `excluded_non_missense` · `excluded_unparseable` · `excluded_invalid_residue` · `excluded_out_of_bounds` · `excluded_wt_mismatch`

Nothing is dropped without a logged reason — see `data/processed/audit_all_records.json` for every one of the 799 records.

Indels, nonsense, frameshift, and splice-site mutations are structurally incompatible with masked-marginal scoring (it assumes a fixed-position, single-residue substitution) and are excluded from scoring, not from the audit trail.

Duplicate mutations across patients (the same substitution recurring in multiple tumors) are not discarded — they're deduplicated to one scoring job per unique (gene, position, WT, mutant) tuple for compute efficiency, with a `recurrence_count` and full `sample_ids` list retained per variant.

## 3. Protein sequences

Canonical, reviewed (Swiss-Prot) human UniProt entries. Gene-symbol resolution includes an explicit safeguard: the returned entry's own primary gene name must exactly match the query. This isn't a hypothetical concern — an earlier, unguarded query for `FBXW7` returned an 85-residue fragment belonging to *FBXW7-AS1* (an antisense RNA locus) as its top hit by search relevance. The safeguard rejects any non-exact or ambiguous match rather than trusting rank order.

| Gene | Accession | Length (aa) |
|---|---|---|
| ALDOB | P05062 | 364 |
| FBXL16 | Q8N461 | 479 |
| IL1RN | P18510 | 177 |
| METTL27 | Q8N6F8 | 245 |
| MSLN | Q13421 | 622 |
| RAC3 | P60763 | 192 |
| SLC38A11 | Q08AI6 | 406 |
| WNT11 | O96014 | 354 |
| TP53 | P04637 | 393 |
| KRAS | P01116 | 189 |
| PIK3CA | P42336 | 1068 |
| FBXW7 | Q969H0 | 707 |
| BRAF | P15056 | 766 |
| SMAD4 | Q13485 | 552 |

## 4. Mutation-to-protein mapping and QC

For every candidate variant, in order: position bounds check → WT-residue check (the amino acid stated in the mutation record must exactly match the UniProt sequence at that position) → valid-amino-acid check. A single MSLN record (R413Q) failed the WT-residue check — UniProt shows threonine, not arginine, at position 413, most likely a transcript-numbering discrepancy specific to that record — and was excluded rather than scored. The rest of MSLN's variants passed cleanly, so this was treated as an isolated case, not a systematic transcript mismatch (which would instead call the whole gene's mapping into question).

## 5. ESM2

`facebook/esm2_t33_650M_UR50D`, via HuggingFace `transformers` (CPU inference; ~3–8 seconds per variant depending on load). Confirmed directly from the model config, not assumed: `max_position_embeddings = 1026`, giving an effective limit of 1,022 usable residue positions, identical across every ESM2 checkpoint size from 35M to 3B parameters — it's a property of the architecture's positional embeddings, not something a bigger checkpoint fixes.

**PIK3CA (1,068 residues) exceeds this limit.** Rather than truncating the sequence, each PIK3CA mutation is scored using a local window of up to 1,023 residues centered on the mutated position, shifted inward near either terminus so the window never runs past the sequence's start or end. The mutation's position is re-indexed relative to the window before masking. Every PIK3CA output row carries `context_window_applied: true`. This follows the same general approach used by Brandes et al. (2023) to extend ESM1b-style scoring to proteins beyond the model's native context length.

## 6. Masked-marginal scoring

For a substitution at position *i*, wild-type residue *wt*, mutant residue *mt*:

1. Tokenize the sequence (or window). Residue *i* occupies token index *i* (index 0 is the model's start-of-sequence token).
2. Replace the token at index *i* with the mask token.
3. Forward pass; take the log-softmax of the logits at the masked position.
4. **Score = log P(mt) − log P(wt)**.

More negative → the model finds the mutant substantially less probable than the wild type at that position ("model-disfavored"). Near zero → the model finds both about equally plausible. Positive → the model finds the mutant *more* probable than the wild type there.

No normalization is applied. Raw scores are retained for every variant. **Scores are not on a universal scale across proteins** — they reflect each sequence's own local probability landscape, not a calibrated severity score, and this is more true still for PIK3CA given its windowed context. Cross-gene comparisons in this project are made at the level of per-gene medians/distributions and within-gene rank, not pooled raw magnitudes.

This is a zero-shot score from a self-supervised sequence model — not a trained classifier, and not calibrated against any pathogenicity label.

## 7. Primary analysis

Per-gene distributions are the primary unit of analysis — not pooled scores across all 259 variants, which would let TP53's 78 unique variants and KRAS's 169 patient instances dominate a picture that Tier A's 3–9-variant genes would otherwise vanish into.

**Median, not mean**, per gene, since these distributions aren't expected to be symmetric. **IQR via Tukey's hinges** (median of each half of the sorted sample) — stated explicitly because for small samples this differs visibly from other valid quartile conventions (e.g. linear-interpolation methods); the choice matters enough for genes with 3–6 variants that it needs to be named, not left implicit.

No percentile-based "strongly disfavored" category is used anywhere in this project. An earlier design draft considered one and dropped it deliberately: a percentile threshold describes the shape of this dataset, not a biological line, and would silently redraw itself if the gene panel ever changed. Where individual variants are called out, it's by fixed rank ("the 5 lowest-scoring variants in this gene"), with ranking language only — never "predicted pathogenic."

No pooled Tier A vs. Tier B hypothesis test is performed. Sample sizes are too small and too unequal (3–10 variants per Tier A gene vs. 8–169 unique variants per Tier B gene) for a formal test to mean much; the comparison here is descriptive.

## 8. Hotspot concordance (secondary)

Hotspot positions were fixed from the literature before any scoring was run: KRAS codons 12/13/61, BRAF codon 600, PIK3CA codons 545/1047. This is a concordance check, not a discovery — these mutations are already among the most functionally characterized substitutions in cancer biology, so ESM2 agreeing with them mainly confirms the pipeline is doing what it's supposed to. It does not validate the method for the genuinely open parts of this project (Tier A, or Tier B's non-hotspot positions).

## 9. ClinVar (secondary, coverage-gated)

Eligible variants require an exact match — same gene, protein position, wild-type residue, and mutant residue — against a ClinVar entry; no fuzzy matching. A gene is only included in the pathogenic-vs-benign comparison if at least 10 matched variants exist on both sides. Five of six Tier B genes failed this gate, mostly because ClinVar's benign-variant coverage for high-profile cancer genes is thin (clinical sequencing overwhelmingly reports and curates pathogenic findings). ClinVar's classifications reflect germline clinical interpretation, which is a different biological question from somatic tumor function — this is treated as an external cross-reference, not a gold-standard label for what these mutations do in a colorectal tumor.

## A note on the −20 to +5 range mentioned in the QC logs

This range is an implementation sanity check only — a bound used to flag any score so extreme it might indicate an indexing bug, not a biological validity threshold. No score in the final dataset fell outside it.
