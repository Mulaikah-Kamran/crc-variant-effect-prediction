# Predicting the functional impact of colorectal cancer mutations with ESM2

A protein language model, applied zero-shot to real somatic mutations from colorectal tumors, asking a simple question: does the model find these substitutions surprising?

## The question

Cancer sequencing turns up thousands of mutations per tumor type, most of them never functionally tested. ESM2 is a protein language model trained on millions of sequences to do one thing: given a protein with one amino acid hidden, guess what belongs there. That's it — no cancer biology, no labeled training data, nothing task-specific.

It turns out that "guess what belongs here" is a decent proxy for "is this position under selective pressure to stay the same." If the model is confident a position should be, say, an arginine, and a real tumor has a mutation swapping in a tryptophan, the model's surprise at that substitution is a signal worth checking against what we already know about the gene.

This project scores real missense mutations from a public colorectal cancer cohort with ESM2 and asks whether that signal lines up with things that are already independently known — and, just as importantly, reports the cases where it doesn't.

## Data

- **Mutations**: TCGA COADREAD (colon + rectal adenocarcinoma, Pan-Cancer Atlas), restricted to the colon-only sample subset (439 patients), pulled via the [cBioPortal API](https://www.cbioportal.org/api).
- **Genes**: 14 total, in two tiers.
  - **Tier A** — the 8-gene early-onset colorectal cancer signature from Marx et al. (2024): `ALDOB`, `FBXL16`, `IL1RN`, `METTL27` (published as *WBSCR27*), `MSLN`, `RAC3`, `SLC38A11`, `WNT11`.
  - **Tier B** — six genes independently established as recurrent CRC somatic drivers by the TCGA colorectal landmark study (Muzny et al., 2012): `TP53`, `KRAS`, `PIK3CA`, `FBXW7`, `BRAF`, `SMAD4`.
- **Sequences**: canonical reviewed human entries from UniProt.

Full provenance, including why these specific genes and not others (APC and TGFBR2 were considered and excluded — see below), is in [`docs/decision_log.md`](docs/decision_log.md).

## What ESM2 actually does here

For each observed mutation (wild-type residue → mutant residue at position *i*):

1. Mask position *i* in the protein sequence.
2. Ask ESM2 for its predicted probability of every amino acid at that position.
3. Score = log P(mutant) − log P(wild-type).

A strongly negative score means the model finds the mutant far less likely than the original residue at that spot — the substitution is *model-disfavored*. That's a statement about sequence conservation, not a diagnosis. It is **not**:

- a pathogenicity call
- a probability of clinical harm
- an experimentally measured effect on protein function
- specific to cancer at all (ESM2 has never seen a label indicating "cancer mutation")

Full method, including how proteins longer than ESM2's ~1,022-residue context window are handled (PIK3CA is 1,068 residues), is in [`docs/methodology.md`](docs/methodology.md).

## Results

259 unique missense variants were scored, from 673 mutation instances passing quality control across 439 patients.

![Per-gene score distributions](results/figures/02_per_gene_distributions.png)

**Tier B is not uniformly negative.** TP53, KRAS, FBXW7, and BRAF cluster strongly model-disfavored. PIK3CA and SMAD4 sit close to zero.

![Tier B heterogeneity](results/figures/03_tier_b_heterogeneity.png)

**Checking against known hotspots gives a mixed answer, not a clean validation.** KRAS's canonical hotspot codons (G12/G13/Q61) score less negative than the gene's other observed positions — the opposite of the naive expectation. BRAF's V600E does score more negative than BRAF's other positions, in the expected direction. PIK3CA shows essentially no separation between hotspot and non-hotspot positions.

![Hotspot concordance](results/figures/04_hotspot_concordance.png)

**A secondary ClinVar cross-check** (matching exact protein positions against germline pathogenic/benign classifications — a different biological context from somatic tumor mutations, discussed in the limitations) was only adequately powered for one gene, KRAS, where pathogenic- and benign-labeled variants had nearly identical median scores. Most of KRAS's ClinVar entries reflect a RASopathy developmental syndrome, not oncogenic function — a reasonable explanation for why a germline label wouldn't track a sequence-conservation score tuned toward cancer relevance.

Full tables: [`results/tables/`](results/tables/).

## What this doesn't show

- ESM2 wasn't trained on cancer data, and correctly flagging a well-known hotspot like BRAF V600E is closer to a sanity check than a discovery — that mutation is already one of the most experimentally characterized substitutions in oncology.
- Tier A's mutation counts are small (3–10 unique variants per gene); anything found there is a lead, not a conclusion.
- Raw scores aren't a single universal scale — PIK3CA's scores come from a local sequence window (see methodology) rather than its full protein, so cross-gene score comparisons are made cautiously, and per-gene rank is generally more meaningful than raw magnitude.
- No wet-lab or clinical data confirms any of these predictions.

Full discussion: [`docs/limitations.md`](docs/limitations.md).

## Reproducing this

```bash
pip install -r requirements.txt
python src/01_pull_mutations.py        # cBioPortal -> data/raw/
python src/02_fetch_sequences.py       # UniProt -> data/raw/
python src/03_filter_and_qc.py         # -> data/processed/
python src/04_score_esm2.py            # -> data/scored/  (checkpoints automatically; ~15-20 min on CPU)
python src/05_primary_analysis.py      # -> results/tables/
python src/06_hotspot_concordance.py   # -> results/tables/
python src/07_clinvar_secondary.py     # -> results/tables/  (hits NCBI E-utilities, rate-limited)
python src/08_generate_figures.py      # -> results/figures/
```

Tests: `pytest tests/`

## Repository map

```
src/                  pipeline scripts, run in numeric order
data/raw/             unmodified API pulls (mutations, sequences)
data/processed/       QC'd and deduplicated variant table, plus the full exclusion audit trail
data/scored/          final ESM2-scored variant table
results/tables/       per-gene summaries, hotspot and ClinVar checks
results/figures/      the four figures above
docs/                 methodology, decision log, limitations
tests/                unit tests for the parsing, QC, and windowing logic
```

## Method reference

Masked-marginal scoring: Meier, J. et al. "Language models enable zero-shot prediction of the effects of mutations on protein function." *NeurIPS* 2021.

Model: Lin, Z. et al. "Evolutionary-scale prediction of atomic-level protein structure with a language model." *Science* 2023. (`facebook/esm2_t33_650M_UR50D`, via HuggingFace `transformers`)

EOCRC signature: Marx, O.M. et al. *Frontiers in Oncology* 2024.

CRC driver frequencies: Muzny, D. et al. "Comprehensive molecular characterization of human colon and rectal cancer." *Nature* 2012.
