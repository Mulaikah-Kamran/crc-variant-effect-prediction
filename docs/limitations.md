# Limitations

**Zero-shot, not trained on this task.** ESM2 has never seen a "cancer driver" or "pathogenic" label. Every score here comes from a model trained only to predict masked amino acids from surrounding sequence context. Agreement with known cancer biology is informative; it is not validation of a purpose-built method, because there isn't one here.

**Sequence only.** No 3D structure, no protein-protein interaction context, no tissue- or cell-type-specific expression is used. A position that looks tolerant in isolated sequence context could still be structurally or functionally critical for reasons ESM2 has no way to see.

**No experimental confirmation.** Every "model-disfavored" or "damaging" statement in this project refers to a sequence probability, not a measured effect on protein stability, activity, or patient outcome.

**Very unequal data per gene.** TP53 contributes 78 unique variants; several Tier A genes contribute 3. Per-gene summaries are reported separately for exactly this reason — a pooled statistic would be dominated by whichever gene happened to have the most mutations in this cohort, not whichever gene is most informative.

**Tier A's coverage is genuinely sparse.** 3–9 unique variants per gene is enough to describe what was observed, not enough to draw a confident conclusion about the gene generally. Anything notable in Tier A (e.g. WNT11's distribution sitting well below zero) is a lead for a larger follow-up, not a finding on its own.

**Known hotspots aren't new information.** KRAS G12/G13/Q61, BRAF V600E, and PIK3CA E545/H1047 are among the most experimentally characterized mutations in oncology. The hotspot-concordance analysis in this project is a sanity check on the pipeline, not a test of anything previously unknown — and even as a sanity check, the result was mixed rather than uniformly clean (see main results).

**TP53's mutation spectrum is mechanistically heterogeneous.** Its real mutations include simple loss-of-function, dominant-negative, and occasional gain-of-function substitutions. A single scalar sequence-disfavor score cannot distinguish between these mechanisms, even where it correctly flags a position as generally intolerant of substitution.

**Cohort vs. population.** TCGA's colon-only Pan-Cancer Atlas cohort (439 patients, largely North American, collected under TCGA-era consent and sequencing protocols) is not necessarily representative of colorectal cancer patients generally, and is a different population from any other cohort this project's findings might eventually be compared against.

**PIK3CA's scores aren't obtained under identical conditions to the rest of the panel.** Its local-window scoring (Section 5 of the methodology) means its raw scores carry an extra source of context-dependence not present for the other 13 genes.

**ClinVar's germline labels answer a different question than this project asks.** A variant's classification as germline-pathogenic (typically from inherited-disease testing) doesn't necessarily track its behavior as a somatic mutation in a colorectal tumor — different selective pressures, different tissue context, sometimes an entirely different clinical syndrome (as appears to be the case for most of KRAS's ClinVar entries, which largely reflect a developmental RASopathy, not oncogenic signaling). The ClinVar check in this project is an external cross-reference, treated accordingly — not a gold-standard label.

**Raw score magnitudes aren't directly comparable across genes** in an absolute sense, since each score reflects a given protein's own local probability landscape rather than a calibrated universal scale. Per-gene distributions and within-gene rank are more defensible than comparing, say, "gene X's most negative score" against "gene Y's most negative score" as if they meant the same thing.
