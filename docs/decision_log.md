# Decision log

Kept because the reasoning behind a choice is usually more interesting, and more defensible, than the choice itself.

## Gene panel

The panel went through three real iterations before landing at 14 genes.

**First draft (rejected): 28 genes.** Combined the 8 published EOCRC genes with two locally-derived sets from an unrelated RNA-seq project on the same cancer type — the 8 genes with the largest fold-change, and 12 more pulled from that project's pathway diagrams. Rejected for two reasons: the fold-change set turned out to be eight variations on one phenomenon (loss of normal colon-epithelial differentiation genes), not eight independent findings, and two of the fold-change hits weren't even protein-coding. The pathway-diagram set was real biology but was selected because it appeared in a different, smaller project's results, not because of any evidence specific to this study.

**Second draft (rejected): 12 genes**, replacing the fold-change genes with proliferation/genome-instability genes (CDK1, BUB1B, ORC1) chosen for having independently documented disease-causing mutations *in an unrelated clinical context* (mosaic variegated aneuploidy, Meier-Gorlin syndrome). This conflated two different kinds of evidence: "mutations in this gene matter somewhere" isn't the same claim as "this gene is a recognized colorectal cancer driver." None of the three appear among the TCGA colorectal landmark study's significantly mutated genes.

**Final panel: 14 genes, two tiers, each justified by one clean criterion:**

- **Tier A** — membership in a published, externally peer-reviewed EOCRC signature (Marx et al. 2024), regardless of this project's own findings.
- **Tier B** — independent status as a recurrent CRC somatic driver in the TCGA colorectal landmark study (Muzny et al. 2012), with documented missense-specific functional evidence (e.g. FBXW7's R465/R479/R505 hotspots have direct experimental proof from knock-in mouse and cell-line studies that these exact substitutions cause loss of substrate degradation).

An optional third tier (BUB1B, ORC1, selected on dysregulation in the earlier RNA-seq project plus unrelated-disease mutation evidence) was considered and dropped. The rationale for it decomposed into three components on inspection — expression-based selection, unrelated-disease germline evidence, and a pathway framing derived from the same expression data — none of which was "recognized CRC driver" or "CRC-relevant missense functional evidence." Stacking three adjacent-but-weaker evidence types doesn't add up to one strong one.

**Excluded despite being obvious candidates**: APC and TGFBR2, both textbook CRC genes, both excluded from this specific study on a concrete, checked basis rather than overlooked. APC's real mutation spectrum in this cohort is 89% nonsense/frameshift (50 missense out of 456 total mutations) — its dominant biological mechanism in CRC is truncation, which a missense-only scoring method can't meaningfully assess. TGFBR2's classic CRC-inactivating event is a specific 1-base-pair frameshift in a poly-adenine tract (the BAT-RII mutation, occurring in >90% of microsatellite-unstable colon cancers) — again, not a substitution this method scores, and its real mutation count in this cohort was thin regardless (18 total).

## Cohort choice

cBioPortal's public instance groups colon and rectal adenocarcinoma into one combined study (`coadread_tcga_pan_can_atlas_2018`). Rather than choosing between the combined cohort and a separately-sourced colon-only project, the study's own built-in colon-only sample list (`..._coad`, 439 patients) was used — a genuine colon-only subset within the standard, citable public study, avoiding both the imprecision of calling the combined cohort "TCGA-COAD" and the extra complexity of pulling from a different source entirely.

## PIK3CA's context window

ESM2's positional embeddings cap effective context at 1,022 residues, confirmed directly from the model config (not assumed) across every checkpoint size from 35M to 3B parameters. PIK3CA is 1,068 residues. Truncating the sequence was rejected outright — it risks silently cutting off a mutation's local context or, for C-terminal positions, removing them from the sequence entirely. A local window (up to 1,023 residues, centered on the mutated position, shifted inward near either terminus) preserves the immediate sequence context that actually drives the model's local prediction. Every PIK3CA score carries an explicit flag noting the window was applied.

## Dropping the percentile threshold

An early plan defined "strongly model-disfavored" as, e.g., the bottom 10% of scores across the dataset. Dropped before implementation: a percentile is a statement about this dataset's shape, not a biological line, and it would silently shift if the panel ever grew or shrank. Fixed-rank callouts ("the 5 lowest-scoring variants in this gene") make the same point without implying a threshold that doesn't exist.
