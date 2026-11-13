# Multiomics feasibility and remaining gates — 2026-11-13

This is a scope/dependency assessment for GOD-15, GOD-16 and GOD-71. It is not an integrated result or a claim that no suitable public data exists.

## Current evidence

The verified local benchmark is GSM8559287, a single mouse E7.0 10x RNA sample. No verified local mass-spectrometry quantitative table, protein biological-replicate map, ATAC matrix/fragment collection or RNA–ATAC barcode-pairing table is available in this project. E7.5 candidate RNA is posterior-region tissue; a whole-embryo protein experiment cannot be treated as matched material without review. Original processing, genome/annotation version and embryo pooling remain unresolved (see `E7.0_source_review.md`).

The [PXD041309 PRIDE record](https://www.ebi.ac.uk/pride/archive/projects/PXD041309) remains the candidate associated with the project's E7.5–E9.5 proteomics paper (DOI 10.1016/j.stem.2024.04.017). This turn confirmed the public project landing page, but it did not obtain and validate a quantitative table or replicate annotations. The landing page alone is not acceptance evidence for protein ingestion. Distinguish real embryos from gastruloids in that study.

## Reuse choice

| Tool | Proposed use | Gate before installation/model fitting |
| --- | --- | --- |
| [scvi-tools](https://github.com/scverse/scvi-tools), BSD-3-Clause | Cross-batch RNA models; MultiVI only with appropriate RNA/ATAC data | Real batch design, barcode pairing evidence, feature/build compatibility, baseline comparison and resource/version check |
| [AlphaPeptStats](https://github.com/MannLabs/alphapeptstats), Apache-2.0 | Native mass-spectrometry downstream quantification workflow | Real quantitative table, units, protein groups, missingness, tissue/stage and biological replicates |

Do not fit totalVI to mass-spectrometry abundance: the protein observation type is different from CITE-seq antibody counts. Do not invent cell-level pairs between unrelated RNA and protein samples. Initial RNA–protein comparisons should be explicitly limited to compatible stage/tissue/replicate units if actual inputs support that level.

## Concrete next work

1. Obtain the candidate's actual protein quantification files and sample/replicate manifest; record file IDs, checksums, search/quantification method and allowed terms. This is data-acquisition work, not necessarily a user-only dependency.
2. Populate the modality contracts and sample mapping with source-supported values; retain unknown IDs.
3. Choose a validated comparison unit and baseline. If a paired RNA/ATAC case is pursued, independently acquire paired matrices and capture metadata.
4. Only then install/fix compatible analysis-tool versions in separate environments and run a real dual-modality baseline.
5. GOD-16 requires three relevant scientific checks, a baseline/sensitivity study and supporting/conflicting reference evidence **after** a genuine dual-modality result exists.

GOD-15 and GOD-16 remain incomplete. GOD-71 has a documented candidate/method assessment, but actual data compatibility and runtime evaluation remain incomplete. GOD-20 cannot become a final release until software review, scientific QC/labels and upstream processing gates pass. There is no formal v0.1.0 tag from this work.
