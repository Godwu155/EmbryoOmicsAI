# Linear Backlog audit — 2026-11-13

Issues were read against their complete original acceptance criteria. States below distinguish implementation review from scientific acceptance; no formal release or PR merge is implied. Current package/CLI is `embryoomics`, not the early planning placeholder `emap`; current root package layout and output filenames follow the accepted SPEC and README.

| Issue | Evidence and work completed | Disposition / remaining gate |
| --- | --- | --- |
| GOD-8 | PRD, SPEC, scope decisions and source-stage boundaries already tracked | In Review; scientific scope review remains separate |
| GOD-9 | Locked Windows environment, packaging, tests and Actions existed. Added `--version`, pinned Ruff basic syntax/undefined-name checks and CI step | In Review; existing Actions run 43371258623 succeeded; verify latest CI separately |
| GOD-10 | `io.py`, synthetic validation tests, real 10x and H5AD `layers/counts` loading | In Review; real H5AD 10,232 × 32,285, counts preserved |
| GOD-11 | QC previews, reasons and optional Scrublet existed. Added explicit mitochondrial gene-ID override and missing-ID rejection | In Review; official biological thresholds/upstream processing remain unresolved |
| GOD-12 | RNA normalization/HVG/PCA/neighbors/Leiden/UMAP, saved seed/parameters and real-data reproducibility | In Review; original output names map to `processed.h5ad` and `figures/umap.png` (cluster/stage/sample panels), no confirmed cell identities |
| GOD-13 | Added `modality_schema.md`, observation/feature/unit/missingness/ID contracts, synthetic input examples and refusal rules | In Review for specification only; protein/ATAC measurement readers are not implemented |
| GOD-14 | Added declared-identity/pairing validator with exact namespace/unit/evidence checks, negative tests and mapping/unmatched/summary outputs | In Review for metadata validation; URLs and declared IDs are not independently verified biological pairing |
| GOD-17 | YAML/CLI pipeline existed. Added output-writing/completed/failed status and JSONL log, temporary ZIP promotion and incomplete-result UI guard | In Review; failures before creating an output directory retain CLI/UI diagnostics rather than fabricating a successful run directory |
| GOD-18 | AI disabled and scoped aggregate payload exists | Keep Backlog. Static fallback/10-case evaluation and provider-neutral client remain unimplemented; no provider, credentials or model assumed |
| GOD-19 | Generated bilingual offline reports, linked artifacts, annotation-refresh tests and ZIP checks | In Review; no invented multiomics or AI sections |
| GOD-20 | Real replay, locked environments, development wheel and documentation exist | In Progress; source processing, final QC/cell-label review, software review/merge and release gate remain open |
| GOD-15 | Added explicit feasibility and data-acquisition gates | Keep Backlog with dependencies; no verified dual-modality input/result. Next acquire candidate quantitative files and replicate metadata |
| GOD-16 | RNA QC/Scrublet sensitivity results exist; no genuine integrated result | Keep Backlog, blocked by GOD-15; cannot claim biological multiomics validation |
| GOD-71 | Recorded separate scvi-tools/AlphaPeptStats assumptions, license candidates and data/runtime gates | In Progress for assessment; installation and real input compatibility remain incomplete |
| GOD-70 | Added offline GSEApy integration, provenance/background/ID audit and global correction; real direct-GO example | In Review; functional associations remain exploratory |

Specifications and planned filenames are not silently treated as implemented capabilities. Current successful runs retain normalized X and raw `counts` in `processed.h5ad`; QC summaries are `qc_preview.csv`, `qc_metrics.csv`, `qc_filter_reasons.csv`. There is no separate materialized `qc_filtered.h5ad`; it is the retained-count layer/observation set of the processed object. These equivalences are documented for reviewers rather than renaming existing artifacts.

See `ENRICHMENT_ACCEPTANCE.md`, `modality_schema.md`, `MULTIOMICS_FEASIBILITY.md`, `SCRUBLET_ACCEPTANCE.md`, `CELLXGENE_ACCEPTANCE.md` and `RELEASE_CHECKLIST.md` for the evidence and limits. Backlog issues with genuine unmet scientific/data requirements remain open; they were not mass-closed.

Final local validation: 69 tests pass in both development and independent locked environments; Ruff basic checks and dependency checks pass. Real GSM8559287 replay produces 10,232 cells and 9 unreviewed clusters, raw counts identical to the previous baseline, ARI=1.0, completed state and matching ZIP status/files. The new wheel builds successfully (SHA256 `5b4a019ebe893531288395c7985d942de9ff0959ffc1a9ecefc9284f528f3ebc`). Evidence remains in Git-ignored `results/GSM8559287_backlog_acceptance_20261113/` and `results/test-temp/`.
