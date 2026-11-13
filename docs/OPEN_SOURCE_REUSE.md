# Open-source reuse plan — 2026-11-13

Keep the local Streamlit workflow, sample/stage identity checks, explicit QC confirmation, human annotation review, and auditable reports. Reuse maintained analysis APIs rather than implementing their algorithms again. Repository license labels below describe the inspected upstream projects; verify the exact selected release and its dependencies before adding a new package. No third-party source has been copied in this batch.

## Implemented: GOD-67

Use the existing pinned Scanpy 1.11.5 dependency (BSD-3-Clause): `sc.pl.dotplot` and `sc.pl.matrixplot`. A bounded gene subset is copied before plotting; the full expression matrix remains sparse. Both plots explicitly use `X` with `use_raw=False`, without scaling. Marker selection, gene IDs, expression semantics and artifacts are saved in `run_manifest.json` and `marker_plot_genes.csv`. All-cluster plots appear in expandable Results panels and the bilingual report/ZIP; they do not follow the individual cluster filter. Existing runs without these artifacts still render normally.

Selection: first three positive-logfoldchange entries per cluster in the existing marker table's rank order, unique gene IDs, maximum 30 genes. This is an exploratory display panel, not a curated biological marker reference or a significance filter. Empty panels are explicitly skipped. Gene labels include IDs so duplicate symbols stay distinct.

Validation: 31 tests passed, including actual Scanpy dot statistics against known means/detection fractions despite a conflicting raw/count layer, duplicate symbols, empty panels, counts preservation, bilingual report links and ZIP inclusion. A new full GSM8559287 run retained 10,232 cells and produced 9 unreviewed clusters in 55.2 seconds; 25 unique marker genes were plotted. Both PNGs were visually inspected. Output: `results/GSM8559287_scanpy_plots_20261113/` (Git-ignored). Scientific QC and cell-type acceptance remain pending.

## Follow-up tasks

GOD-68 now reuses the existing Scanpy Scrublet API without new dependencies. Default is off; scoring and filtering are separate decisions, with explicit rate/threshold and per-run seed. UI preview changes invalidate confirmations and a score digest binds execution to the preview. CLI requires an extra `--confirm-doublets` for filtering. Raw post-QC counts are copied into a fresh AnnData, exact neighbors are used, and `mean_center=False` selects variance-normalized uncentered SVD; this is explicitly different from Scanpy's centered default. Batch, size and unscored-cell checks block unsupported cases. Per-cell audit includes QC-failed cells without invented scores; observed/simulated distributions and decisions survive report refresh/ZIP export. Upstream processing remains unknown for GSM8559287.

| Linear | Upstream | Proposed reuse and acceptance gate |
| --- | --- | --- |
| [GOD-68](https://linear.app/godwu155/issue/GOD-68) | [Scanpy Scrublet](https://github.com/scverse/scanpy), BSD-3-Clause | Review upstream doublet handling; validate count input, dependencies, parameters and score distributions before any optional filtering. |
| [GOD-69](https://linear.app/godwu155/issue/GOD-69) | [CELLxGENE](https://github.com/chanzuckerberg/cellxgene), MIT | Test existing H5AD, embeddings, gene IDs and annotation exchange in an isolated environment; record Windows deployment and resources. |
| [GOD-70](https://linear.app/godwu155/issue/GOD-70) | [GSEApy](https://github.com/zqfang/GSEApy), BSD-3-Clause | Validate mouse IDs, tested-gene background, gene-set provenance/version/terms and multiple-testing correction. |
| [GOD-71](https://linear.app/godwu155/issue/GOD-71) | [scvi-tools](https://github.com/scverse/scvi-tools), BSD-3-Clause; [AlphaPeptStats](https://github.com/MannLabs/alphapeptstats), Apache-2.0 | Evaluate separately after obtaining real batch/replicate designs or proteomic quantification tables. Do not invent cell-level RNA/protein pairing. |

[OmicVerse](https://github.com/omicverse/omicverse) (GPL-3.0) remains a broader module candidate; it has not been installed or integrated. [CellTypist](https://github.com/Teichlab/celltypist) (MIT) needs an embryo-appropriate model/reference before adoption. Software, model and dataset terms must be evaluated separately.
