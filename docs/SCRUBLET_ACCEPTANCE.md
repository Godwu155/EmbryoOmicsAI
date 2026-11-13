# Scrublet software acceptance — 2026-11-13

Scope: optional single-capture droplet RNA scoring via pinned Scanpy 1.11.5. Default remains off. No new dependencies were introduced. `mean_center=False` explicitly selects uncentered variance-normalized SVD; this is not an attempt to replicate the upstream centered default. Exact nearest neighbors avoid requiring an approximate-neighbor dependency; an explicit threshold avoids requiring scikit-image's automatic threshold selection.

## Source evidence

[GSM8559287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8559287) describes a 10x channel per sample and refers processing to Aryamanesh 2022. It does not establish doublet-removal status. The app records `upstream_handling: unknown`; this is not evidence that the authors omitted filtering. Candidate scores are not ground-truth labels.

Method reference: [Scanpy Scrublet API](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html). Our installed, pinned 1.11.5 implementation was inspected before integration. The wrapper accepts raw counts, can internally remove low-information cells/genes, and exposes observed/simulated scores. Nonfinite or unscored observed cells are rejected rather than assigned invented values. Multiple declared capture batches and tiny datasets are rejected with guidance.

## Real-data sensitivity experiment

Reproducible command: `python scripts/validate_scrublet.py --out results/<new_directory>` from the repository root, with the previously validated `results/GSM8559287_scanpy_plots_20261113/` baseline available. The script refuses an existing output directory.

Experimental settings: expected fraction **0.05**, threshold **0.15**, 30 SVD components, seed 42, ordinary QC thresholds left unset as in the baseline. These numbers test software behavior and are **not recommended biological thresholds**. The original data and configuration were preserved.

| Measure | Result |
| --- | --- |
| Score-only retained cells | 10,232 |
| Score-only vs baseline cluster ARI | 1.0 |
| Candidate cells above test threshold | 584 |
| Explicit filter retained cells | 9,648 |
| Clusters, score-only → filter | 9 → 8 (all unreviewed) |
| Filter vs score cluster ARI on shared cells | 0.906244 |
| Unique genes among per-cluster top-10 markers | 83 before, 75 after; 74 shared |
| Preview vs execution retained cells | Matched |
| Retained raw counts vs baseline | Identical |
| Score-only vs filter observed scores | Matched |
| Both ZIPs and manifest artifact membership | Passed |

Cluster IDs are not aligned labels; ARI compares partitions. Marker overlap is descriptive, not a validation of cell identity. Filtering changes downstream clustering and therefore still needs scientific review. Both observed/simulated histogram panels and threshold legends were visually checked.

Local output: `results/GSM8559287_scrublet_acceptance_20261113_v2/`, including `acceptance.json` and separate score/filter runs. These remain Git-ignored. An earlier partial software run failed report rendering before a translation-key fix and is not acceptance evidence.

## Automated checks

45 tests pass in both the development and existing independent locked environment, covering actual Scrublet reproducibility/raw-count isolation, invalid parameters, multiple batches/tiny input rejection, QC-failed cells with blank scores, filter confirmation and stale-digest rejection, score-only retention, preserved counts, bilingual reports, annotation refresh, ZIPs, CLI preview and confirmation, and the full UI review flow including language changes and threshold invalidation. `pip check` passes. The development wheel was built successfully with isolated build dependencies. Wheel SHA256: `08a937dd41a350e9a55801614cdb0f7176921a1ce348c70486dea579c6ed3341`.

PR and Linear links are maintained in [the reuse plan](OPEN_SOURCE_REUSE.md). Scientific validation of upstream processing and final QC/cell labels remains separate from software acceptance.
