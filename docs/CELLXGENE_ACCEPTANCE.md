# CELLxGENE optional viewer — 2026-11-13

GOD-69 adds a verified Windows browsing/export path and a conservative annotation bridge. CELLxGENE remains optional and runs in a separate environment. Source analysis, embeddings, counts and annotations are preserved until a reviewer explicitly applies proposals through the existing review service.

## Installation and launch

Verified: Windows / Python 3.12, CELLxGENE 1.3.0. Upstream is [MIT licensed](https://github.com/chanzuckerberg/cellxgene); its dependencies retain their own terms. The complete tested environment is pinned in `requirements-cellxgene-windows-lock.txt`. In particular, CELLxGENE pins NumPy 2.0.1, while our main analysis environment pins NumPy 2.3.5. **Do not install this lock into the main analysis environment.**

From the repository root, in PowerShell:

```powershell
.venv/Scripts/python.exe -m venv results/cellxgene_env_20261113
results/cellxgene_env_20261113/Scripts/python.exe -m pip install -r requirements-cellxgene-windows-lock.txt
results/cellxgene_env_20261113/Scripts/python.exe -m pip check
.venv/Scripts/python.exe -m embryoomics.cellxgene export --run results/GSM8559287_scanpy_plots_20261113 --out results/my_cellxgene_view
./scripts/start_cellxgene.ps1 -Exchange results/my_cellxgene_view -Port 5057
```

Open the localhost URL printed by the server. Keep its terminal running; Ctrl+C stops it. The helper binds only to `127.0.0.1`, verifies the browsing-copy checksum, and disables on-demand differential expression and gene-set saving so the viewer does not introduce a second untracked analysis workflow. Native Windows launch was tested directly; Docker/WSL was not needed. If running inside a sandbox, serving localhost and browser access may require the host network execution permission.

The Python module can also be invoked as `embryoomics-cellxgene` after installing the updated package. The export command needs the main analysis environment; the viewer command uses the separate CELLxGENE environment.

## Browsing copy

The exporter reads the source H5AD in backed mode, copies normalized `X` into float32 CSC for browsing, preserves the existing `X_umap`, and exposes cluster, stage, sample and available batch/embryo/doublet metadata. No PCA, UMAP or clustering is recomputed. Raw counts and the `raw` slot are omitted from this derived view; the original `processed.h5ad` retains them.

Gene display names are `symbol [gene_id]`, with a gene-ID fallback for missing symbols. This preserves duplicate-symbol genes as distinct entries and allows symbol or ID search. `gene_map.csv` records the mapping. `exchange.json` records source run, source H5AD/annotation SHA256, browsing-file SHA256, dimensions and expression semantics. A new output directory is required; existing views are never overwritten.

## Annotation exchange and review

`cell_annotations.csv` is initialized from the run's candidate labels and contains exactly one row per cell. CELLxGENE edits this copy and saves its own comment header/backups. Do not point CELLxGENE at our per-cluster `annotations.csv`; the schemas differ.

After saving in CELLxGENE, use a second terminal:

```powershell
.venv/Scripts/python.exe -m embryoomics.cellxgene propose --exchange results/my_cellxgene_view --labels results/my_cellxgene_view/cell_annotations.csv --out results/my_cellxgene_proposals
```

Review `annotation_proposals.csv` and `annotation_conflicts.csv`. Only unanimous, nonempty, non-unknown labels can propose a new cluster candidate. Mixed or unresolved cell labels retain the previous candidate and produce a conflict record; no majority-vote label is imposed. Existing supporting/conflicting genes and references are retained, review status becomes `in_review`, and a reviewer note identifies the proposal source. Cell IDs must match exactly; duplicate/missing/extra cells and stale source/view checksums are rejected.

Edit the proposal CSV if necessary, add/check supporting evidence, then explicitly save the reviewed result:

```powershell
.venv/Scripts/python.exe -m embryoomics.cellxgene apply --proposals results/my_cellxgene_proposals --confirm-reviewed
```

Apply uses the existing review lock, stale annotation-version check, validation, before/after history and HTML/ZIP refresh. It does not rerun analysis. An apply changes the source annotation version, so export a new view for further editing. Propose alone never changes source annotations. This initial bridge preserves cluster-level review; it does not introduce new per-cell biological truth labels or split computational clusters.

## Acceptance evidence

Local input: `results/GSM8559287_scanpy_plots_20261113/processed.h5ad`. Browsing output: `results/GSM8559287_cellxgene_20261113/` (Git-ignored).

| Check | Result |
| --- | --- |
| Original H5AD direct launch | Loaded in native Chrome; full UMAP and cluster controls visible |
| Browsing-copy cells × genes | 10,232 × 32,285 |
| Cluster categories | 9, same cell assignments |
| Embedding | Source coordinates preserved in export; live API returns 10,232 × 2 finite values |
| Ttr symbol search | `Ttr [ENSMUSG00000061808]` selected and expression coloring visible |
| Live Ttr expression vs source | Maximum absolute error 0.0 |
| Duplicate-symbol handling | 78 genes share symbols; all display labels unique |
| Browser annotation save | `software_test_only` test label saved; 9 unanimous cluster proposals generated |
| Final browser labels | Restored to `unknown`; UI confirmed All saved |
| Source H5AD/annotation hashes | Unchanged after browser save/proposal round trip |
| Only editable biological field | `candidate_cell_type`; cluster/stage/sample remain read-only |
| Viewer memory after API/browser checks | Working set 476,467,200 bytes (~454 MiB); peak 486,748,160 bytes (~464 MiB) |
| Original direct-viewer memory | Working set ~771 MiB; peak ~786 MiB in this local test |

These measurements describe the tested machine/files, not a general memory guarantee. The browsing file is about 131 MiB and is not smaller than the compressed source in this example. Temporary baseline servers were stopped; the verified local browsing server was left available on port 5057. It must be relaunched if its process or terminal is closed.

`scripts/validate_cellxgene.py` performs read-only live-API comparisons with the source using the isolated CELLxGENE interpreter. Example:

```powershell
results/cellxgene_env_20261113/Scripts/python.exe scripts/validate_cellxgene.py --exchange results/GSM8559287_cellxgene_20261113 --url http://127.0.0.1:5057
```

API evidence: `api_acceptance.json`; process resources: `resources.json`. Browser test proposals in `results/GSM8559287_cellxgene_proposals_20261113/` contain an explicitly synthetic test label and must not be applied to biological results. Automated tests cover source preservation, duplicate/missing symbols, exact cell-ID matching, mixed-label conflicts, stale changes, and explicit apply with history/report refresh.

The final code passed all 53 tests in both the development and existing independent locked environments. Main and viewer dependency checks passed; the PowerShell helper passed syntax validation. The development wheel was built, installed into the independent environment, and its exchange CLI imported successfully from outside the source directory. Wheel SHA256: `a98a1d44a20f336208207c233d41ae695bbfe9ea7d6695103e493dd6ea079d52`.

Scientific sample-processing evidence, QC thresholds and final cell types still require review. Successful browsing and exchange do not validate biological labels.
