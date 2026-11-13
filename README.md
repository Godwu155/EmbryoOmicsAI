# EmbryoOmics AI v0.1

Development priorities and acceptance gates: [next steps](docs/NEXT_STEPS.md).

Local, auditable mouse embryo single-cell RNA analysis from sparse raw counts to QC, Scanpy clustering, marker evidence, human review, and an HTML report. It does not assign definitive cell types, infer a developmental trajectory, or require an API key.

## Setup (PowerShell)

Use Python 3.11 or 3.12. From this repository:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
```

If Python 3.11 is unavailable, use `py -3.12 -m venv .venv`. Pinned direct dependencies are in `pyproject.toml`; `requirements-lock.txt` records the complete installed Python 3.12 environment (`python -m pip install -r requirements-lock.txt` for exact replay). Each run records the installed versions of the main analysis packages. Output, local data, virtual environments, and credentials are Git-ignored.

## Get the verified GEO sample

The [GSM8559287 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8559287) describes an E7.0 mouse embryo sample from [GSE278981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278981). The inspected record did not establish a specific reuse license, so the example records `data_license: unknown`; verify terms before redistribution. Download these three supplementary files into `data/GSM8559287/`, keeping the filenames unchanged:

1. [barcodes TSV](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8559nnn/GSM8559287/suppl/GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_barcodes.tsv.gz)
2. [features TSV](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8559nnn/GSM8559287/suppl/GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_features.tsv.gz)
3. [sparse matrix](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8559nnn/GSM8559287/suppl/GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_matrix.mtx.gz)

PowerShell download commands, when network access is available:

```powershell
New-Item -ItemType Directory -Force data/GSM8559287 | Out-Null
$base = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8559nnn/GSM8559287/suppl/'
$prefix = 'GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_'
foreach ($name in @('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz')) {
  Invoke-WebRequest -Uri ($base + $prefix + $name) -OutFile (Join-Path 'data/GSM8559287' ($prefix + $name))
}
```

Do not unpack the `.gz` files. The loader recognizes this unique common prefix. The sample ID and stage in `examples/GSM8559287.yaml` are taken from GEO; `genome_build: unknown` is deliberately retained. Validate the files and review QC before editing thresholds:

```powershell
embryoomics validate --input data/GSM8559287 --config examples/GSM8559287.yaml
embryoomics run --input data/GSM8559287 --config examples/GSM8559287.yaml --out results/GSM8559287_001 --confirm-qc
streamlit run app.py
```

For a `.h5ad`, set `raw_count_source` to `X`, `raw/X`, or `layers/<name>` in the YAML. The selected matrix must be sparse, finite, nonnegative integer counts. Optional metadata CSV needs exact `cell_id` matches; `sample_id` and `stage` can be supplied by the single-sample config. `embryo_id` and `batch` are never guessed. Optional reviewed annotation CSV can be passed with `--annotations`; use `unknown` when evidence is insufficient. A repeated output path is rejected unless `--overwrite` is explicit; overwrite is restricted to an existing EmbryoOmics run directory inside the current working directory.

## Output and review

Optional offline functional enrichment is available through `python -m embryoomics.enrichment` (or `embryoomics-enrich` after installation), with a local GMT, explicit species/ID/source/version/license metadata and a detected-gene background. See [enrichment usage and limits](docs/ENRICHMENT_ACCEPTANCE.md). It saves a separate report/ZIP and does not modify the reviewed RNA run.

Sample/pairing metadata contracts and the validator are documented in [modality_schema.md](docs/modality_schema.md). Protein and ATAC measurement readers remain planned. The [Backlog audit](docs/BACKLOG_AUDIT.md) records implemented evidence, naming differences and remaining gates for the original Linear tasks.

New RNA runs include `run_status.json` and `run.log`. Output-writing failures are marked failed and the interface refuses to display them as completed results; the offline ZIP is promoted from a temporary file only after packaging finishes. Failures before an output directory exists remain visible as CLI/UI diagnostics. Existing output directories are never modified by a rejected run. `embryoomics --version` reports the installed package version, and CI includes basic Ruff syntax/undefined-name checks.

`mitochondrial_gene_ids` may be an explicit list of unique input gene IDs in the YAML. When supplied, it overrides symbol-prefix detection and all listed IDs must exist. Otherwise the documented case-insensitive `^mt-` symbol rule remains the default. The detected method is recorded in the manifest and shown in the interface.

### Optional Scrublet review

Scrublet is off by default. After input validation, the QC page offers **Score only** or **Filter after explicit confirmation**. Enter an expected doublet fraction and a score threshold explicitly, then calculate the preview. Review observed/simulated score distributions and the cell audit CSV. Score-only keeps candidates; filtering requires a separate confirmation as well as the usual QC confirmation. Changes to input, QC or Scrublet settings invalidate the preview and confirmations. The run recomputes scores and rejects a mismatch with the reviewed preview. Source YAML is preserved.

The method calls pinned Scanpy on raw counts after ordinary QC, uses exact neighbors and variance-normalized, uncentered SVD (`mean_center=False`) to avoid centering the full matrix. This differs from Scanpy's centered default; all parameters are recorded. Analyze one droplet capture batch per run. Tiny data, multiple declared batches, internal unscored cells, or unsuitable PCA dimensions block scoring with instructions; turn Scrublet off when inappropriate. The GEO record does not establish upstream doublet removal for GSM8559287. Candidate labels are exploratory and cannot prove that a cell is a doublet.

For CLI use, add a `doublets` mapping to a copy of your configuration: `mode: score` or `filter`, `expected_doublet_rate` and `threshold` (both explicitly between 0 and 1), `n_prin_comps: 30`, and `upstream_handling: unknown` (or `not_removed` / `removed` when documented). Values depend on the experiment; no scientific threshold is supplied. Preview and review before requesting filtering:

```powershell
embryoomics preview-doublets --input data/GSM8559287 --config path/to/reviewed.yaml --out results/doublet_preview
embryoomics run --input data/GSM8559287 --config path/to/reviewed.yaml --out results/doublet_run --confirm-qc --confirm-doublets
```

The extra `--confirm-doublets` flag is required only for filter mode. Preview output records settings and input identity. Analysis exports `doublet_scores.csv` for **all input cells** (QC-failed cells have blank scores), `doublet_simulated_scores.csv`, and `figures/doublet_scores.png`. The audit distinguishes ordinary QC, candidate status, actual removal and final retention. Retained-cell scores/flags are stored in `processed.h5ad`; the manifest and bilingual report record the method, decision and counts. ZIP and annotation-report refresh include these artifacts.

New runs reuse the pinned Scanpy `dotplot` and `matrixplot` APIs to produce marker dot plots and cluster-mean expression heatmaps. The Results page displays these all-cluster overview plots in expandable panels; the bilingual HTML and offline ZIP include them too. The panel takes the first three positive-logfoldchange markers per cluster in `markers.csv` rank order, deduplicates gene IDs and caps the panel at 30 genes. `marker_plot_genes.csv` and the manifest record the panel and method. Colors show mean log1p normalized expression without scaling; dot size shows the fraction of cells with expression above zero. Labels retain gene IDs, including duplicate symbols. Empty positive-marker panels are recorded and skipped. Older runs remain readable and need a new analysis run to obtain these plots.

Open-source reuse decisions and follow-up tasks: [reuse plan](docs/OPEN_SOURCE_REUSE.md).

Optional CELLxGENE browsing is verified on native Windows in a separate pinned environment. Export a derived view with `.venv/Scripts/python.exe -m embryoomics.cellxgene export --run results/<run> --out results/<new_view>`, then follow the [viewer setup and annotation exchange guide](docs/CELLXGENE_ACCEPTANCE.md). The bridge preserves gene IDs and existing UMAP, converts cell labels into cluster proposals with conflict records, and applies reviewed proposals only with explicit confirmation through the normal annotation-history/report service. CELLxGENE dependencies are not installed into the main analysis environment.

`results/<run>/` contains `report.html` (Data, QC, Results, Review), `run_manifest.json`, `config.yaml`, `qc_metrics.csv`, `qc_preview.csv`, `markers.csv`, `cluster_composition.csv`, `annotations.csv`, `processed.h5ad`, and `figures/`. The manifest includes input file SHA256 hashes, dimensions, thresholds, retained cells, seed, versions, runtime, peak working set on Windows, warnings, and file list. `processed.h5ad` retains sparse raw integer counts in `layers['counts']`; `X` is normalized/log transformed. Marker tests are exploratory cluster-vs-rest Wilcoxon comparisons on log normalized values, with expression fractions derived from raw counts. The report links to its machine readable artifacts.

New runs also produce `report_bundle.zip` containing the report and every linked artifact. Download and extract the ZIP before opening `report.html` to preserve the offline images and file links. In the interface, changes to input files or the source analysis configuration invalidate the earlier QC preview and its confirmation; validate again before running.

After validating input, use the **QC** page to enable or disable minimum/maximum detected-gene and mitochondrial-percentage filters and enter values. The preview updates from the validated per-cell metrics without reloading the count matrix. It shows retained/removed cells and rule-level counts; overlapping rule failures are distinguished from cells newly removed at each step. Changing a threshold cancels the earlier confirmation. Missing values, invalid ranges, and fewer than four retained cells block running; unavailable mitochondrial metrics disable that filter. These controls do not supply recommended biological thresholds.

Confirm the current preview and run into a new results folder. The actual UI-selected values are saved in that run's `config.yaml`, manifest, and `qc_filter_reasons.csv`; the source YAML is preserved. **Download current run configuration** exports the full effective YAML for reuse. CLI runs continue to use the values in their supplied YAML.

The Streamlit interface has Data, QC, Results, and Report tabs. Select **English** or **简体中文** in the sidebar; paths and the current session remain available when switching. Validate input, inspect QC, explicitly confirm thresholds, then run. The exported `report.html` has its own English/简体中文 selector and remembers the browser's choice. Language changes affect presentation only: CSV/JSON field names, gene IDs, user-entered annotations, parameters, and checksums remain unchanged. No AI text is generated in v0.1.

To review existing results, set the sidebar results folder to that run. On **Results**, edit candidate cell types, supporting/conflicting genes, references, review status, and notes, then click **Save annotations and refresh report**. The app saves the CSV, refreshes the HTML and offline ZIP, and records timestamped before/after revisions under `annotation_history/`. Cluster IDs are read-only; editing annotations does not rerun analysis or change the expression matrix or embedding. Changes from another page require reloading before saving. CSV export and CLI `--annotations` import remain available for reuse in other runs.

## Repository layout

The Results page includes a cluster selector that filters the UMAP, marker table, and annotation editor together. Search a gene symbol or ID, then choose a specific gene ID to view log1p normalized expression and detection statistics across clusters. Duplicate symbols retain separate gene IDs. Single-gene reads are chunked rather than loading the entire expression matrix. Save annotation changes before switching clusters; saving one cluster preserves the others, and CSV export includes all clusters.

GitHub validation uses a fresh Windows/Python 3.12 runner to install locked dependencies, run synthetic/interface tests, and build a development wheel. Scientific review and real-data acceptance remain separate from CI.

```text
embryoomics/  # shared config, import, QC, Scanpy, markers, report and CLI
app.py        # local Streamlit interface
examples/     # GEO sample configuration
tests/        # synthetic contract and end-to-end tests
docs/         # PRD, SPEC, decisions
data/         # local inputs; ignored by Git
results/      # local outputs; ignored by Git
```

## Limits

### Experimental human blastocyst path

The separately verified GSE66507 supplementary matrix can now be imported as a declared gene-by-cell TSV and run through the RNA pipeline:

```powershell
python validate_gse66507.py --output data/GSE66507_counts.txt.gz
embryoomics validate --input data/GSE66507_counts.txt.gz --config examples/GSE66507.yaml
embryoomics run --input data/GSE66507_counts.txt.gz --config examples/GSE66507.yaml --out results/GSE66507_trial --confirm-qc
```

This is a software smoke test, not a reviewed human embryo analysis. The 30-cell file supplies Ensembl IDs without gene symbols, so mitochondrial percentage is unavailable and no mitochondrial filter is applied. No embryo IDs, biological replicates, or cell type annotations are inferred from the column names. Clusters remain exploratory until the source metadata and biological evidence are reviewed. The mouse example remains the documented v0.1 benchmark.

The primary v0.1 benchmark is the mouse E7.0 sample by research scope decision. Its larger cell count supports a more useful first RNA workflow test than this 30-cell human file, but a single stage and sample do not represent the full course of mouse embryonic development.

The stage labels and their provenance are recorded in `docs/developmental_stages.yaml`. E7.25 is a candidate for later validation; the E7.5 sample is posterior-region tissue and must not be silently treated as whole-embryo data. Unknown stage labels are retained verbatim and require review before stage comparison. The current single-sample pipeline does not perform cross-stage analysis.

- v0.1 covers one mouse RNA sample per run. No batch integration, FASTQ processing, RNA–ATAC pairing, mass spectrometry import, causal or lineage inference.
- GEO provides a filtered count matrix. The app does not recover cells removed upstream, and no independent embryo ID is fabricated.
- The `^mt-` gene symbol rule is reported; if no symbols match, mitochondrial percentage is unavailable and a mitochondrial threshold is rejected.
- A single sample cannot establish biological replicate variation. Marker p-values should not be read as embryo-level differential expression evidence.
- UMAP is an embedding, not a developmental time axis. Unknown cell types remain `unknown` until human review.

Synthetic test output is for software verification only. A real GSM8559287 result is reported only after the exact files have been obtained and an end-to-end command has completed successfully.

### Verified local example (2026-09-27)

The three GEO files above passed validation as a 10,232-cell × 32,285-gene sparse integer matrix; 13 mitochondrial symbols matched the documented rule. With the example's deliberately unset QC thresholds, 10,232 cells were retained. The end-to-end run produced nine **unreviewed clusters**, not nine confirmed cell types. Its manifest records 53.33 seconds of analysis time, a 1,887,752,192-byte peak Windows working set (about 1.76 GiB), and SHA256 checksums for all three inputs. These figures describe this local machine and these exact files; they are not biological conclusions. The output is in `results/GSM8559287_001/` locally and remains Git-ignored.
