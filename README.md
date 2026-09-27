# EmbryoOmics AI v0.1

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

`results/<run>/` contains `report.html` (Data, QC, Results, Review), `run_manifest.json`, `config.yaml`, `qc_metrics.csv`, `qc_preview.csv`, `markers.csv`, `cluster_composition.csv`, `annotations.csv`, `processed.h5ad`, and `figures/`. The manifest includes input file SHA256 hashes, dimensions, thresholds, retained cells, seed, versions, runtime, peak working set on Windows, warnings, and file list. `processed.h5ad` retains sparse raw integer counts in `layers['counts']`; `X` is normalized/log transformed. Marker tests are exploratory cluster-vs-rest Wilcoxon comparisons on log normalized values, with expression fractions derived from raw counts. The report links to its machine readable artifacts.

The Streamlit interface has Data, QC, Results, and Report tabs. Select **English** or **简体中文** in the sidebar; paths and the current session remain available when switching. Validate input, inspect QC, explicitly confirm thresholds, then run. The exported `report.html` has its own English/简体中文 selector and remembers the browser's choice. Language changes affect presentation only: CSV/JSON field names, gene IDs, user-entered annotations, parameters, and checksums remain unchanged. Export annotation edits and rerun with the exported CSV to preserve review history in a new output directory. No AI text is generated in v0.1.

## Repository layout

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

- v0.1 covers one mouse RNA sample per run. No batch integration, FASTQ processing, RNA–ATAC pairing, mass spectrometry import, causal or lineage inference.
- GEO provides a filtered count matrix. The app does not recover cells removed upstream, and no independent embryo ID is fabricated.
- The `^mt-` gene symbol rule is reported; if no symbols match, mitochondrial percentage is unavailable and a mitochondrial threshold is rejected.
- A single sample cannot establish biological replicate variation. Marker p-values should not be read as embryo-level differential expression evidence.
- UMAP is an embedding, not a developmental time axis. Unknown cell types remain `unknown` until human review.

Synthetic test output is for software verification only. A real GSM8559287 result is reported only after the exact files have been obtained and an end-to-end command has completed successfully.

### Verified local example (2026-09-27)

The three GEO files above passed validation as a 10,232-cell × 32,285-gene sparse integer matrix; 13 mitochondrial symbols matched the documented rule. With the example's deliberately unset QC thresholds, 10,232 cells were retained. The end-to-end run produced nine **unreviewed clusters**, not nine confirmed cell types. Its manifest records 53.33 seconds of analysis time, a 1,887,752,192-byte peak Windows working set (about 1.76 GiB), and SHA256 checksums for all three inputs. These figures describe this local machine and these exact files; they are not biological conclusions. The output is in `results/GSM8559287_001/` locally and remains Git-ignored.
