# Development decisions (v0.1)

## 2026-11-13: interactive QC configuration

The QC page computes threshold previews from validated per-cell metrics, with explicit switches for null/disabled rules. UI-selected thresholds are passed as a validated QC override to the shared service; CLI YAML behavior is retained. The source YAML is preserved, while the run snapshot and manifest record effective values. Filtering uses the same metric-based function in preview and execution, exports overlapping and sequential removal counts, and invalidates confirmation after parameter changes. A stale input/source configuration is rejected before analysis. The real GSM8559287 UI acceptance run retained 8,813 cells in both preview and execution.

## 2026-11-12: direct annotation editing

The Results page now saves candidate annotation edits directly into the selected run. Saving regenerates the HTML report and offline ZIP, with timestamped before/after snapshots under `annotation_history/`. Cluster IDs are read-only; the computational clustering, marker table, expression matrix and embedding are preserved. Stale pages must reload when another page has saved changes. Artifacts are prepared before commit, and ordinary replacement failures restore prior files. The CSV download/import path remains available for exchanging annotations between runs.

## 2026-09-27: workspace and source documents

The supplied working directory `D:\Single-cell analysis` did not exist when work began. The existing PRD and SPEC were found in `D:\EmbryoOmicsAI\docs`, a directory without Git metadata; they were copied into the requested writable workspace. A second, root-level `D:\EmbryoOmicsAI\EmbryoOmicsAI_PRD.md` differs from `docs/EmbryoOmicsAI_PRD.md`. The user explicitly asked to read the two files under `docs/`, so that PRD and the paired SPEC are authoritative for v0.1. The root-level copy was not silently substituted.

## 2026-09-27: GEO distribution vs generic 10x names

The [GSE278981 series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278981) lists three samples and a 246.7 MB `GSE278981_RAW.tar`. The [GSM8559287 sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8559287) identifies an E7.0 mouse embryo sample and separately publishes these exact files:

- `GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_barcodes.tsv.gz`
- `GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_features.tsv.gz`
- `GSM8559287_cellranger_Pmem-E7-0_outs_filtered_feature_bc_matrix_matrix.mtx.gz`

The SPEC implies names exactly `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. The importer accepts the same unique shared prefix on the three filenames and rejects ambiguous or mixed triples. No sample label or embryo ID is inferred from the barcode strings. `sample_id: GSM8559287` and `stage: E7.0` come from the GEO record; genome build remains `unknown` because the record's "Assembly: ARC mice" does not establish a standard reference genome build. This sample is a filtered Cell Ranger feature-barcode matrix, not FASTQ or an independent embryo replication series.

## 2026-09-27: local Python and CSV scope

The SPEC targets Python 3.11. The available bundled runtime is Python 3.12.14, so the package supports 3.11–3.12 and records the exact runtime in each manifest. No compatible 3.11 interpreter was found initially. CSV count matrix import is not implemented because it is optional in the PRD and is unsafe for nontrivial sparse GEO matrices; a descriptive format error is returned. Metadata CSV remains supported.

## 2026-09-27: AI boundary

No LLM provider, citation retrieval, or calibration protocol is specified. v0.1 exposes the allowed aggregate payload shape and a disabled status only. It does not make a live model call, invent literature interpretations, or imply a confidence score. Analysis, reports, and review work without an API key.

## 2026-09-27: downloaded sample validation

The user supplied all three GSM8559287 files. Their measured sizes are 51,872, 290,896, and 101,563,067 bytes for barcodes, features, and matrix respectively. The importer verified a 10,232-cell × 32,285-gene nonnegative integer sparse count matrix, 13 `^mt-` symbols, and unique IDs. No QC thresholds were guessed: the example retains all 10,232 cells and records null thresholds. The successful Scanpy run produced nine unreviewed clusters and all specified report artifacts. Its three SHA256 values are retained in the local run manifest, while the large data and output files are ignored by Git. The original file naming mismatch is resolved in code, without renaming or fabricating metadata.

The inspected GEO sample record does not supply a verified redistribution license. The example config and report therefore state `data_license: unknown` instead of inferring reuse terms from public download availability.

## 2026-11-11: first formal benchmark species

User review confirmed: validate the first release with the mouse E7.0 single sample; extend to other stages later; retain the posterior-region annotation for E7.5.

The user chose mouse as the first formal research benchmark after comparing the two files actually tested. GSM8559287 has 10,232 cells and an already completed local RNA run; GSE66507 has 30 human blastocyst cells and is useful as an import and pipeline smoke test. These datasets cover different developmental stages, so the comparison does not establish that mouse embryo data are generally more complete. The primary v0.1 scientific scope remains the mouse E7.0 sample. Human input support remains experimental until its source metadata, QC assumptions, and biological interpretation are reviewed. One mouse sample does not establish a developmental time series or biological replicate variation.
