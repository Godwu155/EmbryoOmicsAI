# Development decisions (v0.1)

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
