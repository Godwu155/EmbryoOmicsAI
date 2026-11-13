# Modality contracts v1 — GOD-13 (2026-11-13)

RNA ingestion is implemented. Protein and ATAC contracts below define the next supported interfaces; their measurement readers and integration are **not implemented**. Metadata validation does not make them runnable.

| Modality / assay | Observation | Feature and value | Input contract | Missingness / normalization |
| --- | --- | --- | --- | --- |
| RNA / `rna_10x` | Cell within a declared capture sample | Gene ID; raw nonnegative integer UMI counts | Sparse 10x MTX/features/barcodes or H5AD with explicitly selected raw count layer; `gene_symbol` separate from unique gene ID | Zero is an observed zero; no imputation. Preserve counts; normalize/log1p only in working X. |
| Protein / `mass_spec_quantification` | Biological sample, potentially a documented embryo pool | Protein accession or explicit protein group; abundance with named unit | Long CSV: sample_id, protein_id, id_namespace, abundance, unit, missing_reason; companion sample metadata and source search/quantification method | Blank means missing, never silently zero. Record LFQ/intensity/other unit and whether already log transformed. Choose normalization only after reviewing method and missingness. |
| ATAC / `atac_fragments` | Cell within a capture sample | Genomic interval; unique fragment count | `fragments.tsv.gz` + index, or sparse peak-by-cell counts plus peaks and barcodes; exact genome build, coordinate convention and peak-set version required | Zero means no observed fragment in the specified interval. No RNA normalization on ATAC. QC and feature processing remain separate. |

Protein example: `examples/modality_contracts/protein_quantification.csv` is explicitly synthetic and is not an accepted biological dataset. ATAC peak/fragment conventions are documented here rather than supplying fabricated experimental files. CITE-seq antibody counts are not mass-spectrometry protein abundance; methylation requires its own coverage and methylated-count contract and is outside this version.

## IDs and objects

RNA uses AnnData. Each other modality will retain its own object/table until observation identity and measurement units are validated. MuData may be used for verified shared cells, but it cannot create pairing. Sample-level mass-spectrometry tables remain separate from per-cell RNA/ATAC; stage/organism/tissue comparisons are explicitly aggregated and conditional.

Preserve original gene, protein-group and interval identifiers, their namespaces, reference build and annotation release. Map proteins through an explicitly versioned accession→gene table with one-to-many flags; do not choose a gene silently for a protein group. RNA symbols are not unique IDs. Peaks map to a genome/interval first; nearest-gene association is a hypothesis, not an identity mapping. Missing reference versions prevent verified genome/feature linkage.

Reject unknown abundance units, mixed coordinate builds, duplicate primary IDs, nonfinite/negative raw counts, implicit cell↔sample joins, unsupported assays and missing required identity evidence. Report ambiguous/unmatched features instead of dropping them silently. Never broadcast a protein abundance value onto all cells.

## Sample and pairing metadata — GOD-14

Required sample columns: `record_id,dataset_id,sample_id,species,stage,modality,assay,measurement_unit`. `stage` retains the original source label; known-source descriptions live in `docs/developmental_stages.yaml`. No cross-species stage conversion or automatic stage prediction is performed. Optional `embryo_id,cell_id,replicate_id` must come from the source and remain empty when unknown. Cell-resolved links need one metadata row per linked cell; sample-only records may omit cell_id.

`record_id` is unique; `(dataset_id,sample_id,modality,cell_id)` is the composite identity when cell IDs are present. Identical short sample/cell names from different datasets are distinct. Pair CSV columns: `left_record,right_record,pairing_level,evidence_url`.

| Level | Required evidence/consistency | Permitted interpretation |
| --- | --- | --- |
| same_cell | Same dataset/sample, identical nonempty cell_id, cell measurement units, species/stage agreement, explicit reference | Candidate cell-resolved multimodal link, still requires chemistry/barcode/genome audit |
| same_embryo | Same dataset, matching known embryo_id, species/stage and declared replicate agreement, reference | Embryo-level link; not necessarily the same cell |
| same_sample | Same dataset and documented sample_id, compatible species/stage/embryo/replicate fields, reference | Sample-level association only |
| stage_only | Matching known species/stage and reference; can span studies | Stage comparison only; tissue, pooling, batch and replication remain possible confounders |
| unpaired / unknown | Explicitly declared | No pairing inferred |

Validation command:

```powershell
.venv/Scripts/python.exe -m embryoomics.identity --samples examples/modality_contracts/samples.csv --pairs examples/modality_contracts/pairs.csv --out results/sample_mapping_demo
```

Outputs: `sample_mapping.csv`, `unmatched_samples.csv`, `pairing_summary.json`. The validator checks declarations and rejects contradictions; it does not verify a paper merely because its URL was supplied, infer missing IDs, run integration, or assign a developmental stage. The supplied example is synthetic `stage_only` data and cannot establish real biological pairing.
