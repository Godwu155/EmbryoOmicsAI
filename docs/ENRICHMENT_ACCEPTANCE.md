# Offline enrichment — GOD-70 (2026-11-13)

The optional CLI reuses pinned GSEApy 1.3.2 for offline marker over-representation. It accepts only a local GMT plus explicit species, identifier namespace, source URL, version and license metadata. No gene list is sent to Enrichr or another service. Source run artifacts are not modified.

```powershell
.venv/Scripts/python.exe -m embryoomics.enrichment --run results/<run> --gmt data/<sets>.gmt --metadata data/<metadata>.yaml --out results/<new_enrichment>
```

The installed wheel also provides `embryoomics-enrich`. Metadata must contain `species`, `id_type` (`gene_id` or `gene_symbol`), `source_url`, `version` and `license`; optional `gmt_sha256` is enforced. No automatic orthology, symbol case conversion, gene-ID version stripping or genome-build guessing is performed.

Background consists of genes detected in raw counts of retained cells. For symbol-based libraries, blank or ambiguous symbols are excluded from both background and query and listed in `gene_mapping.csv`. Exact IDs can instead be used with a matching GMT. Case-preserving tokenization prevents the library's internal uppercase heuristic from conflating IDs. Positive markers with saved adjusted p ≤ the requested `--marker-fdr` (default 0.05) are taken from the existing **truncated top-30 marker table**, not a complete ranked DE result. Gene-set size defaults to 5–500 after intersection with the background.

GSEApy computes local hypergeometric over-representation p-values. The adapter retains **every** eligible cluster–term test, assigns p=1 to zero-hit or empty-query terms, and applies SciPy BH correction over the entire cluster–term family. It does not rely on a correction over only the hit terms returned by the library. Tests compare with a known hypergeometric result and block network requests.

Outputs include all test results, selected markers and mapping exclusions, the full background, effective GMT, input hashes and method metadata, a bilingual explanatory HTML table and an offline ZIP. New output directories are required. This is a separate optional analysis command; it does not silently alter the RNA pipeline or its reviewed report.

## Official mouse GO software example

[MGI's data report page](https://www.informatics.jax.org/downloads/reports/index.html) links the official mouse MOD GAF. [GO terms of use](https://www.geneontology.org/page/use-and-license) describe attribution and version requirements; record the exact downloaded headers and checksums rather than calling the data perpetually current.

```powershell
New-Item -ItemType Directory -Force data/GO_20261113
Invoke-WebRequest https://current.geneontology.org/annotations/gaf/MOUSE-mod.gaf.gz -OutFile data/GO_20261113/MOUSE-mod.gaf.gz
.venv/Scripts/python.exe scripts/prepare_mouse_go.py --gaf data/GO_20261113/MOUSE-mod.gaf.gz --out data/GO_20261113/direct_BP
.venv/Scripts/python.exe -m embryoomics.enrichment --run results/GSM8559287_scanpy_plots_20261113 --gmt data/GO_20261113/direct_BP/mouse_GO_BP_direct.gmt --metadata data/GO_20261113/direct_BP/metadata.yaml --out results/<new_output>
```

The downloaded GAF states generated-by UniProt, date-generated 2026-08-04 and GO version 2026-08-02. The separately downloaded current ontology stated 2026-07-26; the matching dated ontology URL was unavailable during this task. **The ontology was not used.** `prepare_mouse_go.py` therefore defines an explicit **direct BP annotation library**: MGI mouse taxon 10090, aspect P, exclude NOT/ND, include other evidence codes, exclude any symbol mapping to multiple MGI IDs, no ancestor propagation. GO term IDs retain source links. This is not equivalent to a standard ancestor-propagated GO enrichment, and it does not establish absence of untested parent processes.

Real GSM8559287 software example: 21,231 background symbols; 5,248 eligible direct-annotation sets; 47,232 cluster–term tests across 9 unreviewed clusters. 86 tests have global BH ≤ 0.05 under this explicit selection rule. These are exploratory statistical associations, not 86 confirmed pathways or cell identities. Results depend on the background, truncated marker selection, GO release and lack of ancestor propagation. Unknown source genome/annotation version and symbol nomenclature remain limitations.

Local evidence: `results/GSM8559287_GO_direct_BP_final_20261113/`, including `enrichment_manifest.json`, mapping/background files, full results and ZIP. Raw GO data and results remain Git-ignored. The code and download instructions are tracked; no redistributable dataset is implied by a library license.
