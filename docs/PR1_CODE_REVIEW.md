# PR #1 code review — 2026-11-13

Scope: the complete PR against base `99cc42403962992148a4b399b53a7058d137037f`, initially 61 changed files / 3,465 additions. Reviewed input validation and raw-count handling, QC/confirmation binding, annotation concurrency and rollback, output overwrite/failure handling, sparse exploration, Scanpy/Scrublet/GSEApy reuse, CELLxGENE exchange, metadata pairing, dependency/CI configuration, tests and scientific-boundary documentation.

## Findings resolved before merge

| Finding | Fix and regression evidence |
| --- | --- |
| Displayed annotation CSV and version hash were read separately. A concurrent save between those reads could bind an old table to a new version and allow overwriting newer edits. | `annotation_snapshot` parses and hashes the same byte snapshot. The interface uses both values from that read. A simulated concurrent save confirms the stale draft is rejected. |
| `--overwrite` could remove an existing run that also contained the selected H5AD/config/annotation input. | Reject overwrite when any protected input or configuration resides under the output directory. Regression tests preserve all existing bytes for contained counts and config cases. |
| H5AD numeric counts above int64 maximum could wrap during conversion. | Reject out-of-range raw values before int64 conversion; both unsigned-integer and floating representations are tested. |

Final development and independent locked environments each pass **74 tests**. Ruff basic syntax/undefined-name checks and dependency validation pass. Existing real-data acceptance covers unchanged sparse counts, deterministic baseline clustering, QC/preview agreement, report/ZIP artifacts, optional Scrublet sensitivity, CELLxGENE expression/annotation exchange and offline enrichment. The fixes above do not change the accepted scientific algorithms or baseline parameters.

No unresolved merge-blocking code finding remains from this review. Merge only the tested head after its GitHub CI succeeds; do not bypass repository checks or merge a newer unreviewed head. CI and merge outcomes are recorded in GitHub and Linear.

This approves the software change for development use. Unknown original processing/reference versions, final biological QC/cell-type review and genuine multiomics data remain open. No final scientific release tag is created by merging this PR.
