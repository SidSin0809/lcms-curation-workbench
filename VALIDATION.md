# Validation record

Version 2.1.0 was validated on 2026-08-13 against the complete six-file `Results_Data.zip`, the successful v2.0.1 project-result package, and the supplied GUI screenshots. The regression also retains the earlier v1.1 workbook and SciPy-dependency checks.

## Version 2.1.0 progress, notification, and report regression

The persistent task center receives percentage and named-milestone signals for input scanning, sample-map rebuilding, adaptive filtering, contextual scoring, journal-section generation, CSV writing, XLSX creation, and complete-folder export. Each terminal state writes a visible non-blocking notification and a task-history row containing stage, task, status, local ISO timestamp, duration, and outcome. The history is shown in Step 08 and included as `14_Task_History.csv` plus the `15 Task History` workbook sheet.

The journal generator was tested with the supplied run's active settings: lipidomics, human tissue, MTBE extraction, upper-organic fraction, pooled QC after extraction, Balanced identity cutoffs, and enabled QC/drift/D-ratio rules. It reproduces the final **214 of 1,091** representatives and writes dataset-specific text rather than a static template. Regression assertions cover actual entity counts, shortlist attrition, target-decoy/FDR limitations, the not-recoverable experimental checklist, Markdown export, and the `14 Journal Methods` workbook sheet.

## Version 2.0.1 dependency regression

Version 2.0.0 called `pandas.Series.corr(method="spearman")` for QC injection-order drift and dilution response. Pandas dispatches this path to SciPy, which was intentionally not part of the four-package portable environment. Version 2.0.1 computes Spearman's rho directly as Pearson correlation of deterministic average ranks after pairwise removal of non-finite observations. Constant rank vectors return missing because their correlation is undefined. The result is scientifically equivalent for the coefficient, preserves tie handling, and requires no additional runtime package.

Focused tests cover perfect monotonic and inverse relationships, average-rank ties, pairwise missing values, constant inputs, and an execution guard that rejects every attempted SciPy import.

## Reference full-run signature from the supplied workbook/audit

The supplied v1.1 audit records this complete six-file run:

| Measure | Reference result |
|---|---:|
| CM feature rows | 1,375 |
| AM rows / distinct feature keys | 1,390 / 1,375 |
| CI candidate / accepted rows | 38,850 / 1,375 |
| II compounds / isotope nodes | 1,375 / 142,265 |
| FD spectra / distinct feature keys / peaks | 1,390 / 1,375 / 55,739 |
| ACP rows | 1,375 |
| CM coverage from every source | 100% |
| Unique accepted entities | 1,091 |
| Repeated Accepted-ID groups | 201 |
| Reads in repeated groups | 485 |
| Extra reads collapsed | 284 |

Reference FD fingerprint: size 3.429 MiB; SHA-256 `3fbdcdbe27a98ee0fa5be30ec58154c8eccde5a5904c96356c8bafe32d22ac40`.

## Current attachment-integrity finding

The newest `Results_Data.zip` contains the complete FD source used by the successful v2.0.1 run: 1,390 spectra, 1,375 distinct Feature IDs, 55,739 peaks, 100% CM coverage, zero declared-count mismatches, and SHA-256 `3fbdcdbe27a98ee0fa5be30ec58154c8eccde5a5904c96356c8bafe32d22ac40`. This complete source was used for the v2.1.0 final regression.

An earlier standalone scratch attachment was truncated (821 parsed records and 59.055% CM coverage). The parser's explicit failure/warning behavior for that file remains covered as an input-integrity regression; no missing spectra are reconstructed.

## Entity and threshold regression

The non-FD primary identification metrics reproduce the reference entity signature with the current source set:

| Measure | Result |
|---|---:|
| Accepted feature reads | 1,375 |
| Unique accepted entities | 1,091 |
| Repeated Accepted-ID groups | 201 |
| Repeated-group reads | 485 |
| Extra reads collapsed | 284 |
| Primary-metric ties among selected winners | 0 |
| Fully unresolved manual-review ties | 0 |

Balanced entity-level recommendations remain:

| Metric | Cutoff |
|---|---:|
| Score | ≥ 35.8 |
| Fragmentation Score | ≥ 2.9, or explicit zero state |
| Absolute Mass Error | ≤ 4.8 ppm |
| Isotope Similarity | ≥ 77.9 |

Identity-only Balanced result: **315 of 1,091** representatives. The mass-error boundary remains 4.8318 ppm with ΔBIC 149.37 and Ashman D 2.93.

Version 2 uses 500 deterministic bootstrap resamples. Confidence intervals can differ slightly from the earlier 200-resample export while the recommended point cutoffs remain unchanged.

## Sample and analytical-QC validation

The CM three-row header contains 37 normalized and 37 raw abundance columns:

- 32 samples automatically mapped as Biological (Case or Control);
- 5 pooled-QC injections automatically mapped from the QC condition/sample names;
- 0 mapped blank samples.

Default analytical-QC evaluation over 1,091 representatives produced 941 passes and 150 failures. Applying the successful project's active QC, drift, and D-ratio rules after the Balanced identity filter produced **214 combined passes**. Blank rules were enabled but not evaluable because no blank samples were mapped and therefore did not create observed blank evidence. Because the pooled QC was declared as prepared after extraction, QC dispersion principally assesses injection/instrument and downstream analytical repeatability rather than extraction reproducibility.

## Supplied workbook audit

All 15 supplied workbook sheets were imported, inspected, and rendered. The workbook included run summary, file audit, accepted/repeated/selected/shortlisted/source tables, threshold/method notes, raw CM/AM, accepted CI, ACP, II, and FD summaries. Formula-error scanning found zero `#REF!`, `#DIV/0!`, `#VALUE!`, `#NAME?`, or `#N/A` cells.

The v1.1 workbook's shortlist contains 536 rows because its saved run summary records user-modified settings (Score 36.8 and absolute mass error 10 ppm) rather than the 315-row Balanced preset. Version 2 exports every active manual setting and every entity-level rule decision to remove this ambiguity.

The current CM source also contains 124 compound-link strings ending in the literal vendor placeholder `/NULL`. Version 2 converts these to missing values, emits a run warning, and prevents creation of misleading hyperlinks.

## Automated checks

- Python syntax compilation for all application, launcher, build, analysis, metadata, QC, chemistry, provenance, and export modules.
- Full-source parsing and 1,375/1,091/201/485/284 entity-regression tests.
- Balanced threshold regression tests.
- Sample-role, peak-table integrity, formula-property, and analytical-QC contract tests.
- Contextual lipidomics / human tissue / MTBE / upper-organic scoring.
- Numbered CSV, enhanced XLSX, JSON manifest, and method/limitations export checks.
- Workbook reopening and required-sheet checks with openpyxl.
- Artifact-tool workbook import, formula-error scan, and visual rendering of Read Me, Run Summary, Journal Methods, and Task History. Journal rows use content-aware heights; task columns are widened and wrapped so long text is not clipped.
- Pinned dependency, launcher, release-version, light-theme visibility, and WCAG contrast contracts.
- Final regression result: **16 automated tests passed** against the supplied full-size CM, AM, CI, II, FD, and ACP inputs, covering the SciPy-blocked QC statistics path, dynamic report counts/limitations, progress/notification contracts, task-history outputs, and workbook reopening with the added journal/task sheets.

The Linux execution image used for source validation lacks the host `libEGL.so.1` needed to instantiate Qt widgets, so GUI construction is covered by static compilation and widget/theme contracts here. The runtime itself is pinned to the official PySide6 6.11.1 Windows wheel and is launched in an isolated environment by `RUN_WINDOWS.bat`.
