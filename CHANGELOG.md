# Changelog

## 2.1.0 — 2026-08-13

- Added a persistent background-pipeline progress bar visible at every workflow stage.
- Moved sample-map rebuilding, filtering, contextual scoring, journal drafting, CSV writing, XLSX creation, and complete-folder export into responsive worker tasks with named progress milestones.
- Added non-blocking success/failure notifications plus a timestamped, duration-aware task-history audit in the GUI, CSV package, and XLSX workbook.
- Added a dataset-specific high-impact-journal Methods and curation-analysis generator driven by declared assay/matrix/extraction settings, mapped samples, file audits, active thresholds/QC rules, shortlist evidence, contextual scores, and dataset warnings.
- Added actual attrition and evidence summaries, equations, enabled/evaluable/disabled rule states, primary references, a reproducibility statement, and a do-not-fabricate manuscript completion checklist.
- Added `JOURNAL_READY_METHODS_AND_ANALYSIS.md`, its text companion, `14_Task_History.csv`, and dedicated Journal Methods and Task History workbook sheets.
- Extended the machine-readable manifest with curation outcomes and critical interpretation limitations.
- Added regression coverage for the 214-compound supplied-run shortlist, dynamic report contents, background UI contracts, and new workbook sheets.

## 2.0.1 — 2026-08-13

- Fixed the Windows scan failure caused by pandas dynamically importing the unlisted SciPy package for Spearman correlation.
- Replaced both QC injection-order and dilution-response correlations with a NumPy implementation of Spearman's rho using pairwise-finite observations and deterministic average ranks for ties.
- Added regression coverage for monotonic, inverse, tied, missing, constant, and explicitly SciPy-blocked execution paths.
- Kept the portable environment lean: SciPy is not required and the four pinned runtime dependencies remain complete.

## 2.0.0 — 2026-08-13

- Reorganized the GUI into eight auditable stages: project/files, sample map, structure/MS², identity audit, representative resolution, QC/filtering, contextual interpretation, and export.
- Added editable/importable/exportable sample metadata with pooled QC, multiple blank types, biological samples, reference/system-suitability/calibration roles, exclusions, groups, subjects, technical replicates, batches, injection order, dilution, matrix, extraction, and phase.
- Added pooled-QC detection, conventional CV, robust MAD-based CV, D-ratio, QC drift, dilution response, technical-replicate CV, biological detection, and blank-contribution metrics.
- Added peak-level MSP export with relative intensity, neutral loss, entropy, normalized entropy, effective peak count, low-intensity/top-five shares, and integrity flags.
- Added formula exact mass, element counts, DBE, elemental ratios, formula-polarity proxy, lipid-like flag, Kendrick mass defect, chemical-family cues, adduct polarity, and ion-mode consistency.
- Added method-specific MTBE/Folch/Bligh–Dyer/generic LLE, protein-precipitation, HLB/C18/HILIC/normal-phase/ion-exchange SPE, and monophasic extraction-location models.
- Expanded source scoring for human, microbiome, host–microbe, food, drug, and environmental evidence across human biofluid/tissue, microbiome, cell-line, bacterial, viral-exposure, and environmental contexts; added entropy, margin, direct-evidence count, and provenance.
- Added rule-level filter decisions for every representative, 500-resample threshold intervals, project save/load JSON, complete numbered CSVs, enhanced XLSX, manifest, and method/limitations text.
- Added explicit detection of failed source coverage and incomplete/truncated MSP terminal records.
- Added sanitization and dataset warnings for vendor compound-link placeholders ending in `NULL`/`NA`.

## 1.1.0 — 2026-08-13

- Fixed invisible Windows text caused by a dark native palette inheriting into white stylesheet panels.
- Forced Qt's cross-platform Fusion style and applied explicit Active, Inactive, Disabled, selection, placeholder, tooltip, and status-bar palette roles.
- Added complete foreground/background styling for labels, group boxes, inputs, dropdown lists, tables, progress bars, checkboxes, and disabled controls.
- Split the header `Method & version` and page `Clear paths` button styles so both remain readable.
- Added a Qt-free theme regression and contrast test suite.
- Added Accepted Description to the identity, repeated-read winner, shortlist, and source/phase GUI previews.
- Improved the Windows launcher error message, UTF-8 handling, and dependency-install diagnostics.

## 1.0.0 — 2026-08-13

- Initial six-file LC–MS curation, ranking, adaptive-threshold, source/phase scoring, CSV, and XLSX release.
