# LC–MS Curation Workbench

**Auditable compound curation, analytical QC, MS² evidence processing, and reproducible reporting for untargeted LC–MS workflows**

[![Version](https://img.shields.io/badge/version-2.1.0-blue.svg)](https://github.com/SidSin0809/lcms-curation-workbench/releases)
[![Python](https://img.shields.io/badge/Python-3.11--3.14-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)](#installation)
[![Research Software](https://img.shields.io/badge/type-research%20software-purple.svg)](#citation)

**LC–MS Curation Workbench** is a local Python desktop application for systematic reconciliation, quality control, curation, interpretation, and export of compound-identification results from **Progenesis QI LC–MS/MS workflows**.

The software converts six complementary Progenesis QI exports into an auditable entity-level workflow integrating:

- cross-file feature and identity reconciliation;
- editable biological/QC/blank sample metadata;
- repeated Accepted Compound ID resolution;
- peak-level MS² extraction and spectrum-quality descriptors;
- analytical-QC assessment;
- dataset-adaptive identity filtering;
- chemical-property and ion-mode consistency checks;
- contextual source and extraction-phase interpretation;
- provenance-aware filtering decisions; and
- dataset-specific Methods and curation-analysis reporting for manuscript preparation.

> **Current release: v2.1.0**

---

## Why this tool?

Untargeted LC–MS workflows often distribute relevant identification evidence across multiple export files while leaving several scientifically important decisions-duplicate identity resolution, QC assessment, threshold selection, missing-evidence handling, and reporting-to manual post-processing.

LC–MS Curation Workbench provides a single reproducible workflow in which each compound can be traced from the original feature-level evidence to the final shortlist.

The emphasis is not simply on filtering compounds, but on making **why a compound was retained, rejected, collapsed, flagged, or interpreted explicitly auditable**.

---

## Core capabilities

### 1. Six-file reconciliation

The workflow jointly analyzes:

- `CM.csv`
- `AM.csv`
- `CI.csv`
- `II.xml`
- `FD.msp`
- `ACP.csv`

1. Compound Measurements (CM.csv)
2. Adduct Measurements (AM.csv)
3. Compound Identifications (CI.csv)
4. Isotope Information (II.xml)
5. Fragment Database (FD.msp)
6. Additional Compound Properties (ACP.csv)

Files are structurally scanned, cross-reconciled through Progenesis feature identifiers, and recorded with file-level audit information and SHA-256 fingerprints.

Missing or incomplete source coverage is explicitly reported rather than silently inferred.

### 2. Editable sample and QC metadata

Every normalized CM abundance column is represented in an editable sample map.

Supported roles include:

- Biological sample
- Pooled QC
- Process/extraction blank
- Solvent blank
- Reference material
- System-suitability sample
- Calibration sample
- Excluded sample
- Unknown

Additional metadata include biological group, subject, technical replicate, analytical batch, injection order, dilution, sample matrix, extraction method, and analyzed phase.

### 3. MS² evidence extraction

MSP spectra are retained at both spectrum and fragment-peak levels.

Derived descriptors include:

- fragment m/z and intensity;
- relative intensity;
- precursor neutral loss;
- base peak;
- spectral entropy;
- normalized spectral entropy;
- effective peak count;
- low-intensity peak fraction;
- top-five intensity share; and
- declared-versus-observed peak-count integrity.

These are descriptive evidence metrics and are not presented as identification FDR.

### 4. Repeated-identity resolution

Accepted Compound IDs occurring in multiple chromatographic features are retained and audited before entity-level collapse.

Representative selection follows a deterministic lexicographic hierarchy:

1. higher Progenesis Score;
2. higher Fragmentation Score;
3. lower absolute mass error;
4. higher Isotope Similarity.

Residual ties use evidence completeness, MS² information, pooled-QC behavior, biological detection/abundance, and stable Feature ID.

Unresolved cases remain explicitly flagged for manual review.

### 5. Analytical QC

Where appropriate mapped samples are available, the workbench calculates:

- pooled-QC detection rate;
- conventional CV;
- robust MAD-based CV;
- D-ratio;
- biological detection rate;
- blank contribution;
- biological/blank abundance ratio;
- injection-order drift;
- technical-replicate CV; and
- dilution-response correlation.

Rules that cannot be evaluated because relevant QC, blank, dilution, or replicate samples are absent are classified as **not evaluable**, rather than automatically passed.

### 6. Dataset-adaptive identity filtering

One representative per accepted identity is used to estimate dataset-specific evidence distributions.

The **Balanced** strategy uses lower-quartile evidence thresholds with a validated or robust absolute-mass-error boundary.

Mass-error distribution separation is accepted only when a two-component model provides:

- ΔBIC > 10; and
- Ashman D ≥ 2.

Recommended thresholds additionally receive deterministic bootstrap intervals.

Inclusive, Balanced, and Stringent presets remain editable.

> These thresholds are curation/triage rules. They are not target-decoy FDR estimates and do not constitute structural confirmation.

### 7. Formula and chemical-context analysis

For supported molecular formulae, the software calculates:

- monoisotopic exact mass;
- elemental composition;
- double-bond equivalents;
- H/C, O/C, and N/C ratios;
- heteroatom/carbon ratio;
- Kendrick mass and Kendrick mass defect;
- formula-polarity proxy;
- lipid-like formula cues;
- chemical-family text cues;
- adduct polarity; and
- ion-mode consistency.

Formula-derived descriptors are explicitly treated as proxies and not as measurements of structure, pKa, logP, recovery, or biological origin.

### 8. Contextual source and extraction interpretation

Evidence can be summarized across six contextual classes:

- human endogenous;
- microbiome-derived;
- host–microbe co-metabolic;
- food-derived;
- drug-derived;
- environmental-derived.

The model exports the underlying evidence, top-class margin, normalized entropy, evidence strength, confidence, provenance, and manual-review priority.

Extraction-location models support contexts including MTBE, Folch/Bligh–Dyer, generic liquid–liquid extraction, protein precipitation, reversed-phase/HLB SPE, HILIC/normal-phase SPE, ion-exchange/mixed-mode SPE, and monophasic extraction.

These outputs represent **evidence-weighted plausibility**, not experimentally measured recovery or causal origin.

---

## Eight-stage workflow

| Stage | Function |
|---|---|
| **01 — Project + files** | Experimental metadata, six input files, hashing and structural scan |
| **02 — Sample map** | Biological/QC/blank roles and sample metadata |
| **03 — Structure + MS²** | Cross-file audit and peak-level MS² processing |
| **04 — Identity audit** | Accepted IDs, formulae, adducts and evidence reconciliation |
| **05 — Read resolution** | Deterministic representative selection |
| **06 — QC + filter** | Identity evidence and analytical-QC filtering |
| **07 — Context model** | Chemistry, source evidence and extraction plausibility |
| **08 — Journal report + export** | Reproducible outputs and manuscript-oriented reporting |

---

## Background processing and progress tracking

Version 2.1.0 introduces a persistent task center for computationally intensive operations.

Background workers are used for:

- full-file scanning;
- sample-map recomputation;
- filtering;
- contextual scoring;
- manuscript drafting;
- CSV generation;
- workbook generation; and
- complete-folder export.

Each task records:

- scientific/technical stage;
- progress milestone;
- completion percentage;
- outcome;
- start/end timestamps; and
- elapsed time.

Completion and failure notifications are non-blocking, and task history is exported with the analysis.

---

## Journal-oriented Methods and analysis reporting

The final workflow stage generates a dataset-specific manuscript draft from the **actual recorded analysis state**, including:

- project and extraction context;
- mapped sample/QC/blank counts;
- six-file reconciliation;
- identity attrition;
- repeated-read resolution;
- MS² calculations;
- selected thresholds;
- analytical-QC rules;
- enabled/disabled/not-evaluable criteria;
- final shortlist statistics;
- source/context distributions;
- equations;
- limitations;
- references; and
- a manuscript-completion checklist.

The generator does not invent unavailable experimental parameters.

The resulting document is intended to be verified and integrated with the laboratory, chromatography, acquisition, and downstream statistical methods of the study.

---

## Installation

### Windows — recommended

1. Download the latest portable source package or release.
2. Extract it to a writable directory.
3. Do **not** execute it from inside the ZIP.
4. Double-click:

```text
RUN_WINDOWS.bat
