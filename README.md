# LC–MS Compound Curation Workbench

Portable local Python desktop application for complete six-file Progenesis QI reconciliation, editable sample/QC metadata, peak-level MS² extraction, repeated-identity resolution, analytical-QC review, dataset-adaptive evidence filtering, chemical-property calculation, and contextual source/extraction interpretation.

Version 2.1.0 adds a persistent background-task center, stage-specific progress milestones, non-blocking completion/failure notifications, and an exported task audit. The final panel now generates a dataset-specific, journal-ready Methods and curation-analysis draft from the selected options, mapped samples, active rules, pipeline state, and final shortlist. The Windows visibility fix remains enforced through a platform-independent Fusion style and explicit foreground/background palette contracts. QC Spearman correlations are computed internally from average ranks, so SciPy is not required.

## Start on Windows

1. Extract the ZIP to a new writable folder. Do not run it from inside the ZIP.
2. Install a 64-bit Python 3.11–3.14 if Python is not already available.
3. Double-click **`RUN_WINDOWS.bat`**.
4. The launcher creates `.venv` inside the application folder and installs every pinned runtime dependency with Python/pip. Administrator access is not required.
5. In the first panel, select `CM.csv`, `AM.csv`, `CI.csv`, `II.xml`, `FD.msp`, and `ACP.csv`; define the experiment; then click **Scan all six files**.

First launch requires internet access to install the pinned dependencies. Later launches reuse the local environment. If installation was interrupted, run:

```powershell
py -3 install_and_run.py --repair
```

Runtime dependencies are pinned in `requirements.txt`:

- PySide6 6.11.1 — desktop interface and accessible deterministic palette;
- pandas 3.0.5 — structured joins, sample-level metrics, and CSV handling;
- NumPy 2.3.5 — robust numerical/statistical calculations;
- openpyxl 3.1.5 — formatted multi-sheet XLSX output.

## Eight-stage workflow

1. **Project + files** — captures project, analyst, metabolomics/lipidomics, biological system, matrix, extraction, phase, QC pooling point, ion/acquisition mode, search tolerance, and pooled-QC regex. It hashes and scans every input.
2. **Sample map** — creates an editable row for every normalized CM sample. Roles include Biological, Pooled QC, process/extraction/solvent blank, reference material, system suitability, calibration, excluded, and unknown. Group, subject, technical replicate, batch, injection order, dilution, matrix, extraction, and phase can be imported/exported.
3. **Structure + MS²** — audits row counts, keys, duplicate rows, coverage, file sizes, and SHA-256 hashes. Every MSP spectrum and peak is retained with relative intensity, neutral loss, base peak, entropy, normalized entropy, effective peak count, and integrity flag.
4. **Identity audit** — reconciles Accepted Compound ID, Accepted Description, Formula, Adducts, links, CI counts, isotope information, spectra, chemical properties, and all sample intensities. Repeated accepted IDs are independently searchable/exportable.
5. **Read resolution** — passes unique entities unchanged and selects one representative from every repeated group. The primary hierarchy is Score descending, Fragmentation Score descending, absolute Mass Error ascending, and Isotope Similarity descending. Complete ties use evidence completeness, fragment peaks, pooled-QC detection/CV, biological detection/abundance, and stable Feature ID. Unresolved ties are flagged.
6. **QC + filter** — combines identity thresholds with mapped-sample evidence: pooled-QC detection, standard and robust CV, biological detection, sample/blank ratio, D-ratio, injection-order drift, technical-replicate CV, and optional MS² support. Every representative receives individual pass/fail columns and reasons.
7. **Context model** — calculates formula exact mass, elemental counts, DBE, H/C, O/C, N/C, heteroatom ratios, Kendrick mass defect, lipid/formula cues, adduct polarity, and ion-mode consistency. It then scores human endogenous, microbiome-derived, host–microbe co-metabolic, food-, drug-, and environmental-derived evidence plus method-specific LLE/SPE extraction-location plausibility.
8. **Journal report + export** — generates a dataset-specific Methods and curation-analysis section with actual sample/file/entity/attrition/QC/source counts, equations, rule states, limitations, references, and a completion checklist. It writes numbered stage CSVs, task history, a formatted XLSX workbook, JSON manifest, Markdown/text report, and method records.

## Background progress and task notifications

The persistent task center remains visible below every workflow panel. Full-file scanning, sample-map recomputation, rule application, contextual scoring, manuscript drafting, CSV output, workbook creation, and complete-folder export run in worker threads so the GUI remains responsive. Each task reports named scientific/technical milestones, percentage completion, final outcome, timestamps, and elapsed time. A visible completion or failure notification appears when processing ends; the task history is shown in Step 08 and exported to CSV/XLSX.

Only one mutating/background pipeline task runs at a time, preventing sample maps, thresholds, or outputs from being changed while another stage is using them.

## Dataset-specific journal section

The report generator uses only recorded run evidence. It includes study and extraction context; mapped sample/QC/blank counts; six-file reconciliation; repeated-read selection; MS² processing equations; active adaptive thresholds and bootstrap intervals; QC formulas and evaluability; entity attrition; shortlist evidence distributions; contextual source/phase summaries; limitations; primary references; and a checklist of experimental details that cannot be recovered from these six exports.

Enabled, disabled, and enabled-but-not-evaluable rules are stated separately. Missing blanks, RT/CCS, standards, target-decoy results, or acquisition details are never presented as observed evidence and are never fabricated. The text is a rigorous draft to merge with laboratory and instrument methods—not an automatic claim of MSI Level 1, identification FDR, causal origin, or measured extraction recovery.

Project settings and the reviewed sample map can be saved to or restored from a JSON project file. The project file stores paths and settings, not copies of the LC–MS data.

## Threshold strategy

Recommendations are calculated after one representative has been selected per identity, preventing identities with multiple chromatographic reads from overweighting the distributions.

- **Score** and **Isotope Similarity** use entity-level lower quartiles for the Balanced preset.
- **Fragmentation Score** uses the lower quartile among positive values. Zero remains an explicit unavailable/disabled/unmatched evidence state in Inclusive and Balanced modes; Stringent requires positive support.
- **Absolute Mass Error** uses a deterministic two-cluster boundary only when a two-Gaussian model has ΔBIC > 10 and Ashman D ≥ 2. Otherwise it uses a robust upper quartile capped at 5 ppm and at the declared search tolerance.
- Each recommended cutoff has a deterministic 500-resample bootstrap interval.
- Inclusive, Balanced, and Stringent presets remain editable. These are triage rules, not target-decoy FDR or compound confirmation.

## Pooled QC and blank logic

The tool never treats absent QC/blank evidence as observed quality.

- QC detection, mean/median, conventional CV, robust MAD-based CV, D-ratio, and injection-order drift require mapped pooled-QC injections.
- Pooled QC prepared before extraction includes preparation plus analytical variation; pooling after extraction mainly measures analytical variation. The declared pooling point is preserved.
- Biological detection and abundance use only included Biological samples.
- Blank detection, blank abundance, biological/blank median ratio, and blank contribution require mapped blank samples.
- Batch, injection order, dilution, and technical-replicate identifiers enable additional diagnostics when sufficient data exist.
- Default QC values are editable starting points (80% QC detection, 30% QC CV, 20% biological detection, 3× sample/blank ratio). They are not universal acceptance limits.

## Source and extraction interpretation

The source model combines explicit text/ontology evidence with declared biological-system, sample-matrix, assay, and exposure priors. It exports class probabilities, top-class margin, normalized entropy, direct-evidence count, evidence strength, provenance, confidence, and review priority. Prior-only or association-not-observed results remain low confidence.

The phase model distinguishes:

- MTBE and other upper-organic LLE;
- Folch/Bligh–Dyer lower-organic LLE;
- aqueous/polar fractions;
- protein-precipitation supernatant versus pellet/interphase;
- reversed-phase/HLB, HILIC/normal-phase, and ion-exchange/mixed-mode SPE retained/eluate versus flow-through;
- unknown/uncertain conditions.

Formula polarity is a transparent proxy. When a description contains MolLogP or TPSA, these values are parsed and preserved. None of these outputs is a measured partition coefficient or recovery probability; pH, pKa, solvent ratio, salt, matrix, temperature, emulsions, sorbent chemistry, wash/elution conditions, and ionization remain experimentally important.

## Complete output package

**Export complete results folder** writes:

- `00_Project_Metadata.csv`
- `01_Sample_Metadata.csv`
- `02_File_Audit.csv`
- `03_All_Accepted_Feature_Reads.csv`
- `04_Repeated_Accepted_ID_Reads.csv`
- `05_MS2_Spectra.csv`
- `06_MS2_Fragment_Peaks.csv`
- `07_Selected_Compound_Representatives.csv`
- `08_Analytical_QC.csv`
- `09_Filter_Decisions_All_Entities.csv`
- `10_Filtered_Compound_Shortlist.csv`
- `11_Chemistry_Source_and_Phase.csv`
- `12_Threshold_Method.csv`
- `13_Method_Notes.csv`
- `14_Task_History.csv`
- `JOURNAL_READY_METHODS_AND_ANALYSIS.md`
- `JOURNAL_READY_METHODS_AND_ANALYSIS.txt`
- `LCMS_Compound_Curation_Results.xlsx`
- `analysis_manifest.json`
- `METHOD_AND_LIMITATIONS.txt`
- `SCIENTIFIC_METHOD.md`

The XLSX workbook also includes dataset-specific Journal Methods and Task History sheets, a Read Me sheet, raw CM and AM, accepted CI rows, ACP, II summary, and optionally every CI candidate. CSV/XLSX text beginning with spreadsheet formula characters is neutralized to reduce formula-injection risk.

## Scientific limitations

- Progenesis Score depends on enabled evidence dimensions. RT/CCS/fragmentation values of zero can mean unused evidence, so a universal overall cutoff is not defensible.
- Exact formula or accurate mass does not resolve structural isomers.
- Repeated Accepted IDs can be isomers, in-source products, alternative adduct behavior, or false matches. Keep the repeated-read audit.
- Spectrum entropy and peak counts describe complexity; without a reference/decoy search they do not establish identity or FDR.
- MSI Level-1 identification requires appropriate orthogonal evidence and an authentic standard analyzed under comparable conditions.
- Source scores are evidence-weighted heuristics, not calibrated causal or clinical probabilities.
- Phase scores estimate plausibility, not extraction recovery.

Primary references and URLs are embedded in every complete export, including Progenesis QI scoring, MSI reporting standards, mQACC pooled-QC guidance, spectral entropy, MiMeDB, MTBE extraction, extraction-method comparison, adduct complexity, and target-decoy context.

## Build a standalone Windows folder/executable

Double-click `BUILD_WINDOWS_PORTABLE.bat` on Windows. It creates a dedicated build environment and produces:

- `dist/LCMS_Compound_Curation/LCMS_Compound_Curation.exe` with runtime files;
- `dist/LCMS_Compound_Curation_Windows_Portable.zip`.

Build the Windows executable on Windows; PyInstaller bundles separately for each operating system.

## Validation

From the application folder:

```powershell
set LCMS_TEST_DATA=C:\path\to\six_files
.venv\Scripts\python -m unittest discover -s tests -v
```

See `VALIDATION.md` and `DATA_AUDIT.md` for the supplied-data signature, workbook audit, current attachment integrity finding, and automated checks.
