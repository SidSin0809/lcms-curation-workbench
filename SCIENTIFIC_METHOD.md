# Scientific and statistical method

## 1. Evidence model

The workbench treats the Progenesis feature identifier as the cross-file key and keeps three conceptual levels separate:

1. **Feature/read level** — one chromatographic feature with its accepted candidate, adduct, isotope, MS², and sample intensities.
2. **Entity level** — all reads sharing an Accepted Compound ID. If the ID is absent, Formula + canonical name is used; feature ID is the final fallback.
3. **Sample level** — biological, pooled-QC, blank, reference, system-suitability, calibration, or excluded injections.

This separation prevents repeated reads from inflating threshold distributions and prevents pooled-QC/blank injections from being treated as biological replicates.

## 2. Representative selection

Within each entity, reads are ordered lexicographically:

1. highest Score;
2. highest Fragmentation Score;
3. smallest absolute Mass Error (ppm);
4. highest Isotope Similarity;
5. highest evidence completeness;
6. greatest linked fragment-peak count;
7. highest pooled-QC detection;
8. lowest pooled-QC CV;
9. highest biological detection;
10. highest biological median abundance;
11. stable Feature ID.

Lexicographic selection is intentionally used instead of a weighted sum: a large change on one vendor-specific scale cannot cancel a poor result on a differently scaled dimension. Score is rounded to its exported precision (0.1), Fragmentation Score to 0.01, absolute mass error to 0.001 ppm, and Isotope Similarity to 0.01 for deterministic tie recognition. The complete rank trace remains attached to every alternative.

## 3. Adaptive identity thresholds

Thresholds are estimated only from the selected entity table.

### Balanced preset

- Score cutoff: entity-level first quartile.
- Isotope cutoff: entity-level first quartile.
- Fragmentation cutoff: first quartile among positive values. Zero is retained as an explicit missing/disabled/unmatched evidence state.
- Absolute mass-error cutoff: a validated distribution boundary when defensible; otherwise a robust upper bound.

### Mass-error separation

A one-Gaussian model is compared with a deterministic two-Gaussian mixture. The two-component boundary is accepted only when:

- ΔBIC = BIC(one component) − BIC(two components) > 10; and
- Ashman D ≥ 2.

The accepted boundary is capped by the user-declared search tolerance. If separation is inadequate, the cutoff is `min(search tolerance, 5 ppm, Q3 absolute ppm)`.

### Stability

Each Balanced cutoff has a deterministic 500-resample nonparametric bootstrap interval. These intervals estimate dataset-specific cutoff stability; they do not estimate identification FDR.

### Presets

- **Inclusive:** discovery-oriented 10th-percentile lower evidence bounds and 90th-percentile mass-error bound; zero fragmentation allowed.
- **Balanced:** lower-quartile evidence with validated/robust mass-error bound; zero fragmentation allowed.
- **Stringent:** distribution separation or upper evidence quantiles; positive fragmentation required.

## 4. Analytical-QC calculations

Calculations use only samples whose roles are included and explicitly mapped.

For positive measured intensities `x`:

- detection rate = number of values > 0 / number of non-missing values;
- CV (%) = `100 × sample SD(x) / mean(x)`;
- robust CV (%) = `100 × 1.4826 × median(|x − median(x)|) / median(x)`;
- D-ratio (%) = `100 × SD(QC) / sqrt(SD(QC)^2 + SD(biological)^2)`;
- biological/blank ratio = biological positive median / blank positive median;
- blank contribution (%) = `100 × blank median / max(blank median, biological median)`.

QC drift is fitted on log2 pooled-QC intensity against mapped injection order. The reported across-run change is `100 × (2^(slope × order span) − 1)`, accompanied by Spearman rho. It is diagnostic by default because nonlinear drift and few QC injections may make a linear summary incomplete.

Technical-replicate CV is calculated within mapped biological units and summarized by the median across eligible units. Dilution response uses Spearman correlation only when at least three distinct mapped dilution factors exist.

Default editable starting points are 80% pooled-QC detection, 30% pooled-QC CV, 20% biological detection, and a 3× biological/blank ratio. The software does not assert that these are universal acceptance limits.

## 5. MS² representation

For each spectrum, intensities are normalized as `p_i = I_i / sum(I)` and spectral entropy is:

`H = −sum(p_i × ln(p_i))`.

Normalized entropy is `H / ln(n)` for `n > 1`, and effective peak count is `exp(H)`. Base peak, peak range, top-five intensity share, low-intensity peak share, declared/observed count integrity, relative intensity, and precursor neutral loss are also exported.

These values describe spectrum complexity and integrity. A reference/decoy spectral-library search is required to translate similarity into a defensible identification error estimate.

## 6. Formula-derived chemistry

For supported elemental formulae, the workbench calculates monoisotopic mass and element counts. Double-bond equivalents are approximated as:

`DBE = 1 + C − (H + F + Cl + Br + I)/2 + N/2`.

H/C, O/C, N/C, heteroatom/C, Kendrick mass and Kendrick mass defect (CH₂ base), chemical-family text cues, lipid-like formula cues, and adduct/ion-mode consistency are exported.

The formula-polarity proxy increases with O, N, P, and S relative to carbon. It is deliberately labeled a proxy: formula alone cannot determine structure, pKa, tautomerism, measured logP, or recovery.

## 7. Source-evidence scoring

Six mutually competing classes are scored:

- human endogenous;
- microbiome-derived;
- host–microbe co-metabolic;
- food-derived;
- drug-derived;
- environmental-derived.

The initial prior changes with biological system, matrix, and assay. Explicit description/ontology rules then add log-evidence weights. Microbial-only, human-only, co-metabolic, food/dietary, drug/pharmaceutical, and pollutant cues receive stronger weights than secondary associations. `associated_not_observed` and expected/predicted-only statuses reduce confidence and remove direct-evidence credit.

Scores are converted to percentage-like values with a temperature-scaled softmax. The export includes normalized entropy, top-two margin, direct-evidence count, rule strength, prior context, evidence list, confidence, and review priority. These values are evidence likelihoods, not calibrated causal probabilities.

## 8. Extraction-location model

The model predicts relative plausibility across organic-rich extract, aqueous/polar extract, SPE retained/eluate, SPE flow-through, interphase/pellet, and uncertain states.

Evidence order:

1. explicit method/behavior labels in the accepted description;
2. parsed MolLogP/TPSA when present;
3. formula-polarity proxy and lipid/class cues;
4. extraction method and analyzed fraction.

Method-specific logic includes upper-organic MTBE, lower-organic Folch/Bligh–Dyer, generic LLE, protein precipitation, reversed-phase/HLB SPE, HILIC/normal-phase SPE, ion-exchange/mixed-mode SPE, and monophasic polar extraction. Ion-exchange results are given higher uncertainty without pKa and loading/elution pH.

The analyzed-phase likelihood is not recovery. Recovery requires matrix-matched standards or spike/recovery experiments because pH, solvent/sample ratio, salt, temperature, emulsion/interphase handling, sorbent loading, wash/elution, and ion suppression are not fully observed in these files.

## 9. Identification claims

The overall Progenesis Score is a mean of enabled/available similarity dimensions; unavailable RT, CCS, or fragmentation can be represented as zero. The supplied data have no accepted-identification RT error and no ACP CCS evidence. Therefore:

- do not treat a universal overall Score as instrument-independent;
- do not treat the adaptive shortlist as target-decoy FDR;
- do not treat exact mass/formula as structural confirmation;
- do not treat source/phase likelihoods as experimentally proven origin/recovery;
- retain repeated reads and manual-review flags.

## 10. Dataset-specific manuscript reporting

The final-stage report is generated after filtering from the actual project settings, sample map, six-file audit, selected representatives, active thresholds, QC-rule states, shortlist, contextual results, and dataset warnings. It reports enabled, disabled, and enabled-but-not-evaluable rules separately. The report also includes formulas, numerical attrition, evidence distributions, contextual class counts, limitations, references, and a checklist of LC, MS, extraction, blank, standard, library-search, and downstream-statistical details that are not recoverable from the six exports.

The generator does not infer missing experimental settings. Its output is a manuscript-ready curation draft that must be verified and merged with laboratory/acquisition records. It cannot create MSI Level-1 evidence, target-decoy FDR, causal source probabilities, or recovery measurements that were not supplied.

## 11. Primary references

- Progenesis QI identification scoring: <https://www.nonlinear.com/progenesis/qi/v3.1/faq/identifications-scoring-algorithm.aspx>
- MSI minimum reporting standards: <https://pubmed.ncbi.nlm.nih.gov/24039616/>
- Lipidomics Minimal Reporting Checklist: <https://www.nature.com/articles/s42255-022-00628-3>
- Pooled-QC guidance: <https://pmc.ncbi.nlm.nih.gov/articles/PMC5960010/>
- mQACC pooled-QC practices: <https://doi.org/10.1021/acs.analchem.3c02924>
- Spectral entropy: <https://www.nature.com/articles/s41592-021-01331-z>
- MiMeDB 2.0: <https://doi.org/10.1093/nar/gkaf1272>
- MTBE lipid extraction: <https://pubmed.ncbi.nlm.nih.gov/18281723/>
- Plasma extraction comparison: <https://pubmed.ncbi.nlm.nih.gov/28000704/>
- Electrospray adduct complexity: <https://pubmed.ncbi.nlm.nih.gov/35058016/>
- Target-decoy context: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6252074/>
