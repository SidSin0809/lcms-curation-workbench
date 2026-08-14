# Supplied-data and interface audit

## Sources reviewed

- Original six filenames: CM.csv, AM.csv, CI.csv, II.xml, FD.msp, ACP.csv.
- Nine stage/export files plus `analysis_manifest.json` from the earlier run.
- Eight GUI screenshots covering input, identity, representative, threshold, and source/phase panels.
- The complete multi-sheet XLSX workbook, inspected sheet-by-sheet and formula-error scanned.

## Structural understanding

- CM and AM use three header rows. Sample abundance fields encode section, condition/group, and sample leaf name. CM supplies 37 normalized and 37 raw sample columns.
- CI is a 38,850-row candidate table; 1,375 rows are user-accepted. Candidate multiplicity is retained separately from accepted feature reads.
- II is hierarchical XML containing compound, condition, sample, adduct, and isotope nodes; the reference file contains 142,265 isotope nodes.
- FD is MSP with spectrum metadata followed by m/z–intensity peak pairs. Version 2 exports both spectrum-level and peak-level tables.
- ACP supplies accepted-compound property rows. In the supplied data, CCS and ACP adduct are absent for every row.
- Cross-file joins use the Progenesis feature identifier. Accepted Compound ID is used only after feature-level evidence is reconciled.

## Data-quality findings addressed in v2

1. Sample roles were previously inferred only by a regex and were not editable/exported as a complete metadata map.
2. Pooled-QC CV existed, but robust CV, D-ratio, drift, technical-replicate, dilution, biological detection, and blank metrics were absent.
3. MSP information was summarized but individual fragment peaks were not exported.
4. Wide tables exposed all columns at once, making long headers difficult to navigate. V2 uses focused panel views, full tooltips, movable headers, capped initial widths, and complete exports.
5. Threshold settings in the supplied workbook differ from the Balanced preset; entity-level pass/fail reasons were not exported. V2 adds a complete filter-decision table.
6. Source scores were predominantly low-confidence and many rows repeated prior-driven probabilities. V2 exports entropy, class margin, direct-evidence count, evidence strength, prior context, and review priority.
7. Extraction scoring lacked explicit SPE retained/flow-through logic and comprehensive formula-derived chemistry. V2 adds method-specific LLE/SPE logic and transparent chemical-property proxies.
8. An earlier standalone scratch copy of FD.msp was truncated; V2 correctly flagged failed coverage and incomplete terminal records. The newest `Results_Data.zip` contains the complete 1,390-spectrum/55,739-peak FD source and was used for the final v2.1 regression.
9. CM contains 124 vendor-generated compound URLs whose terminal identifier is the literal placeholder `NULL`. V2 blanks these values, records a dataset warning, and does not export them as usable hyperlinks.
10. The successful attached project used lipidomics/human-tissue/MTBE/upper-organic context, pooled QC after extraction, Balanced identity evidence, and enabled QC/drift/D-ratio rules. This produced 214 final compounds from 1,091 representatives. V2.1 generates the manuscript draft from those actual values and states that after-extraction QC does not assess extraction reproducibility.

## Interpretation safeguards

- Missing blank, RT, CCS, standard, or target-decoy evidence is explicitly unavailable rather than treated as a pass.
- A pooled QC is not treated as a biological replicate.
- Representative selection is lexicographic so differently scaled metrics cannot cancel each other.
- Formula descriptors are labeled as formula-derived proxies.
- Spectral entropy is descriptive unless a reference/decoy comparison is supplied.
- Source and extraction outputs are evidence likelihoods/plausibilities, not calibrated causal probabilities or measured recoveries.
