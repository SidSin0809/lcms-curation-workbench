from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd

from .engine import AnalysisBundle, evaluate_filters


ProgressCallback = Callable[[int, str], None]


REFERENCES = (
    (
        "1",
        "Sumner LW et al. Proposed minimum reporting standards for chemical analysis: Chemical Analysis Working Group, Metabolomics Standards Initiative. Metabolomics. 2007;3:211–221.",
        "https://doi.org/10.1007/s11306-007-0082-2",
    ),
    (
        "2",
        "McDonald JG et al. Introducing the Lipidomics Minimal Reporting Checklist. Nature Metabolism. 2022;4:1086–1088.",
        "https://doi.org/10.1038/s42255-022-00628-3",
    ),
    (
        "3",
        "mQACC. Current practices in LC–MS untargeted metabolomics: a scoping review on the use of pooled quality-control samples. Analytical Chemistry. 2023;95:18645–18654.",
        "https://doi.org/10.1021/acs.analchem.3c02924",
    ),
    (
        "4",
        "Broadhurst D et al. Guidelines and considerations for system-suitability and quality-control samples in untargeted clinical metabolomics. Metabolomics. 2018;14:72.",
        "https://doi.org/10.1007/s11306-018-1367-3",
    ),
    (
        "5",
        "Progenesis QI identification-scoring algorithm (vendor technical documentation).",
        "https://www.nonlinear.com/progenesis/qi/v3.1/faq/identifications-scoring-algorithm.aspx",
    ),
    (
        "6",
        "Li Y et al. Spectral entropy outperforms MS/MS dot product similarity for small-molecule identification. Nature Methods. 2021.",
        "https://doi.org/10.1038/s41592-021-01331-z",
    ),
    (
        "7",
        "Scheubert K et al. Significance estimation for large-scale metabolomics annotations by spectral matching. Nature Communications. 2017;8:1494.",
        "https://doi.org/10.1038/s41467-017-01318-5",
    ),
    (
        "8",
        "Matyash V et al. Lipid extraction by methyl-tert-butyl ether for high-throughput lipidomics. Journal of Lipid Research. 2008;49:1137–1146.",
        "https://doi.org/10.1194/jlr.D700041-JLR200",
    ),
    (
        "9",
        "MiMeDB 2.0: the Human Microbial Metabolome Database.",
        "https://doi.org/10.1093/nar/gkaf1272",
    ),
)


def _emit(progress: ProgressCallback | None, percent: int, message: str) -> None:
    if progress is not None:
        progress(max(0, min(100, int(percent))), message)


def _text(value: Any, fallback: str = "not declared") -> str:
    if value is None:
        return fallback
    try:
        if pd.isna(value):
            return fallback
    except (TypeError, ValueError):
        pass
    cleaned = str(value).strip()
    return cleaned if cleaned else fallback


def _humanize(value: Any) -> str:
    text = _text(value)
    labels = {
        "mtbe": "MTBE biphasic liquid–liquid extraction",
        "folch": "Folch chloroform/methanol liquid–liquid extraction",
        "bligh-dyer": "Bligh–Dyer liquid–liquid extraction",
        "dia-mse": "DIA/MSE",
        "hdmse": "HDMSE",
        "dda": "data-dependent MS/MS",
        "aif": "all-ion fragmentation",
        "upper-organic": "upper-organic",
        "lower-organic": "lower-organic",
        "upper-aqueous": "upper-aqueous",
        "lower-aqueous": "lower-aqueous",
    }
    return labels.get(text.casefold(), text.replace("-", " ").replace("_", " "))


def _pct(numerator: int | float, denominator: int | float) -> str:
    return f"{(float(numerator) / float(denominator) * 100):.1f}%" if denominator else "not estimable"


def _quantiles(frame: pd.DataFrame, column: str, absolute: bool = False) -> str:
    if column not in frame or frame.empty:
        return "not estimable"
    values = pd.to_numeric(frame[column], errors="coerce")
    if absolute:
        values = values.abs()
    values = values[np.isfinite(values)]
    if values.empty:
        return "not estimable"
    q1, median, q3 = values.quantile([0.25, 0.5, 0.75]).tolist()
    return f"median {median:.2f} (IQR {q1:.2f}–{q3:.2f})"


def _numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _preset_name(bundle: AnalysisBundle, settings: dict[str, Any]) -> str:
    keys = ("score", "fragmentation", "abs_mass_error", "isotope", "allow_zero_fragmentation")
    for preset in bundle.thresholds.get("presets", []):
        if all(
            bool(settings.get(key)) == bool(preset.get(key))
            if key == "allow_zero_fragmentation"
            else abs(float(settings.get(key, np.nan)) - float(preset.get(key, np.nan))) < 1e-9
            for key in keys
        ):
            return str(preset.get("name", "preset"))
    return "manual/custom"


def _rule_state(enabled: bool, evaluable: bool) -> str:
    if not enabled:
        return "Disabled"
    return "Enabled and evaluable" if evaluable else "Enabled but not evaluable"


def _failure_counts(decisions: pd.DataFrame, column: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    if column not in decisions:
        return counts
    for value in decisions[column].dropna().astype(str):
        if value in {"", "No identification-evidence rule failed", "No enabled analytical-QC rule failed"}:
            continue
        for reason in value.split(";"):
            reason = reason.strip()
            if reason:
                counts[reason] += 1
    return counts


def _markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    def safe(value: Any) -> str:
        return _text(value, "—").replace("|", "\\|").replace("\n", " ")

    output = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    output.extend("| " + " | ".join(safe(value) for value in row) + " |" for row in rows)
    return "\n".join(output)


def build_journal_sections(
    bundle: AnalysisBundle,
    shortlist: pd.DataFrame,
    provenance: pd.DataFrame,
    filter_settings: dict[str, Any],
    provenance_settings: dict[str, Any] | None,
    progress: ProgressCallback | None = None,
) -> list[dict[str, str]]:
    """Create dataset-specific, manuscript-ready methods and curation analysis."""

    _emit(progress, 3, "Reading project, sample, and file-audit metadata…")
    provenance_settings = provenance_settings or {}
    parameters = {**bundle.parameters, **provenance_settings}
    decisions = evaluate_filters(bundle.selected_reads, filter_settings) if filter_settings else pd.DataFrame()
    included = bundle.sample_metadata.loc[bundle.sample_metadata["Include"].astype(bool)].copy()
    roles = included["Sample role"].value_counts()
    biological_n = int(roles.get("Biological", 0))
    qc_n = int(roles.get("Pooled QC", 0))
    blank_n = int(sum(roles.get(role, 0) for role in ("Process blank", "Extraction blank", "Solvent blank")))
    biological_units = int(
        included.loc[included["Sample role"].eq("Biological"), "Subject / biological unit"].replace("", np.nan).nunique()
    )
    technical_replicates = int(
        included.loc[included["Sample role"].eq("Biological"), "Technical replicate"].astype(bool).sum()
    )
    group_counts = (
        included.loc[included["Sample role"].eq("Biological"), "Condition / group"]
        .replace("", "Unspecified")
        .value_counts()
        .to_dict()
    )
    group_text = ", ".join(f"{key}: {value}" for key, value in group_counts.items()) or "not declared"
    pooling = _text(parameters.get("qc_pooling_point"), "unknown").casefold()
    pooling_interpretation = (
        "Because the pooled QC was declared as pooled after extraction, its variance principally monitors injection/instrument and downstream analytical repeatability; it does not estimate extraction reproducibility."
        if "after" in pooling
        else "Because the pooled QC was declared as pooled before extraction, its variance includes sample-preparation and analytical components."
        if "before" in pooling
        else "The QC pooling point was not declared; the relative contribution of extraction and analytical variance therefore cannot be separated."
    )

    audit_rows = []
    for _, row in bundle.audits.iterrows():
        audit_rows.append(
            [
                row.get("Role"),
                f"{int(row.get('Records', 0)):,}",
                f"{int(row.get('Distinct feature keys', 0)):,}",
                f"{float(row.get('CM coverage (%)', np.nan)):.3f}%" if pd.notna(row.get("CM coverage (%)")) else "not estimable",
                row.get("Status"),
            ]
        )
    audit_table = _markdown_table(["Source", "Records", "Distinct feature keys", "CM coverage", "Audit"], audit_rows)

    sections: list[dict[str, str]] = []
    sections.append(
        {
            "section": "Methods",
            "subsection": "Study context and scope",
            "text": (
                f"An untargeted {_humanize(parameters.get('assay_type'))} dataset from {_humanize(parameters.get('sample_matrix'))} "
                f"({ _humanize(parameters.get('biological_system')) } biological context) was curated in {_humanize(parameters.get('ion_mode'))} ion mode. "
                f"The declared acquisition mode was {_humanize(parameters.get('acquisition_mode'))}, the precursor-mass search tolerance was "
                f"{float(parameters.get('mass_tolerance', np.nan)):.1f} ppm, and the declared extraction context was {_humanize(parameters.get('extraction_method'))} "
                f"with analysis of the {_humanize(parameters.get('analyzed_phase'))} fraction. The workflow curated the six exported Progenesis QI evidence files; "
                "it did not reconstruct chromatographic integration or rerun database searching. Reporting was structured to retain the metadata, sample-preparation, "
                "quality-control, identification, and preprocessing information emphasized by the Metabolomics Standards Initiative and, for lipidomics studies, the "
                "Lipidomics Minimal Reporting Checklist [1,2]."
            ),
        }
    )
    sections.append(
        {
            "section": "Methods",
            "subsection": "Sample metadata and quality-control design",
            "text": (
                f"The three-row CM abundance header yielded {len(bundle.normalized_columns):,} normalized and {len(bundle.raw_columns):,} raw abundance columns. "
                f"After explicit sample-role review, {biological_n:,} biological injections representing {biological_units:,} biological units "
                f"({technical_replicates:,} injections marked as technical replicates), {qc_n:,} pooled-QC injections, and {blank_n:,} blank injections were included. "
                f"Biological-group membership was {group_text}. Pooled QC, biological samples, blanks, reference materials, calibration samples, system-suitability injections, "
                "and excluded injections were never treated as interchangeable. "
                f"{pooling_interpretation} This explicit role mapping follows recommendations that pooled QCs be reported in sufficient detail and used to quantify repeatability, "
                "drift, and feature quality [3,4]."
            ),
        }
    )
    sections.append(
        {
            "section": "Methods",
            "subsection": "Input reconciliation and provenance",
            "text": (
                "CM (compound measurements), AM (adduct/ion measurements), CI (all identification candidates), II (isotope/ion XML), FD (MSP fragmentation spectra), "
                "and ACP (accepted-compound properties) were scanned in full. SHA-256 digests, record/field counts, duplicate keys, declared MS² peak-count integrity, and "
                "coverage relative to CM were retained. Cross-file joins used the Progenesis feature identifier; no row was joined on compound name alone. The source audit was:\n\n"
                + audit_table
            ),
        }
    )

    _emit(progress, 24, "Documenting identity reconciliation and MS² processing…")
    sections.append(
        {
            "section": "Methods",
            "subsection": "Accepted identities and repeated-read resolution",
            "text": (
                f"The CI table contained {len(bundle.ci):,} candidate assignments, of which {len(bundle.accepted_reads):,} feature-level reads carried a user-accepted identity. "
                "Entities were keyed by Accepted Compound ID; missing IDs fell back to formula plus canonicalized accepted description and finally to Feature ID. "
                f"This produced {len(bundle.selected_reads):,} accepted entities. {bundle.duplicate_groups:,} Accepted Compound IDs mapped to more than one chromatographic feature "
                f"({len(bundle.duplicate_reads):,} reads), and {bundle.extra_reads_collapsed:,} nonrepresentative reads were collapsed only in the representative matrix while being preserved in a separate audit table. "
                "For each repeated entity, reads were ranked lexicographically by decreasing overall Score, decreasing Fragmentation Score, increasing absolute mass error, and decreasing Isotope Similarity. "
                "Exact primary-metric ties were resolved by decreasing evidence completeness, MS² peak support, pooled-QC detection, biological detection and abundance, then by increasing QC CV and a stable Feature ID. "
                "Lexicographic ordering was used instead of an additive composite so that vendor metrics on different scales could not compensate for one another. The selected read, rank, reason, full comparison trace, and unresolved-manual-review flag were retained."
            ),
        }
    )
    ms2_supported = int(_numeric_series(bundle.selected_reads, "FD observed peaks").fillna(0).gt(0).sum())
    sections.append(
        {
            "section": "Methods",
            "subsection": "Tandem-MS representation",
            "text": (
                f"All {len(bundle.fd):,} MSP spectra and {len(bundle.fd_peaks):,} fragment peaks were exported at spectrum and peak levels. "
                f"Among entity representatives, {ms2_supported:,}/{len(bundle.selected_reads):,} had at least one linked fragment peak. For spectrum intensities Iᵢ, pᵢ = Iᵢ/ΣI, "
                "spectral entropy was calculated as H = −Σpᵢln(pᵢ), normalized entropy as H/ln(n) for n > 1, and effective peak count as exp(H). Base peak, top-five intensity share, "
                "low-intensity share, relative intensity, neutral loss, and declared-versus-observed peak-count integrity were also retained [6]. These descriptors characterize spectrum complexity and integrity; "
                "without a reference/decoy spectral-library search they do not estimate identification significance or false-discovery rate [7]."
            ),
        }
    )

    _emit(progress, 42, "Summarizing adaptive thresholds and analytical-QC rules…")
    preset = _preset_name(bundle, filter_settings) if filter_settings else "not applied"
    threshold_rows = []
    for key, label, direction in (
        ("score", "Overall Score", "≥"),
        ("fragmentation", "Fragmentation Score", "≥"),
        ("abs_mass_error", "Absolute mass error", "≤"),
        ("isotope", "Isotope Similarity", "≥"),
    ):
        recommended = bundle.thresholds.get(key, {})
        active_value = filter_settings.get(key, "not applied") if filter_settings else "not applied"
        unit = " ppm" if key == "abs_mass_error" else ""
        interval = (
            f"{recommended.get('ci_low'):g}–{recommended.get('ci_high'):g}{unit}"
            if recommended.get("ci_low") is not None and recommended.get("ci_high") is not None
            else "not estimable"
        )
        threshold_rows.append([label, f"{direction} {active_value}{unit}", interval, recommended.get("method", "not estimated")])
    threshold_table = _markdown_table(["Metric", "Active criterion", "Recommended-cutoff bootstrap interval", "Recommendation method"], threshold_rows)
    diagnostics = bundle.thresholds.get("diagnostics", {})
    mass_statement = (
        f"The absolute-mass-error distribution supported a two-component boundary at {float(diagnostics.get('mass_intersection')):.3f} ppm "
        f"(ΔBIC {float(diagnostics.get('mass_delta_bic')):.2f}; Ashman D {float(diagnostics.get('mass_ashman_d')):.2f})."
        if diagnostics.get("mass_mixture_used")
        else "The absolute-mass-error distribution did not satisfy the predeclared two-component separation criteria; a conservative robust bound was used."
    )
    sections.append(
        {
            "section": "Methods",
            "subsection": "Dataset-adaptive identification-evidence thresholds",
            "text": (
                f"Threshold recommendations were estimated after identity collapse, preventing multiply detected entities from overweighting the empirical distributions. The active identity setting was {preset}. "
                "For the balanced strategy, overall Score and Isotope Similarity used entity-level first quartiles; Fragmentation Score used the first quartile among positive values, while zero could be retained as an explicit unavailable/disabled/unmatched state; "
                "and absolute mass error used a deterministic two-component boundary only when ΔBIC > 10 and Ashman D ≥ 2, otherwise a robust upper bound capped by the declared search tolerance. "
                "Cutoff stability was summarized using 500 deterministic nonparametric bootstrap resamples. "
                f"{mass_statement} The active criteria were:\n\n{threshold_table}\n\n"
                "The overall Progenesis Score combines the enabled similarity dimensions; unavailable RT, CCS, or fragmentation dimensions can contribute zero, so these empirical cutoffs are dataset-specific triage rules rather than universal identification thresholds [5]."
            ),
        }
    )

    qc_available = qc_n >= 2
    blank_available = blank_n >= 1
    drift_available = qc_n >= 4 and included.loc[included["Sample role"].eq("Pooled QC"), "Injection order"].notna().sum() >= 4
    d_ratio_available = qc_n >= 2 and biological_n >= 2
    ms2_available = ms2_supported > 0
    qc_rule_rows = [
        ["Pooled-QC detection", _rule_state(bool(filter_settings.get("apply_qc_filter")), qc_n >= 1), f"≥ {float(filter_settings.get('min_qc_detection', 0)):.0%}"],
        ["Pooled-QC CV", _rule_state(bool(filter_settings.get("apply_qc_filter")), qc_available), f"≤ {float(filter_settings.get('max_qc_cv', np.nan)):g}%"],
        ["Biological detection", _rule_state(True, biological_n >= 1), f"≥ {float(filter_settings.get('min_biological_detection', 0)):.0%}"],
        ["Biological/blank median ratio", _rule_state(bool(filter_settings.get("apply_blank_filter")), blank_available and biological_n >= 1), f"≥ {float(filter_settings.get('min_biological_blank_ratio', np.nan)):g}"],
        ["Absolute across-run QC drift", _rule_state(bool(filter_settings.get("apply_drift_filter")), drift_available), f"≤ {float(filter_settings.get('max_abs_qc_drift', np.nan)):g}%"],
        ["D-ratio", _rule_state(bool(filter_settings.get("apply_d_ratio_filter")), d_ratio_available), f"≤ {float(filter_settings.get('max_d_ratio', np.nan)):g}%"],
        ["Linked MS² peak count", _rule_state(bool(filter_settings.get("require_ms2")), ms2_available), f"≥ {int(filter_settings.get('min_fragment_peaks', 0))}"],
    ]
    sections.append(
        {
            "section": "Methods",
            "subsection": "Analytical-QC metrics and feature filtering",
            "text": (
                "Positive intensities were used for abundance summaries. Detection rate was the fraction of nonmissing measurements > 0. Conventional CV was 100×SD/mean and robust CV was "
                "100×1.4826×MAD/median. The dispersion ratio was D-ratio = 100×SD(QC)/√[SD(QC)² + SD(biological)²], where lower values indicate that technical dispersion contributes less to total dispersion [4]. "
                "Biological/blank ratio was the ratio of positive medians. Across-run drift was estimated by linear regression of log₂ pooled-QC intensity on mapped injection order and expressed as "
                "100×[2^(slope×order span)−1], with Spearman rank correlation reported as a monotonicity diagnostic. Technical-replicate CV and dilution-response Spearman rho were reported when their metadata requirements were met. "
                "Missing evidence was classified as not evaluable and could not create an observed pass. Active states were:\n\n"
                + _markdown_table(["Rule", "State", "Criterion"], qc_rule_rows)
            ),
        }
    )

    _emit(progress, 63, "Writing dataset-specific curation analysis…")
    if not decisions.empty:
        id_pass = int(decisions["Identification evidence pass"].sum())
        qc_pass = int(decisions["Analytical QC pass"].sum())
        final_pass = int(decisions["Filter pass"].sum())
        id_failures = _failure_counts(decisions, "Identification filter fail reasons")
        qc_failures = _failure_counts(decisions, "Analytical QC fail reasons")
        id_failure_text = "; ".join(f"{reason}: {count:,}" for reason, count in id_failures.most_common()) or "none"
        qc_failure_text = "; ".join(f"{reason}: {count:,}" for reason, count in qc_failures.most_common()) or "none"
    else:
        id_pass = qc_pass = final_pass = 0
        id_failure_text = qc_failure_text = "filters not applied"
    sections.append(
        {
            "section": "Curation analysis",
            "subsection": "Data integrity and attrition",
            "text": (
                f"All six sources were structurally readable. The minimum CM-key coverage across sources was {pd.to_numeric(bundle.audits['CM coverage (%)'], errors='coerce').min():.3f}%. "
                f"From {len(bundle.accepted_reads):,} accepted feature reads, identity collapse yielded {len(bundle.selected_reads):,} representatives. "
                f"The identity-evidence criteria retained {id_pass:,}/{len(bundle.selected_reads):,} entities ({_pct(id_pass, len(bundle.selected_reads))}); "
                f"{qc_pass:,}/{len(bundle.selected_reads):,} did not fail an enabled, evaluable analytical-QC rule; and their intersection yielded {final_pass:,} shortlisted compounds "
                f"({_pct(final_pass, len(bundle.selected_reads))} of representatives; {_pct(final_pass, len(bundle.accepted_reads))} of accepted feature reads). "
                f"Identity-rule failures were: {id_failure_text}. Analytical-QC failures were: {qc_failure_text}. Because rules overlap, failure counts are not additive."
            ),
        }
    )
    sections.append(
        {
            "section": "Curation analysis",
            "subsection": "Evidence profile of the final shortlist",
            "text": (
                f"Within the {len(shortlist):,}-compound shortlist, overall Score was {_quantiles(shortlist, 'Score')}, Fragmentation Score was {_quantiles(shortlist, 'Fragmentation Score')}, "
                f"absolute mass error was {_quantiles(shortlist, 'Absolute Mass Error (ppm)', absolute=True)}, and Isotope Similarity was {_quantiles(shortlist, 'Isotope Similarity')}. "
                f"Pooled-QC CV was {_quantiles(shortlist, 'QC CV%')}; biological detection was {_quantiles(shortlist, 'Biological detection rate')}. "
                f"Linked fragment peaks were present for {int(_numeric_series(shortlist, 'FD observed peaks').fillna(0).gt(0).sum()):,}/{len(shortlist):,} shortlisted compounds. "
                "These summaries characterize the selected evidence distribution and must not be interpreted as an identification error rate."
            ),
        }
    )

    _emit(progress, 79, "Summarizing source and extraction-context outputs…")
    if provenance.empty:
        context_text = "Contextual source and extraction-location scoring was not run; no source-origin or phase summary was generated."
    else:
        source_counts = provenance["Primary source class"].value_counts()
        source_text = ", ".join(
            f"{label}: {int(source_counts.get(label, 0)):,} ({_pct(int(source_counts.get(label, 0)), len(provenance))})"
            for label in (
                "Human endogenous",
                "Microbiome-derived",
                "Host–microbe co-metabolic",
                "Food-derived",
                "Drug-derived",
                "Environmental-derived",
            )
            if int(source_counts.get(label, 0))
        )
        confidence_counts = provenance["Source confidence"].value_counts()
        confidence_text = ", ".join(f"{key}: {int(value):,}" for key, value in confidence_counts.items())
        phase_counts = provenance["Dominant predicted phase"].value_counts()
        phase_text = ", ".join(f"{key}: {int(value):,}" for key, value in phase_counts.items())
        analyzed_values = _numeric_series(provenance, "Analyzed-phase likelihood (%)")
        analyzed_values = analyzed_values[np.isfinite(analyzed_values)]
        analyzed_text = (
            f"median {analyzed_values.median():.1f}% (IQR {analyzed_values.quantile(.25):.1f}–{analyzed_values.quantile(.75):.1f}%)"
            if not analyzed_values.empty
            else "not estimable"
        )
        context_text = (
            f"Using the declared {_humanize(parameters.get('assay_type'))}, {_humanize(parameters.get('biological_system'))}, {_humanize(parameters.get('sample_matrix'))}, "
            f"{_humanize(parameters.get('extraction_method'))}, and {_humanize(parameters.get('analyzed_phase'))} context, primary source classes were {source_text}. "
            f"Source-confidence labels were {confidence_text}. Dominant extraction-location predictions were {phase_text}; analyzed-phase likelihood was {analyzed_text}. "
            "The six source values are temperature-scaled softmax-normalized rule scores combining declared context with description/ontology cues, including MiMeDB relationships [9]. "
            "The phase model combines explicit description labels, parsed MolLogP/TPSA when present, formula-derived polarity/lipid cues, and method-specific extraction logic. "
            "For MTBE, the organic-rich phase is expected to be the upper layer [8]. These outputs are hypothesis-ranking scores, not calibrated causal-origin probabilities, measured partition coefficients, or extraction recoveries."
        )
    sections.append({"section": "Curation analysis", "subsection": "Contextual source and extraction interpretation", "text": context_text})

    warnings = " ".join(f"{index + 1}) {warning}" for index, warning in enumerate(bundle.warnings)) or "No automated dataset warning was generated."
    sections.append(
        {
            "section": "Curation analysis",
            "subsection": "Identification confidence and limitations",
            "text": (
                "The shortlist is a transparent curation result, not an MSI Level-1 identification set and not a target–decoy-controlled FDR. Exact mass, formula, adduct, isotope similarity, and vendor fragmentation scores do not by themselves resolve structural or stereochemical isomers. "
                "Repeated-ID alternatives, source links, selection traces, raw/normalized intensities, and all rule-level decisions should accompany downstream statistics. Dataset-specific limitations were: "
                + warnings
            ),
        }
    )

    reporting_items = [
        "LC column chemistry, dimensions, particle size, guard column, mobile-phase compositions/additives, gradient, flow rate, column temperature, injection volume, and autosampler temperature.",
        "Mass-spectrometer manufacturer/model, source settings, mass-calibration and lock-mass procedure, scan range/rate, DIA/MSE collision-energy functions, and acquisition software version.",
        "Sample mass/volume, solvent identities and ratios, extraction repetitions, mixing/sonication, temperature, centrifugation, phase volumes, drying/reconstitution, storage, and freeze–thaw history.",
        "Run-order randomization/blocking, QC-conditioning injections and frequency, carryover strategy, system-suitability criteria, and batch-correction or normalization procedure used upstream of these exports.",
        "Blank design and results, authentic-standard/retention-time confirmation, spectral-library name/version/search parameters, target–decoy strategy, CCS reference, and identification-level nomenclature.",
        "Downstream missing-value handling, transformation, normalization, batch correction, univariate/multivariate models, covariates, multiple-testing correction, effect sizes, confidence intervals, and software versions.",
    ]
    sections.append(
        {
            "section": "Journal completion checklist",
            "subsection": "Information not recoverable from the six exports",
            "text": (
                "The following information must be completed from laboratory records and acquisition/processing methods before manuscript submission:\n\n"
                + "\n".join(f"- {item}" for item in reporting_items)
                + "\n\nDo not infer or fabricate these items from feature tables. This checklist is aligned with MSI and Lipidomics Minimal Reporting Checklist expectations [1,2]."
            ),
        }
    )

    _emit(progress, 94, "Adding primary references and reproducibility statement…")
    reference_text = "\n".join(f"{number}. {citation} {url}" for number, citation, url in REFERENCES)
    sections.append({"section": "References", "subsection": "Primary methods and reporting sources", "text": reference_text})
    sections.append(
        {
            "section": "Reproducibility",
            "subsection": "Audit trail",
            "text": (
                "The complete export contains input SHA-256 hashes, project and sample metadata, raw source summaries, every accepted and repeated read, spectrum- and peak-level MS² tables, the representative-selection trace, adaptive-threshold estimates and bootstrap intervals, "
                "entity-level QC/filter decisions, the shortlist, contextual scores, a machine-readable manifest, and this dataset-specific report. These files define the computational curation state but do not replace the original raw instrument files or acquisition method."
            ),
        }
    )
    _emit(progress, 100, "Journal-ready methodology and curation analysis complete")
    return sections


def journal_report_markdown(
    bundle: AnalysisBundle,
    shortlist: pd.DataFrame,
    provenance: pd.DataFrame,
    filter_settings: dict[str, Any],
    provenance_settings: dict[str, Any] | None,
    progress: ProgressCallback | None = None,
) -> str:
    sections = build_journal_sections(bundle, shortlist, provenance, filter_settings, provenance_settings, progress)
    blocks = [
        "# Manuscript-ready LC–MS curation methodology and analysis",
        "",
        "> Dataset-specific text generated from the selected project options, audited pipeline state, active thresholds, mapped samples, and final shortlisted compounds. Verify and merge this section with the experimental LC–MS method before submission.",
    ]
    current_section = ""
    for item in sections:
        if item["section"] != current_section:
            blocks.extend(["", f"## {item['section']}"])
            current_section = item["section"]
        blocks.extend(["", f"### {item['subsection']}", "", item["text"]])
    return "\n".join(blocks).strip() + "\n"


def journal_report_table(
    bundle: AnalysisBundle,
    shortlist: pd.DataFrame,
    provenance: pd.DataFrame,
    filter_settings: dict[str, Any],
    provenance_settings: dict[str, Any] | None,
) -> pd.DataFrame:
    return pd.DataFrame(
        build_journal_sections(bundle, shortlist, provenance, filter_settings, provenance_settings),
        columns=["section", "subsection", "text"],
    ).rename(columns={"section": "Section", "subsection": "Subsection", "text": "Manuscript-ready text"})
