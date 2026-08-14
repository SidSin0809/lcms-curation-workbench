from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


BLANK_ROLES = {"Process blank", "Extraction blank", "Solvent blank"}


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Return one-based average ranks with deterministic handling of ties."""
    numeric = np.asarray(values, dtype=float)
    order = np.argsort(numeric, kind="mergesort")
    sorted_values = numeric[order]
    ranks = np.empty(len(numeric), dtype=float)
    start = 0
    while start < len(numeric):
        stop = start + 1
        while stop < len(numeric) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        average_rank = ((start + 1) + stop) / 2.0
        ranks[order[start:stop]] = average_rank
        start = stop
    return ranks


def _spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation without pandas' optional SciPy dependency.

    Spearman's rho is the Pearson correlation of paired average ranks. Missing
    or infinite pairs are removed together; a constant rank vector is not
    informative and returns NaN.
    """
    x_values = np.asarray(x, dtype=float).reshape(-1)
    y_values = np.asarray(y, dtype=float).reshape(-1)
    if x_values.size != y_values.size:
        raise ValueError("Spearman inputs must have equal length")
    mask = np.isfinite(x_values) & np.isfinite(y_values)
    if int(mask.sum()) < 2:
        return float("nan")
    x_ranks = _average_ranks(x_values[mask])
    y_ranks = _average_ranks(y_values[mask])
    x_centered = x_ranks - x_ranks.mean()
    y_centered = y_ranks - y_ranks.mean()
    denominator = float(np.linalg.norm(x_centered) * np.linalg.norm(y_centered))
    if denominator == 0.0:
        return float("nan")
    rho = float(np.dot(x_centered, y_centered) / denominator)
    return float(np.clip(rho, -1.0, 1.0))


def _positive_values(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(pd.to_numeric, errors="coerce").where(lambda values: values > 0)


def _cv(values: pd.DataFrame) -> pd.Series:
    return values.std(axis=1, ddof=1) / values.mean(axis=1).abs() * 100


def _robust_cv(values: pd.DataFrame) -> pd.Series:
    median = values.median(axis=1)
    deviation = values.sub(median, axis=0).abs().median(axis=1)
    return 1.4826 * deviation / median.abs() * 100


def _detection(values: pd.DataFrame) -> pd.Series:
    numeric = values.apply(pd.to_numeric, errors="coerce")
    return numeric.gt(0).sum(axis=1) / numeric.notna().sum(axis=1).clip(lower=1)


def _role_columns(sample_metadata: pd.DataFrame, roles: set[str]) -> list[str]:
    if sample_metadata.empty:
        return []
    included = sample_metadata["Include"].astype(bool)
    return sample_metadata.loc[included & sample_metadata["Sample role"].isin(roles), "Normalized abundance column"].tolist()


def _technical_replicate_cv(frame: pd.DataFrame, metadata: pd.DataFrame) -> pd.Series:
    biological = metadata.loc[
        metadata["Include"].astype(bool) & metadata["Sample role"].eq("Biological") & metadata["Technical replicate"].astype(bool)
    ]
    groups: list[list[str]] = []
    for _, group in biological.groupby("Subject / biological unit", sort=False):
        columns = group["Normalized abundance column"].tolist()
        if len(columns) >= 2:
            groups.append(columns)
    if not groups:
        return pd.Series(np.nan, index=frame.index)
    group_cvs = pd.concat([_cv(_positive_values(frame[columns])) for columns in groups], axis=1)
    return group_cvs.median(axis=1)


def compute_abundance_metrics(frame: pd.DataFrame, sample_metadata: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    qc_columns = [column for column in _role_columns(sample_metadata, {"Pooled QC"}) if column in output.columns]
    bio_columns = [column for column in _role_columns(sample_metadata, {"Biological"}) if column in output.columns]
    blank_columns = [column for column in _role_columns(sample_metadata, BLANK_ROLES) if column in output.columns]

    def values(columns: list[str]) -> pd.DataFrame:
        return _positive_values(output[columns]) if columns else pd.DataFrame(index=output.index)

    qc = values(qc_columns)
    bio = values(bio_columns)
    blank = values(blank_columns)
    output["Pooled QC sample count"] = len(qc_columns)
    output["Biological sample count"] = len(bio_columns)
    output["Blank sample count"] = len(blank_columns)
    output["QC detection rate"] = _detection(output[qc_columns]) if qc_columns else np.nan
    output["QC mean abundance"] = qc.mean(axis=1) if qc_columns else np.nan
    output["QC median abundance"] = qc.median(axis=1) if qc_columns else np.nan
    output["QC CV%"] = _cv(qc) if len(qc_columns) >= 2 else np.nan
    output["QC robust CV% (MAD)"] = _robust_cv(qc) if len(qc_columns) >= 3 else np.nan
    output["Biological detection rate"] = _detection(output[bio_columns]) if bio_columns else np.nan
    output["Biological mean abundance"] = bio.mean(axis=1) if bio_columns else np.nan
    output["Biological median abundance"] = bio.median(axis=1) if bio_columns else np.nan
    output["Biological CV%"] = _cv(bio) if len(bio_columns) >= 2 else np.nan
    output["Blank detection rate"] = _detection(output[blank_columns]) if blank_columns else np.nan
    output["Blank median abundance"] = blank.median(axis=1) if blank_columns else np.nan
    output["Blank maximum abundance"] = blank.max(axis=1) if blank_columns else np.nan
    if qc_columns and bio_columns:
        qc_sd = qc.std(axis=1, ddof=1)
        bio_sd = bio.std(axis=1, ddof=1)
        output["D-ratio (%)"] = qc_sd / np.sqrt(qc_sd.pow(2) + bio_sd.pow(2)) * 100
    else:
        output["D-ratio (%)"] = np.nan
    if blank_columns and bio_columns:
        blank_median = output["Blank median abundance"].fillna(0)
        biological_median = output["Biological median abundance"].fillna(0)
        output["Biological/blank median ratio"] = np.where(
            blank_median > 0,
            biological_median / blank_median,
            np.where(biological_median > 0, np.inf, np.nan),
        )
        output["Blank contribution (%)"] = blank_median / np.maximum(blank_median, biological_median) * 100
    else:
        output["Biological/blank median ratio"] = np.nan
        output["Blank contribution (%)"] = np.nan
    output["Technical replicate median CV%"] = _technical_replicate_cv(output, sample_metadata)

    qc_order = sample_metadata.loc[
        sample_metadata["Include"].astype(bool) & sample_metadata["Sample role"].eq("Pooled QC"),
        ["Normalized abundance column", "Injection order"],
    ]
    drift = pd.Series(np.nan, index=output.index)
    drift_rho = pd.Series(np.nan, index=output.index)
    if len(qc_order) >= 4:
        columns = [column for column in qc_order["Normalized abundance column"] if column in output.columns]
        orders = pd.to_numeric(qc_order.set_index("Normalized abundance column").loc[columns, "Injection order"], errors="coerce")
        if orders.notna().sum() >= 4:
            x = orders.to_numpy(float)
            span = float(np.nanmax(x) - np.nanmin(x))
            for index, row in _positive_values(output[columns]).iterrows():
                y = np.log2(row.to_numpy(float))
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() >= 4 and span > 0:
                    slope = float(np.polyfit(x[mask], y[mask], 1)[0])
                    drift.at[index] = (2 ** (slope * span) - 1) * 100
                    drift_rho.at[index] = _spearman_rho(x[mask], y[mask])
    output["QC drift across run (%)"] = drift
    output["QC injection-order Spearman rho"] = drift_rho

    dilution = sample_metadata.loc[
        sample_metadata["Include"].astype(bool) & sample_metadata["Sample role"].isin({"Pooled QC", "Calibration"}),
        ["Normalized abundance column", "Dilution factor"],
    ]
    dilution_rho = pd.Series(np.nan, index=output.index)
    if dilution["Dilution factor"].nunique(dropna=True) >= 3:
        columns = [column for column in dilution["Normalized abundance column"] if column in output.columns]
        x = pd.to_numeric(dilution.set_index("Normalized abundance column").loc[columns, "Dilution factor"], errors="coerce")
        for index, row in _positive_values(output[columns]).iterrows():
            values = row.to_numpy(float)
            mask = np.isfinite(x.to_numpy(float)) & np.isfinite(values)
            if mask.sum() >= 3:
                dilution_rho.at[index] = _spearman_rho(x.to_numpy(float)[mask], values[mask])
    output["Dilution-response Spearman rho"] = dilution_rho
    return output


def default_qc_settings() -> dict[str, Any]:
    return {
        "min_qc_detection": 0.80,
        "max_qc_cv": 30.0,
        "min_biological_detection": 0.20,
        "min_biological_blank_ratio": 3.0,
        "max_abs_qc_drift": 40.0,
        "max_d_ratio": 50.0,
        "require_ms2": False,
        "min_fragment_peaks": 3,
        "apply_blank_filter": True,
        "apply_qc_filter": True,
        "apply_drift_filter": False,
        "apply_d_ratio_filter": False,
    }


def evaluate_analytical_qc(frame: pd.DataFrame, settings: dict[str, Any] | None = None) -> pd.DataFrame:
    settings = {**default_qc_settings(), **(settings or {})}
    output = frame.copy()
    decisions: list[str] = []
    reasons: list[str] = []
    coverage: list[str] = []
    for _, row in output.iterrows():
        failures: list[str] = []
        available: list[str] = []
        if settings["apply_qc_filter"] and pd.notna(row.get("QC detection rate")):
            available.append("pooled-QC detection")
            if float(row["QC detection rate"]) < settings["min_qc_detection"]:
                failures.append(f"QC detection < {settings['min_qc_detection']:.0%}")
        if settings["apply_qc_filter"] and pd.notna(row.get("QC CV%")):
            available.append("pooled-QC CV")
            if float(row["QC CV%"]) > settings["max_qc_cv"]:
                failures.append(f"QC CV > {settings['max_qc_cv']:g}%")
        if pd.notna(row.get("Biological detection rate")):
            available.append("biological detection")
            if float(row["Biological detection rate"]) < settings["min_biological_detection"]:
                failures.append(f"biological detection < {settings['min_biological_detection']:.0%}")
        if settings["apply_blank_filter"] and pd.notna(row.get("Biological/blank median ratio")):
            available.append("blank contribution")
            if float(row["Biological/blank median ratio"]) < settings["min_biological_blank_ratio"]:
                failures.append(f"sample/blank ratio < {settings['min_biological_blank_ratio']:g}")
        if settings["apply_drift_filter"] and pd.notna(row.get("QC drift across run (%)")):
            available.append("QC drift")
            if abs(float(row["QC drift across run (%)"])) > settings["max_abs_qc_drift"]:
                failures.append(f"absolute QC drift > {settings['max_abs_qc_drift']:g}%")
        if settings["apply_d_ratio_filter"] and pd.notna(row.get("D-ratio (%)")):
            available.append("D-ratio")
            if float(row["D-ratio (%)"]) > settings["max_d_ratio"]:
                failures.append(f"D-ratio > {settings['max_d_ratio']:g}%")
        if settings["require_ms2"]:
            available.append("MS2 peak support")
            if float(row.get("FD observed peaks") or 0) < settings["min_fragment_peaks"]:
                failures.append(f"fragment peaks < {settings['min_fragment_peaks']}")
        decisions.append("Fail" if failures else "Pass" if available else "Not evaluable")
        reasons.append("; ".join(failures) if failures else "No enabled analytical-QC rule failed" if available else "No mapped QC/blank evidence")
        coverage.append("; ".join(available) if available else "None")
    output["Analytical QC decision"] = decisions
    output["Analytical QC fail reasons"] = reasons
    output["Analytical QC evidence available"] = coverage
    output["Analytical QC pass"] = output["Analytical QC decision"].ne("Fail")
    return output
