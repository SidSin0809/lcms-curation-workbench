from __future__ import annotations

import csv
import hashlib
import math
import re
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

from .chemistry import enrich_chemical_properties
from .metadata import build_sample_metadata, normalize_sample_metadata
from .qc import compute_abundance_metrics, evaluate_analytical_qc

FILE_ROLES = ("CM", "AM", "CI", "II", "FD", "ACP")

ROLE_LABELS = {
    "CM": "Compound measurements (three-row header CSV)",
    "AM": "Adduct measurements (three-row header CSV)",
    "CI": "Compound identifications (all candidates CSV)",
    "II": "Isotope/ion information (XML)",
    "FD": "Fragmentation data (MSP)",
    "ACP": "Accepted compound properties (CSV)",
}

REQUIRED_FIELDS: dict[str, list[str]] = {
    "CM": [
        "Metadata|Compound",
        "Metadata|m/z",
        "Metadata|Retention time (min)",
        "Tags|Accepted Compound ID",
        "Tags|Accepted Description",
        "Tags|Score",
        "Tags|Fragmentation Score",
        "Tags|Mass Error (ppm)",
        "Tags|Isotope Similarity",
    ],
    "AM": [
        "Metadata|Compound",
        "Metadata|Ion",
        "Metadata|m/z",
        "Metadata|Charge",
        "Metadata|Adduct type",
    ],
    "CI": [
        "Compound",
        "Compound ID",
        "Accepted?",
        "Adducts",
        "Formula",
        "Score",
        "Fragmentation Score",
        "Mass Error (ppm)",
        "Isotope Similarity",
        "Description",
    ],
    "II": ["compound@identifier", "condition/sample/adduct/isotope"],
    "FD": ["DATABASE_ID", "Comment", "Num Peaks", "peak pairs"],
    "ACP": ["Compound", "Compound ID", "Retention time (min)"],
}


@dataclass(slots=True)
class AnalysisBundle:
    generated_at: str
    parameters: dict[str, Any]
    files: dict[str, Path]
    audits: pd.DataFrame
    warnings: list[str]
    cm: pd.DataFrame
    am: pd.DataFrame
    ci: pd.DataFrame
    accepted_ci: pd.DataFrame
    acp: pd.DataFrame
    ii: pd.DataFrame
    fd: pd.DataFrame
    fd_peaks: pd.DataFrame
    sample_metadata: pd.DataFrame
    accepted_reads: pd.DataFrame
    duplicate_reads: pd.DataFrame
    selected_reads: pd.DataFrame
    analytical_qc: pd.DataFrame
    duplicate_groups: int
    extra_reads_collapsed: int
    thresholds: dict[str, Any]
    normalized_columns: list[str]
    raw_columns: list[str]


def clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return re.sub(r"\s+", " ", str(value)).strip()


def clean_link(value: Any) -> str:
    """Return only usable source links; vendor NULL/NA URL placeholders become missing."""
    link = clean(value)
    if not link:
        return ""
    terminal = link.rstrip("/").rsplit("/", 1)[-1].casefold()
    if terminal in {"null", "none", "na", "n/a", "nan", "unknown"}:
        return ""
    return link


def finite(value: Any) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def safe_join(values: Iterable[Any]) -> str:
    return "; ".join(dict.fromkeys(text for value in values if (text := clean(value))))


def canonical_name(description: str) -> str:
    text = clean(description)
    if not text:
        return ""
    text = re.split(r"\s+\([A-Z][A-Za-z0-9]*; represented-formula", text, maxsplit=1, flags=re.IGNORECASE)[0]
    text = text.split(" is included because ", 1)[0]
    return text.rstrip(". ")


def _unique_columns(columns: Iterable[str]) -> list[str]:
    seen: Counter[str] = Counter()
    output: list[str] = []
    for column in columns:
        seen[column] += 1
        output.append(column if seen[column] == 1 else f"{column} [{seen[column]}]")
    return output


def _read_text_encoding(path: Path) -> str:
    raw = path.read_bytes()[:65536]
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return "utf-16"
    try:
        raw.decode("utf-8-sig")
        return "utf-8-sig"
    except UnicodeDecodeError:
        return "cp1252"


def read_three_header_csv(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    encoding = _read_text_encoding(path)
    with path.open("r", encoding=encoding, newline="", errors="replace") as handle:
        reader = csv.reader(handle)
        try:
            header_rows = [next(reader) for _ in range(3)]
        except StopIteration as exc:
            raise ValueError(f"{path.name} requires three header rows plus data rows.") from exc

    width = max(map(len, header_rows))
    header_rows = [row + [""] * (width - len(row)) for row in header_rows]
    section = "Metadata"
    condition = ""
    keys: list[str] = []
    structure: list[dict[str, Any]] = []
    normalized: list[str] = []
    raw: list[str] = []

    for index in range(width):
        top = clean(header_rows[0][index])
        middle = clean(header_rows[1][index])
        leaf = clean(header_rows[2][index]) or f"Unnamed column {index + 1}"
        if top:
            section = top
            condition = ""
        if middle:
            if middle.casefold() == "tags":
                section = "Tags"
                condition = ""
            elif section in {"Normalised abundance", "Raw abundance"}:
                condition = middle
        if section in {"Normalised abundance", "Raw abundance"}:
            key = f"{section}|{condition or 'Unassigned'}|{leaf}"
            (normalized if section == "Normalised abundance" else raw).append(key)
        elif section == "Tags":
            key = f"Tags|{leaf}"
        else:
            key = f"Metadata|{leaf}"
        keys.append(key)
        structure.append({"index": index, "section": section, "condition": condition, "label": leaf, "key": key})

    keys = _unique_columns(keys)
    frame = pd.read_csv(
        path,
        skiprows=3,
        header=None,
        names=keys,
        dtype=str,
        encoding=encoding,
        keep_default_na=False,
        na_filter=False,
        low_memory=False,
    )
    return frame, {"columns": structure, "normalized": normalized, "raw": raw, "encoding": encoding}


def _detect_header_row(path: Path, required: set[str], maximum: int = 8) -> tuple[int, str]:
    encoding = _read_text_encoding(path)
    with path.open("r", encoding=encoding, newline="", errors="replace") as handle:
        reader = csv.reader(handle)
        for index, row in enumerate(reader):
            values = {clean(value) for value in row}
            if required.issubset(values):
                return index, encoding
            if index + 1 >= maximum:
                break
    raise ValueError(f"Could not locate the expected CSV header in {path.name}.")


def read_standard_csv(path: Path, required: set[str]) -> pd.DataFrame:
    header_row, encoding = _detect_header_row(path, required)
    return pd.read_csv(
        path,
        skiprows=header_row,
        dtype=str,
        encoding=encoding,
        keep_default_na=False,
        na_filter=False,
        low_memory=False,
    )


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_isotope_xml(path: Path, progress: Callable[[int, str], None] | None = None) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, compound in ET.iterparse(path, events=("end",)):
        if _local_name(compound.tag) != "compound":
            continue
        condition_nodes = [node for node in list(compound) if _local_name(node.tag) == "condition"]
        sample_nodes = [sample for condition in condition_nodes for sample in list(condition) if _local_name(sample.tag) == "sample"]
        adduct_nodes = [adduct for sample in sample_nodes for adduct in list(sample) if _local_name(adduct.tag) == "adduct"]
        isotope_nodes = [isotope for adduct in adduct_nodes for isotope in list(adduct) if _local_name(isotope.tag) == "isotope"]
        isotope_mz: list[float] = []
        isotope_abundance: list[float] = []
        for isotope in isotope_nodes:
            for child in list(isotope):
                value = finite(child.text)
                if value is None:
                    continue
                if _local_name(child.tag) == "mz":
                    isotope_mz.append(value)
                elif _local_name(child.tag) == "abundance":
                    isotope_abundance.append(value)
        normalized = [value for node in sample_nodes if (value := finite(node.get("normalizedAbundance"))) is not None]
        statistics = next((node for node in list(compound) if _local_name(node.tag) == "statistics"), None)
        stat_values: dict[str, str] = {}
        mean_lowest = ""
        mean_highest = ""
        if statistics is not None:
            for child in list(statistics):
                name = _local_name(child.tag)
                stat_values[name] = clean(child.text)
                if name == "mean":
                    mean_lowest = clean(child.get("lowest"))
                    mean_highest = clean(child.get("highest"))
        rows.append(
            {
                "Feature ID": clean(compound.get("identifier")),
                "II retention time": finite(compound.get("retentionTime")),
                "II neutral mass": finite(compound.get("neutralMass")),
                "II ANOVA p": finite(stat_values.get("anova")),
                "II max fold change": finite(stat_values.get("maxFoldChange")),
                "II lowest mean": mean_lowest,
                "II highest mean": mean_highest,
                "II conditions": safe_join(node.get("name") for node in condition_nodes),
                "II sample count": len(sample_nodes),
                "II distinct samples": len({clean(node.get("name")) for node in sample_nodes if clean(node.get("name"))}),
                "II adduct-node count": len(adduct_nodes),
                "II adduct descriptions": safe_join(node.get("description") for node in adduct_nodes),
                "II isotope-node count": len(isotope_nodes),
                "II isotope m/z min": min(isotope_mz) if isotope_mz else np.nan,
                "II isotope m/z max": max(isotope_mz) if isotope_mz else np.nan,
                "II isotope abundance sum": sum(isotope_abundance),
                "II normalized abundance sum": sum(normalized),
                "II nonzero normalized samples": sum(value > 0 for value in normalized),
            }
        )
        compound.clear()
        if progress and len(rows) % 250 == 0:
            progress(min(58, 46 + len(rows) // 125), f"Scanning isotope XML… {len(rows):,} compounds")
    if not rows:
        raise ValueError("II XML contains no compound nodes with an identifier.")
    return pd.DataFrame(rows)


def parse_msp_detailed(
    path: Path,
    progress: Callable[[int, str], None] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    encoding = _read_text_encoding(path)
    records: list[dict[str, Any]] = []
    peak_rows: list[dict[str, Any]] = []
    metadata: dict[str, str] = {}
    peaks: list[tuple[float, float]] = []

    def emit() -> None:
        nonlocal metadata, peaks
        if not metadata and not peaks:
            return
        declared = finite(metadata.get("Num Peaks"))
        base_peak = max(peaks, key=lambda pair: pair[1]) if peaks else (np.nan, np.nan)
        feature_id = clean(metadata.get("Comment"))
        spectrum_index = len(records) + 1
        precursor = finite(metadata.get("PrecursorMZ"))
        positive_intensities = np.asarray([max(0.0, pair[1]) for pair in peaks], dtype=float)
        total_intensity = float(positive_intensities.sum())
        probabilities = positive_intensities / total_intensity if total_intensity > 0 else np.asarray([], dtype=float)
        entropy = float(-(probabilities[probabilities > 0] * np.log(probabilities[probabilities > 0])).sum()) if len(probabilities) else np.nan
        normalized_entropy = entropy / math.log(len(probabilities)) if len(probabilities) > 1 else 0.0 if len(probabilities) == 1 else np.nan
        base_intensity = float(base_peak[1]) if peaks else np.nan
        low_intensity_share = (
            float(sum(intensity < 0.01 * base_intensity for _, intensity in peaks) / len(peaks) * 100)
            if peaks and base_intensity > 0
            else np.nan
        )
        top_five_share = (
            float(np.sort(positive_intensities)[-5:].sum() / total_intensity * 100)
            if total_intensity > 0
            else np.nan
        )
        integrity = bool(declared is not None and int(declared) == len(peaks))
        quality_flag = "Fail" if not peaks else "Review" if not integrity or len(peaks) < 3 else "Pass"
        records.append(
            {
                "FD database ID": clean(metadata.get("DATABASE_ID")),
                "Feature ID": feature_id,
                "FD spectrum index": spectrum_index,
                "FD name": clean(metadata.get("Name")),
                "FD precursor type": clean(metadata.get("Precursor_type")),
                "FD precursor m/z": precursor,
                "FD formula": clean(metadata.get("Formula")),
                "FD declared peaks": int(declared) if declared is not None else np.nan,
                "FD observed peaks": len(peaks),
                "FD peak-count integrity": integrity,
                "FD fragment m/z min": min((pair[0] for pair in peaks), default=np.nan),
                "FD fragment m/z max": max((pair[0] for pair in peaks), default=np.nan),
                "FD base peak m/z": base_peak[0],
                "FD base peak intensity": base_peak[1],
                "FD total fragment intensity": total_intensity,
                "FD spectral entropy": entropy,
                "FD normalized spectral entropy": normalized_entropy,
                "FD effective peak count": math.exp(entropy) if math.isfinite(entropy) else np.nan,
                "FD top-five intensity share (%)": top_five_share,
                "FD low-intensity peak share (%)": low_intensity_share,
                "FD spectrum quality flag": quality_flag,
                "FD spectrum-quality limitation": "Entropy and peak counts describe spectrum complexity; without a reference/decoy library they do not establish identity or FDR.",
            }
        )
        for peak_index, (fragment_mz, intensity) in enumerate(peaks, 1):
            relative = intensity / base_intensity * 100 if base_intensity and base_intensity > 0 else np.nan
            neutral_loss = precursor - fragment_mz if precursor is not None and precursor >= fragment_mz else np.nan
            peak_rows.append(
                {
                    "Feature ID": feature_id,
                    "FD spectrum index": spectrum_index,
                    "FD database ID": clean(metadata.get("DATABASE_ID")),
                    "FD precursor type": clean(metadata.get("Precursor_type")),
                    "FD precursor m/z": precursor,
                    "Peak index": peak_index,
                    "Fragment m/z": fragment_mz,
                    "Fragment intensity": intensity,
                    "Relative intensity (%)": relative,
                    "Neutral loss from precursor (Da)": neutral_loss,
                    "Base peak": peak_index == int(np.argmax(positive_intensities)) + 1 if len(positive_intensities) else False,
                }
            )
        metadata = {}
        peaks = []

    with path.open("r", encoding=encoding, errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                emit()
                continue
            match = re.match(r"^([^:]+):\s*(.*)$", line)
            if match and not re.match(r"^[+-]?\d", line):
                metadata[match.group(1).strip()] = match.group(2).strip()
                continue
            parts = re.split(r"\s+", line)
            if len(parts) >= 2:
                try:
                    peaks.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
            if progress and len(records) and len(records) % 300 == 0:
                progress(min(68, 59 + len(records) // 150), f"Scanning fragmentation spectra… {len(records):,} records")
        emit()
    if not records:
        raise ValueError("FD MSP contains no valid spectrum records.")
    return pd.DataFrame(records), pd.DataFrame(peak_rows)


def parse_msp(path: Path, progress: Callable[[int, str], None] | None = None) -> pd.DataFrame:
    summary, _peaks = parse_msp_detailed(path, progress)
    return summary


def _assert_columns(role: str, frame: pd.DataFrame) -> None:
    if role in {"II", "FD"}:
        return
    missing = [field for field in REQUIRED_FIELDS[role] if field not in frame.columns]
    if missing:
        raise ValueError(f"{role} is missing required field{'s' if len(missing) != 1 else ''}: {', '.join(missing)}")


def _series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column in frame.columns:
        return frame[column].map(clean)
    return pd.Series([""] * len(frame), index=frame.index, dtype=str)


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace("", np.nan), errors="coerce")


def _sample_cv(row: pd.Series) -> float:
    values = pd.to_numeric(row, errors="coerce").dropna()
    values = values[values > 0]
    if len(values) < 2 or float(values.mean()) == 0:
        return np.nan
    return float(values.std(ddof=1) / abs(values.mean()) * 100)


def _positive_median(row: pd.Series) -> float:
    values = pd.to_numeric(row, errors="coerce").dropna()
    values = values[values > 0]
    return float(values.median()) if len(values) else np.nan


def _detection_rate(row: pd.Series) -> float:
    values = pd.to_numeric(row, errors="coerce").dropna()
    return float((values > 0).mean()) if len(values) else np.nan


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _selection_reason(best: pd.Series, second: pd.Series | None) -> str:
    if second is None:
        return "Unique accepted compound ID"
    comparisons = [
        ("_score_sort", "Highest overall Score"),
        ("_fragment_sort", "Score tie → highest Fragmentation Score"),
        ("_mass_sort", "Score + fragmentation tie → smallest absolute mass error"),
        ("_isotope_sort", "Score + fragmentation + mass tie → highest Isotope Similarity"),
        ("Evidence completeness", "Four-metric tie → most complete orthogonal evidence"),
        ("FD observed peaks", "Four-metric tie → richer linked fragmentation spectrum"),
        ("QC detection rate", "Four-metric tie → stronger QC detection"),
        ("QC CV%", "Four-metric tie → lower pooled-QC CV"),
    ]
    for field, label in comparisons:
        left = best.get(field)
        right = second.get(field)
        if pd.isna(left) and pd.isna(right):
            continue
        if left != right:
            return label
    return "Unresolved evidence tie → stable Feature ID order; manual review flagged"


def _selection_trace(group: pd.DataFrame) -> str:
    entries: list[str] = []
    for rank, (_, row) in enumerate(group.iterrows(), 1):
        ppm = row.get("Absolute Mass Error (ppm)")
        ppm_text = "NA" if pd.isna(ppm) else f"{float(ppm):.3f}"
        entries.append(
            f"{rank}. {row['Feature ID']} [Score {clean(row.get('Score')) or 'NA'}; "
            f"Frag {clean(row.get('Fragmentation Score')) or 'NA'}; |ppm| {ppm_text}; "
            f"Isotope {clean(row.get('Isotope Similarity')) or 'NA'}]"
        )
    return " | ".join(entries)


def _percentile(values: Iterable[float], probability: float) -> float:
    array = np.asarray([value for value in values if value is not None and math.isfinite(float(value))], dtype=float)
    return float(np.quantile(array, probability)) if len(array) else 0.0


def _numeric_summary(values: Iterable[float]) -> dict[str, float | int]:
    array = np.asarray([value for value in values if value is not None and math.isfinite(float(value))], dtype=float)
    if not len(array):
        return {"n": 0, "min": 0, "q10": 0, "q25": 0, "median": 0, "q75": 0, "q90": 0, "max": 0}
    q = np.quantile(array, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1])
    return dict(zip(("min", "q10", "q25", "median", "q75", "q90", "max"), map(float, q))) | {"n": len(array)}


def _normal_pdf(values: np.ndarray, mean: float, variance: float) -> np.ndarray:
    return np.exp(-((values - mean) ** 2) / (2 * variance)) / math.sqrt(2 * math.pi * variance)


def _fit_two_gaussian(values: np.ndarray) -> dict[str, Any] | None:
    values = np.asarray(values, dtype=float)
    if len(values) < 30:
        return None
    means = np.quantile(values, [0.3, 0.7]).astype(float)
    variance = max(float(np.var(values)), 1e-6)
    variances = np.array([variance, variance], dtype=float)
    weights = np.array([0.5, 0.5], dtype=float)
    for _ in range(120):
        p1 = weights[0] * _normal_pdf(values, means[0], variances[0])
        p2 = weights[1] * _normal_pdf(values, means[1], variances[1])
        total = np.maximum(p1 + p2, np.finfo(float).tiny)
        resp = p1 / total
        n1 = float(resp.sum())
        n2 = len(values) - n1
        if n1 < 2 or n2 < 2:
            return None
        next_means = np.array([(resp * values).sum() / n1, ((1 - resp) * values).sum() / n2])
        next_variances = np.array([
            max(float((resp * (values - next_means[0]) ** 2).sum() / n1), 1e-6),
            max(float(((1 - resp) * (values - next_means[1]) ** 2).sum() / n2), 1e-6),
        ])
        delta = float(np.max(np.abs(next_means - means)))
        means, variances, weights = next_means, next_variances, np.array([n1 / len(values), n2 / len(values)])
        if delta < 1e-7:
            break
    if means[0] > means[1]:
        means = means[::-1]
        variances = variances[::-1]
        weights = weights[::-1]
    mixture = np.maximum(weights[0] * _normal_pdf(values, means[0], variances[0]) + weights[1] * _normal_pdf(values, means[1], variances[1]), np.finfo(float).tiny)
    log_likelihood = float(np.log(mixture).sum())
    bic = 5 * math.log(len(values)) - 2 * log_likelihood
    grid = np.linspace(means[0], means[1], 2001)[1:-1]
    density = weights[0] * _normal_pdf(grid, means[0], variances[0]) + weights[1] * _normal_pdf(grid, means[1], variances[1])
    valley = float(grid[int(np.argmin(density))])
    ashman_d = float(abs(means[1] - means[0]) / math.sqrt(float((variances[0] + variances[1]) / 2)))
    return {"bic": bic, "means": means.tolist(), "variances": variances.tolist(), "weights": weights.tolist(), "intersection": valley, "ashman_d": ashman_d}


def _one_gaussian_bic(values: np.ndarray) -> float:
    mean = float(np.mean(values))
    variance = max(float(np.var(values)), 1e-6)
    likelihood = np.maximum(_normal_pdf(values, mean, variance), np.finfo(float).tiny)
    return 2 * math.log(len(values)) - 2 * float(np.log(likelihood).sum())


def _two_means_boundary(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if not len(values):
        return 0.0
    left, right = map(float, np.quantile(values, [0.25, 0.75]))
    for _ in range(100):
        mask = np.abs(values - left) <= np.abs(values - right)
        if mask.all() or (~mask).all():
            return float(np.median(values))
        next_left = float(values[mask].mean())
        next_right = float(values[~mask].mean())
        delta = max(abs(next_left - left), abs(next_right - right))
        left, right = next_left, next_right
        if delta < 1e-9:
            break
    return (left + right) / 2


def _mulberry32(seed: int) -> Callable[[], float]:
    state = seed & 0xFFFFFFFF

    def random_value() -> float:
        nonlocal state
        state = (state + 0x6D2B79F5) & 0xFFFFFFFF
        value = state
        value = ((value ^ (value >> 15)) * (value | 1)) & 0xFFFFFFFF
        value ^= (value + (((value ^ (value >> 7)) * (value | 61)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        return ((value ^ (value >> 14)) & 0xFFFFFFFF) / 4294967296

    return random_value


def _bootstrap(values: np.ndarray, estimator: Callable[[np.ndarray], float], seed: int) -> tuple[float, float]:
    if len(values) < 5:
        estimate = estimator(values)
        return estimate, estimate
    random_value = _mulberry32(seed)
    estimates = np.empty(500, dtype=float)
    for index in range(500):
        sample = np.asarray([values[int(random_value() * len(values))] for _ in range(len(values))], dtype=float)
        estimates[index] = estimator(sample)
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def _passes(frame: pd.DataFrame, settings: dict[str, Any]) -> pd.Series:
    score = pd.to_numeric(frame["Score"], errors="coerce")
    fragment = pd.to_numeric(frame["Fragmentation Score"], errors="coerce")
    mass = pd.to_numeric(frame["Mass Error (ppm)"], errors="coerce").abs()
    isotope = pd.to_numeric(frame["Isotope Similarity"], errors="coerce")
    fragment_pass = fragment.ge(settings["fragmentation"])
    if settings.get("allow_zero_fragmentation", False):
        fragment_pass = fragment_pass | fragment.eq(0)
    return score.ge(settings["score"]) & fragment_pass & mass.le(settings["abs_mass_error"]) & isotope.ge(settings["isotope"])


QC_FILTER_KEYS = {
    "min_qc_detection",
    "max_qc_cv",
    "min_biological_detection",
    "min_biological_blank_ratio",
    "max_abs_qc_drift",
    "max_d_ratio",
    "require_ms2",
    "min_fragment_peaks",
    "apply_blank_filter",
    "apply_qc_filter",
    "apply_drift_filter",
    "apply_d_ratio_filter",
}


def evaluate_filters(selected: pd.DataFrame, settings: dict[str, Any]) -> pd.DataFrame:
    output = selected.copy()
    score = pd.to_numeric(output["Score"], errors="coerce")
    fragment = pd.to_numeric(output["Fragmentation Score"], errors="coerce")
    mass = pd.to_numeric(output["Mass Error (ppm)"], errors="coerce").abs()
    isotope = pd.to_numeric(output["Isotope Similarity"], errors="coerce")
    output["Pass Score cutoff"] = score.ge(settings["score"])
    output["Pass Fragmentation cutoff"] = fragment.ge(settings["fragmentation"])
    if settings.get("allow_zero_fragmentation", False):
        output["Pass Fragmentation cutoff"] = output["Pass Fragmentation cutoff"] | fragment.eq(0)
    output["Pass absolute mass-error cutoff"] = mass.le(settings["abs_mass_error"])
    output["Pass Isotope Similarity cutoff"] = isotope.ge(settings["isotope"])
    identity_columns = [
        "Pass Score cutoff",
        "Pass Fragmentation cutoff",
        "Pass absolute mass-error cutoff",
        "Pass Isotope Similarity cutoff",
    ]
    output["Identification evidence pass"] = output[identity_columns].all(axis=1)
    identification_reasons: list[str] = []
    for _, row in output.iterrows():
        failed = [column.removeprefix("Pass ") for column in identity_columns if not bool(row[column])]
        identification_reasons.append("; ".join(failed) if failed else "No identification-evidence rule failed")
    output["Identification filter fail reasons"] = identification_reasons
    if QC_FILTER_KEYS.intersection(settings):
        output = evaluate_analytical_qc(output, settings)
    else:
        output["Analytical QC decision"] = "Not applied"
        output["Analytical QC fail reasons"] = "Analytical-QC rules were not enabled for this filter operation"
        output["Analytical QC evidence available"] = "Not applied"
        output["Analytical QC pass"] = True
    output["Filter pass"] = output["Identification evidence pass"] & output["Analytical QC pass"]
    output["Filter decision"] = np.where(output["Filter pass"], "Pass", "Fail")
    output["Filter fail reasons"] = [
        "; ".join(
            part
            for part in (
                "" if row["Identification evidence pass"] else str(row["Identification filter fail reasons"]),
                "" if row["Analytical QC pass"] else str(row["Analytical QC fail reasons"]),
            )
            if part
        )
        or "None"
        for _, row in output.iterrows()
    ]
    return output


def recommend_thresholds(selected: pd.DataFrame, mass_tolerance: float) -> dict[str, Any]:
    scores = pd.to_numeric(selected["Score"], errors="coerce").dropna().to_numpy(float)
    fragments = pd.to_numeric(selected["Fragmentation Score"], errors="coerce").dropna().to_numpy(float)
    positive_fragments = fragments[fragments > 0]
    abs_mass = pd.to_numeric(selected["Mass Error (ppm)"], errors="coerce").abs().dropna().to_numpy(float)
    isotopes = pd.to_numeric(selected["Isotope Similarity"], errors="coerce").dropna().to_numpy(float)
    mass_gmm = _fit_two_gaussian(abs_mass)
    delta_bic = _one_gaussian_bic(abs_mass) - mass_gmm["bic"] if mass_gmm else 0.0
    mixture_used = bool(mass_gmm and delta_bic > 10 and mass_gmm["ashman_d"] >= 2)
    mass_boundary = _two_means_boundary(abs_mass)
    mass_value = min(mass_tolerance, mass_boundary if mixture_used else min(5.0, _percentile(abs_mass, 0.75)))

    score_ci = _bootstrap(scores, lambda array: _percentile(array, 0.25), 24011)
    fragment_ci = _bootstrap(positive_fragments, lambda array: _percentile(array, 0.25), 24012)
    isotope_ci = _bootstrap(isotopes, lambda array: _percentile(array, 0.25), 24013)
    mass_ci = _bootstrap(abs_mass, _two_means_boundary if mixture_used else lambda array: _percentile(array, 0.75), 24014)

    balanced = {
        "score": round(_percentile(scores, 0.25), 1),
        "fragmentation": round(_percentile(positive_fragments, 0.25), 1),
        "abs_mass_error": round(mass_value, 1),
        "isotope": round(_percentile(isotopes, 0.25), 1),
        "allow_zero_fragmentation": True,
    }
    inclusive = {
        "score": round(_percentile(scores, 0.10), 1),
        "fragmentation": round(_percentile(positive_fragments, 0.10), 1),
        "abs_mass_error": round(min(mass_tolerance, _percentile(abs_mass, 0.90)), 1),
        "isotope": round(_percentile(isotopes, 0.10), 1),
        "allow_zero_fragmentation": True,
    }
    score_gmm = _fit_two_gaussian(scores)
    fragment_gmm = _fit_two_gaussian(positive_fragments)
    isotope_gmm = _fit_two_gaussian(isotopes)
    stringent = {
        "score": round(score_gmm["intersection"] if score_gmm else _percentile(scores, 0.75), 1),
        "fragmentation": round(fragment_gmm["intersection"] if fragment_gmm else _percentile(positive_fragments, 0.50), 1),
        "abs_mass_error": round(mass_value, 1),
        "isotope": round(isotope_gmm["intersection"] if isotope_gmm else _percentile(isotopes, 0.50), 1),
        "allow_zero_fragmentation": False,
    }
    preset_descriptions = {
        "Inclusive": "Discovery-oriented: 10th-percentile evidence bounds; zero fragmentation retained as a review state.",
        "Balanced": "Default shortlist: entity-level lower-quartile evidence plus a validated mass-error boundary.",
        "Stringent": "High-evidence triage: distribution separation and mandatory positive fragmentation support.",
    }
    presets: list[dict[str, Any]] = []
    for name, settings in (("Inclusive", inclusive), ("Balanced", balanced), ("Stringent", stringent)):
        presets.append({"name": name, **settings, "estimated_passes": int(_passes(selected, settings).sum()), "description": preset_descriptions[name]})
    return {
        "score": {
            "value": balanced["score"], "ci_low": round(score_ci[0], 1), "ci_high": round(score_ci[1], 1),
            "method": "Selected-entity Q1",
            "rationale": "Robust lower quartile after identity collapse; avoids assuming a universal Progenesis Score cutoff when RT/CCS terms are unused.",
        },
        "fragmentation": {
            "value": balanced["fragmentation"], "ci_low": round(fragment_ci[0], 1), "ci_high": round(fragment_ci[1], 1),
            "method": "Q1 among positive scores",
            "rationale": "Zero is retained explicitly because it may encode absent, disabled, or unmatched fragmentation evidence; strict mode requires positive support.",
        },
        "abs_mass_error": {
            "value": balanced["abs_mass_error"], "ci_low": round(min(mass_ci[0], mass_tolerance), 1), "ci_high": round(min(mass_ci[1], mass_tolerance), 1),
            "method": "Stable two-cluster boundary, Gaussian separation-validated" if mixture_used else "Robust Q3, capped at 5 ppm and search tolerance",
            "rationale": "Absolute ppm is split by a deterministic minimum-variance boundary only when ΔBIC > 10 and Ashman D ≥ 2; otherwise a conservative robust bound is used.",
        },
        "isotope": {
            "value": balanced["isotope"], "ci_low": round(isotope_ci[0], 1), "ci_high": round(isotope_ci[1], 1),
            "method": "Selected-entity Q1",
            "rationale": "Retains the central 75% of entity-level isotope evidence without claiming an instrument-independent universal cutoff.",
        },
        "presets": presets,
        "profiles": {
            "Score": _numeric_summary(scores),
            "Fragmentation Score": _numeric_summary(fragments),
            "Absolute Mass Error (ppm)": _numeric_summary(abs_mass),
            "Isotope Similarity": _numeric_summary(isotopes),
        },
        "diagnostics": {
            "mass_mixture_used": mixture_used,
            "mass_delta_bic": delta_bic,
            "mass_ashman_d": mass_gmm["ashman_d"] if mass_gmm else 0.0,
            "mass_intersection": mass_boundary if mixture_used else None,
            "zero_fragmentation_share": float((fragments == 0).mean()) if len(fragments) else 0.0,
        },
    }


def apply_filters(selected: pd.DataFrame, settings: dict[str, Any]) -> pd.DataFrame:
    decisions = evaluate_filters(selected, settings)
    return decisions.loc[decisions["Filter pass"]].copy().reset_index(drop=True)


def _make_audit(role: str, path: Path, frame: pd.DataFrame, key_field: str, cm_keys: set[str], note: str) -> dict[str, Any]:
    values = _series(frame, key_field)
    values = values[values.ne("")]
    source_keys = set(values)
    duplicates = int(len(values) - values.nunique())
    coverage = len(source_keys & cm_keys) / len(cm_keys) if cm_keys else 0.0
    status = "Pass" if coverage == 1 else "Review" if coverage >= 0.95 else "Fail"
    if duplicates and role in {"AM", "FD"} and status == "Pass":
        status = "Review"
    return {
        "Role": role,
        "File": path.name,
        "Size (MiB)": round(path.stat().st_size / 1024 / 1024, 3),
        "SHA-256": _sha256(path),
        "Records": len(frame),
        "Fields": len(frame.columns),
        "Key field": key_field,
        "Distinct feature keys": len(source_keys),
        "Duplicate rows": duplicates,
        "CM coverage (%)": round(coverage * 100, 3),
        "Status": status,
        "Finding": note,
    }


def _resolve_representatives(base: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int, int]:
    output = base.copy()
    prior_columns = {
        "Read count for entity",
        "Selection rank",
        "Selected representative",
        "Selection reason",
        "Selection trace",
        "Read type",
        "Primary-metric tie",
        "Manual identity review",
    }
    output = output.drop(columns=[column for column in prior_columns if column in output.columns])
    for field in ("Score", "Fragmentation Score", "Mass Error (ppm)", "Isotope Similarity"):
        output[field] = pd.to_numeric(output[field], errors="coerce")
    output["Absolute Mass Error (ppm)"] = output["Mass Error (ppm)"].abs()
    output["_score_sort"] = output["Score"].round(1).fillna(-np.inf)
    output["_fragment_sort"] = output["Fragmentation Score"].round(2).fillna(-np.inf)
    output["_mass_sort"] = output["Absolute Mass Error (ppm)"].round(3).fillna(np.inf)
    output["_isotope_sort"] = output["Isotope Similarity"].round(2).fillna(-np.inf)
    output["_qc_detection_sort"] = pd.to_numeric(output.get("QC detection rate"), errors="coerce").fillna(-1)
    output["_qc_cv_sort"] = pd.to_numeric(output.get("QC CV%"), errors="coerce").fillna(np.inf)
    output["_bio_detection_sort"] = pd.to_numeric(output.get("Biological detection rate"), errors="coerce").fillna(-1)
    output["_bio_abundance_sort"] = pd.to_numeric(output.get("Biological median abundance"), errors="coerce").fillna(-1)
    output = output.sort_values(
        [
            "Entity key",
            "_score_sort",
            "_fragment_sort",
            "_mass_sort",
            "_isotope_sort",
            "Evidence completeness",
            "FD observed peaks",
            "_qc_detection_sort",
            "_qc_cv_sort",
            "_bio_detection_sort",
            "_bio_abundance_sort",
            "Feature ID",
        ],
        ascending=[True, False, False, True, False, False, False, False, True, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    output["Read count for entity"] = output.groupby("Entity key")["Entity key"].transform("size")
    output["Selection rank"] = output.groupby("Entity key").cumcount() + 1
    output["Selected representative"] = output["Selection rank"].eq(1)
    reasons: dict[str, str] = {}
    traces: dict[str, str] = {}
    primary_ties: dict[str, bool] = {}
    manual_review: dict[str, bool] = {}
    primary_fields = ["_score_sort", "_fragment_sort", "_mass_sort", "_isotope_sort"]
    all_tie_fields = primary_fields + ["Evidence completeness", "FD observed peaks", "_qc_detection_sort", "_qc_cv_sort", "_bio_detection_sort", "_bio_abundance_sort"]
    for entity, group in output.groupby("Entity key", sort=False):
        second = group.iloc[1] if len(group) > 1 else None
        reasons[entity] = _selection_reason(group.iloc[0], second)
        traces[entity] = _selection_trace(group)
        primary_ties[entity] = bool(second is not None and all(group.iloc[0][field] == second[field] for field in primary_fields))
        manual_review[entity] = bool(second is not None and all(group.iloc[0][field] == second[field] for field in all_tie_fields))
    output["Selection reason"] = [
        reasons[entity] if rank == 1 else f"Not selected; evidence rank {rank}"
        for entity, rank in zip(output["Entity key"], output["Selection rank"])
    ]
    output["Selection trace"] = output["Entity key"].map(traces)
    output["Primary-metric tie"] = output["Entity key"].map(primary_ties)
    output["Manual identity review"] = output["Entity key"].map(manual_review)
    output["Read type"] = np.where(output["Read count for entity"].gt(1), "Repeated accepted ID", "Unique accepted ID")
    helper_columns = [column for column in output.columns if column.startswith("_")]
    accepted_reads = output.drop(columns=helper_columns).copy()
    duplicate_reads = accepted_reads.loc[accepted_reads["Read count for entity"].gt(1)].copy()
    selected_reads = accepted_reads.loc[accepted_reads["Selected representative"]].copy()
    selected_reads = selected_reads.sort_values(["Score", "Accepted Compound ID"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    duplicate_groups = int((accepted_reads.groupby("Entity key").size() > 1).sum())
    extra_reads_collapsed = len(accepted_reads) - len(selected_reads)
    return accepted_reads, duplicate_reads, selected_reads, duplicate_groups, extra_reads_collapsed


def analyze_files(
    files: dict[str, str | Path],
    parameters: dict[str, Any] | None = None,
    progress: Callable[[int, str], None] | None = None,
) -> AnalysisBundle:
    parameters = {
        "project_name": "LCMS Curation Project",
        "assay_type": "metabolomics",
        "biological_system": "human",
        "sample_matrix": "",
        "extraction_method": "unknown",
        "analyzed_phase": "unknown",
        "acquisition_mode": "unknown",
        "ion_mode": "negative",
        "mass_tolerance": 10.0,
        "qc_pattern": "QC|Pool",
        **(parameters or {}),
    }
    paths = {role: Path(files[role]).expanduser().resolve() for role in FILE_ROLES}
    missing = [role for role, path in paths.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing input file(s): {', '.join(missing)}")
    update = progress or (lambda _percent, _message: None)

    update(2, "Reading compound measurements…")
    cm, _cm_structure = read_three_header_csv(paths["CM"])
    _assert_columns("CM", cm)
    update(12, "Reading adduct measurements…")
    am, _ = read_three_header_csv(paths["AM"])
    _assert_columns("AM", am)
    update(22, "Reading every identification candidate…")
    ci = read_standard_csv(paths["CI"], {"Compound", "Compound ID", "Accepted?"})
    _assert_columns("CI", ci)
    update(38, "Reading accepted-compound properties…")
    acp = read_standard_csv(paths["ACP"], {"Compound", "Compound ID"})
    _assert_columns("ACP", acp)
    update(45, "Scanning isotope XML…")
    ii = parse_isotope_xml(paths["II"], update)
    update(59, "Scanning every fragmentation spectrum…")
    fd, fd_peaks = parse_msp_detailed(paths["FD"], update)
    update(69, "Reconciling feature identifiers across all six files…")

    accepted_ci = ci.loc[_series(ci, "Accepted?").ne("")].copy()
    accepted_first = accepted_ci.drop_duplicates("Compound", keep="first").set_index("Compound", drop=False)
    acp_first = acp.drop_duplicates("Compound", keep="first").set_index("Compound", drop=False)

    base = pd.DataFrame(
        {
            "Feature ID": _series(cm, "Metadata|Compound"),
            "Accepted Compound ID": _series(cm, "Tags|Accepted Compound ID"),
            "Accepted Description": _series(cm, "Tags|Accepted Description"),
            "Formula": _series(cm, "Tags|Formula"),
            "Adducts": _series(cm, "Tags|Adducts"),
            "Compound Link": _series(cm, "Tags|Compound Link"),
            "Score": _series(cm, "Tags|Score"),
            "Fragmentation Score": _series(cm, "Tags|Fragmentation Score"),
            "Mass Error (ppm)": _series(cm, "Tags|Mass Error (ppm)"),
            "Isotope Similarity": _series(cm, "Tags|Isotope Similarity"),
            "Retention Time Error (min)": _series(cm, "Tags|Retention Time Error (mins)"),
            "m/z": _series(cm, "Metadata|m/z"),
            "Charge": _series(cm, "Metadata|Charge"),
            "Neutral mass (Da)": _series(cm, "Metadata|Neutral mass (Da)"),
            "Retention time (min)": _series(cm, "Metadata|Retention time (min)"),
            "Peak width (min)": _series(cm, "Metadata|Chromatographic peak width (min)"),
            "Identifications": _series(cm, "Metadata|Identifications"),
            "ANOVA p": _series(cm, "Metadata|Anova (p)"),
            "q value": _series(cm, "Metadata|q Value"),
            "Max fold change": _series(cm, "Metadata|Max Fold Change"),
            "Highest mean group": _series(cm, "Metadata|Highest Mean"),
            "Lowest mean group": _series(cm, "Metadata|Lowest Mean"),
            "Isotope distribution": _series(cm, "Metadata|Isotope Distribution"),
            "Maximum abundance": _series(cm, "Metadata|Maximum Abundance"),
            "Minimum CV%": _series(cm, "Metadata|Minimum CV%"),
        }
    )
    base["Canonical name"] = base["Accepted Description"].map(canonical_name)
    placeholder_link_count = int(base["Compound Link"].map(lambda value: bool(clean(value)) and not bool(clean_link(value))).sum())
    base["Compound Link"] = base["Compound Link"].map(clean_link)

    for index, feature_id in enumerate(base["Feature ID"]):
        accepted = accepted_first.loc[feature_id] if feature_id in accepted_first.index else None
        accepted_id = base.at[index, "Accepted Compound ID"] or (clean(accepted.get("Compound ID")) if accepted is not None else "")
        description = base.at[index, "Accepted Description"] or (clean(accepted.get("Description")) if accepted is not None else "")
        formula = base.at[index, "Formula"] or (clean(accepted.get("Formula")) if accepted is not None else "")
        adducts = base.at[index, "Adducts"] or (clean(accepted.get("Adducts")) if accepted is not None else "")
        base.at[index, "Accepted Compound ID"] = accepted_id
        base.at[index, "Accepted Description"] = description
        base.at[index, "Canonical name"] = canonical_name(description)
        base.at[index, "Formula"] = formula
        base.at[index, "Adducts"] = adducts
        if not base.at[index, "Compound Link"] and accepted is not None:
            base.at[index, "Compound Link"] = clean_link(accepted.get("Link"))
        if accepted is not None:
            for target, source in (
                ("Score", "Score"),
                ("Fragmentation Score", "Fragmentation Score"),
                ("Mass Error (ppm)", "Mass Error (ppm)"),
                ("Isotope Similarity", "Isotope Similarity"),
            ):
                if not clean(base.at[index, target]):
                    base.at[index, target] = clean(accepted.get(source))
        if not accepted_id and feature_id in acp_first.index:
            base.at[index, "Accepted Compound ID"] = clean(acp_first.loc[feature_id].get("Compound ID"))

    base["Entity key"] = base["Accepted Compound ID"]
    fallback = "FORMULA_NAME::" + base["Formula"] + "::" + base["Canonical name"].str.casefold()
    empty_entity = base["Entity key"].eq("")
    base.loc[empty_entity, "Entity key"] = fallback.loc[empty_entity]
    no_annotation = empty_entity & base["Formula"].eq("") & base["Canonical name"].eq("")
    base.loc[no_annotation, "Entity key"] = "FEATURE::" + base.loc[no_annotation, "Feature ID"]

    ci_counts = ci.groupby("Compound", dropna=False).size().rename("CI candidate count")
    accepted_counts = accepted_ci.groupby("Compound", dropna=False).size().rename("CI accepted count")
    am_summary = am.groupby(_series(am, "Metadata|Compound"), sort=False).agg(
        **{
            "AM row count": ("Metadata|Compound", "size"),
            "AM ions": ("Metadata|Ion", safe_join),
            "AM adduct types": ("Metadata|Adduct type", safe_join),
            "AM charges": ("Metadata|Charge", safe_join),
        }
    )
    base = base.join(ci_counts, on="Feature ID").join(accepted_counts, on="Feature ID").join(am_summary, on="Feature ID")
    base["CI candidate count"] = base["CI candidate count"].fillna(0).astype(int)
    base["CI accepted count"] = base["CI accepted count"].fillna(0).astype(int)
    base["AM row count"] = base["AM row count"].fillna(0).astype(int)

    acp_small = acp.copy()
    acp_small["Feature ID"] = _series(acp_small, "Compound")
    acp_small = acp_small.drop_duplicates("Feature ID").rename(
        columns={
            "Retention time (min)": "ACP retention time",
            "CCS (angstrom^2)": "ACP CCS (angstrom^2)",
            "Adduct": "ACP adduct",
        }
    )
    acp_columns = [column for column in ("Feature ID", "ACP retention time", "ACP CCS (angstrom^2)", "ACP adduct") if column in acp_small.columns]
    base = base.merge(acp_small[acp_columns], on="Feature ID", how="left")
    base = base.merge(ii, on="Feature ID", how="left")

    fd_numeric = fd.copy()
    fd_numeric["FD observed peaks"] = pd.to_numeric(fd_numeric["FD observed peaks"], errors="coerce").fillna(0)
    fd_summary = fd_numeric.groupby("Feature ID", sort=False).agg(
        **{
            "FD spectrum count": ("Feature ID", "size"),
            "FD precursor types": ("FD precursor type", safe_join),
            "FD precursor m/z values": ("FD precursor m/z", safe_join),
            "FD observed peaks": ("FD observed peaks", "sum"),
            "FD peak-count integrity": ("FD peak-count integrity", "all"),
            "FD base peaks m/z": ("FD base peak m/z", safe_join),
            "FD database IDs": ("FD database ID", safe_join),
            "FD median spectral entropy": ("FD spectral entropy", "median"),
            "FD median normalized spectral entropy": ("FD normalized spectral entropy", "median"),
            "FD median effective peak count": ("FD effective peak count", "median"),
            "FD minimum top-five intensity share (%)": ("FD top-five intensity share (%)", "min"),
            "FD spectrum quality flags": ("FD spectrum quality flag", safe_join),
        }
    )
    base = base.join(fd_summary, on="Feature ID")
    for column in ("FD spectrum count", "FD observed peaks"):
        base[column] = base[column].fillna(0).astype(int)
    base["FD peak-count integrity"] = (
        base["FD peak-count integrity"].map(lambda value: bool(value) if pd.notna(value) else False).astype(bool)
    )

    update(73, "MS² spectra and fragment peaks linked to compound features…")
    normalized_columns = [column for column in cm.columns if column.startswith("Normalised abundance|")]
    raw_columns = [column for column in cm.columns if column.startswith("Raw abundance|")]
    abundance = cm[["Metadata|Compound", *normalized_columns, *raw_columns]].rename(columns={"Metadata|Compound": "Feature ID"})
    base = base.merge(abundance, on="Feature ID", how="left")
    sample_metadata_value = parameters.get("sample_metadata")
    if isinstance(sample_metadata_value, pd.DataFrame):
        sample_metadata = normalize_sample_metadata(sample_metadata_value, normalized_columns)
    else:
        sample_metadata = build_sample_metadata(
            normalized_columns,
            raw_columns,
            str(parameters.get("qc_pattern") or "QC|Pool"),
            parameters,
        )
    update(76, "Computing biological, pooled-QC, blank, replicate, and drift metrics…")
    base = compute_abundance_metrics(base, sample_metadata)
    update(81, "Deriving formula chemistry, polarity, family, and adduct evidence…")
    base = enrich_chemical_properties(base, str(parameters.get("ion_mode") or "mixed"))

    presence_columns = [
        "Accepted Compound ID", "Accepted Description", "Formula", "Adducts", "Score", "Fragmentation Score",
        "Mass Error (ppm)", "Isotope Similarity", "Retention Time Error (min)", "ACP CCS (angstrom^2)",
    ]
    presence = pd.DataFrame({column: base.get(column, "") for column in presence_columns})
    present = presence.apply(lambda column: column.map(lambda value: clean(value) != ""))
    base["Evidence completeness"] = present.sum(axis=1) + base["FD spectrum count"].gt(0).astype(int)

    update(84, "Reconciling accepted identities and resolving repeated chromatographic reads…")
    accepted_reads, duplicate_reads, selected_reads, duplicate_groups, extra_reads_collapsed = _resolve_representatives(base)
    update(86, "Evaluating one representative per entity against analytical-QC evidence…")
    analytical_qc = evaluate_analytical_qc(selected_reads)

    cm_keys = set(_series(cm, "Metadata|Compound")) - {""}
    notes = {
        "CM": f"{len(cm):,} feature rows; {len(normalized_columns)} normalised and {len(raw_columns)} raw sample columns.",
        "AM": f"{len(am) - _series(am, 'Metadata|Compound').nunique():,} extra ion/adduct rows retained as one-to-many evidence.",
        "CI": f"{len(ci):,} candidate rows; {len(accepted_ci):,} user-accepted rows.",
        "II": f"{int(ii['II isotope-node count'].sum()):,} isotope nodes scanned.",
        "FD": f"{int(fd['FD observed peaks'].sum()):,} fragment peaks; {int((~fd['FD peak-count integrity']).sum()):,} declared-count mismatches.",
        "ACP": f"{int(_series(acp, 'CCS (angstrom^2)').eq('').sum()):,} rows lack CCS; {int(_series(acp, 'Adduct').eq('').sum()):,} rows lack ACP adduct.",
    }
    audits = pd.DataFrame(
        [
            _make_audit("CM", paths["CM"], cm.rename(columns={"Metadata|Compound": "Feature ID"}), "Feature ID", cm_keys, notes["CM"]),
            _make_audit("AM", paths["AM"], am.rename(columns={"Metadata|Compound": "Feature ID"}), "Feature ID", cm_keys, notes["AM"]),
            _make_audit("CI", paths["CI"], ci.rename(columns={"Compound": "Feature ID"}), "Feature ID", cm_keys, notes["CI"]),
            _make_audit("II", paths["II"], ii, "Feature ID", cm_keys, notes["II"]),
            _make_audit("FD", paths["FD"], fd, "Feature ID", cm_keys, notes["FD"]),
            _make_audit("ACP", paths["ACP"], acp.rename(columns={"Compound": "Feature ID"}), "Feature ID", cm_keys, notes["ACP"]),
        ]
    )

    warnings: list[str] = []
    failed_roles = audits.loc[audits["Status"].eq("Fail"), "Role"].tolist()
    if failed_roles:
        warnings.append(f"Cross-file coverage failed for: {', '.join(failed_roles)}. Results remain inspectable, but missing source records must be resolved before the run is treated as complete.")
    incomplete_fd = int((fd["Feature ID"].eq("") | fd["FD declared peaks"].isna()).sum())
    if incomplete_fd:
        warnings.append(f"FD contains {incomplete_fd:,} incomplete terminal/metadata record indicators; the MSP file may be truncated or malformed.")
    if _series(cm, "Tags|Retention Time Error (mins)").eq("").all():
        warnings.append("Retention-time error is absent for every accepted identification; RT cannot act as orthogonal identity evidence.")
    if _series(acp, "CCS (angstrom^2)").eq("").all():
        warnings.append("CCS is absent for every accepted compound; ion-mobility confirmation cannot be assessed.")
    if _series(acp, "Adduct").eq("").all():
        warnings.append("ACP adduct is empty; adduct evidence is taken from CM/CI/AM and checked for consistency.")
    if placeholder_link_count:
        warnings.append(
            f"{placeholder_link_count:,} vendor compound-link placeholders ending in NULL/NA were blanked; "
            "they are missing links, not database records."
        )
    if duplicate_groups:
        warnings.append(f"{duplicate_groups:,} accepted IDs map to multiple chromatographic features. One representative is selected, but all reads remain preserved because isomers, in-source products, adduct behavior, or false matches remain possible.")
    mapped_qc = int((sample_metadata["Sample role"] == "Pooled QC").sum())
    mapped_blanks = int(sample_metadata["Sample role"].isin({"Process blank", "Extraction blank", "Solvent blank"}).sum())
    if mapped_qc < 3:
        warnings.append(f"Only {mapped_qc} pooled-QC injections are mapped. QC CV and drift estimates require replicate pooled-QC measurements and should be interpreted cautiously.")
    if mapped_blanks == 0:
        warnings.append("No blank samples are mapped; blank contribution, carryover, and extraction-background filtering are not evaluable.")
    warnings.append("No target-decoy search, authentic-standard run, or externally validated labels were supplied; adaptive thresholds are triage rules, not an estimated identification FDR.")

    update(88, "Estimating robust dataset-adaptive thresholds…")
    thresholds = recommend_thresholds(selected_reads, float(parameters["mass_tolerance"]))
    update(100, "Analysis ready")
    return AnalysisBundle(
        generated_at=datetime.now(UTC).isoformat(),
        parameters=parameters,
        files=paths,
        audits=audits,
        warnings=warnings,
        cm=cm,
        am=am,
        ci=ci,
        accepted_ci=accepted_ci,
        acp=acp,
        ii=ii,
        fd=fd,
        fd_peaks=fd_peaks,
        sample_metadata=sample_metadata,
        accepted_reads=accepted_reads,
        duplicate_reads=duplicate_reads,
        selected_reads=selected_reads,
        analytical_qc=analytical_qc,
        duplicate_groups=duplicate_groups,
        extra_reads_collapsed=extra_reads_collapsed,
        thresholds=thresholds,
        normalized_columns=normalized_columns,
        raw_columns=raw_columns,
    )


def rebuild_with_sample_metadata(
    bundle: AnalysisBundle,
    metadata: pd.DataFrame,
    progress: Callable[[int, str], None] | None = None,
) -> AnalysisBundle:
    update = progress or (lambda _percent, _message: None)
    update(5, "Validating and normalizing the sample map…")
    normalized = normalize_sample_metadata(metadata, bundle.normalized_columns)
    update(18, "Recomputing pooled-QC, biological, blank, drift, and replicate metrics…")
    base = compute_abundance_metrics(bundle.accepted_reads, normalized)
    update(52, "Re-resolving one representative per accepted identity…")
    accepted, duplicates, selected, duplicate_groups, collapsed = _resolve_representatives(base)
    bundle.sample_metadata = normalized
    bundle.accepted_reads = accepted
    bundle.duplicate_reads = duplicates
    bundle.selected_reads = selected
    bundle.duplicate_groups = duplicate_groups
    bundle.extra_reads_collapsed = collapsed
    update(67, "Re-evaluating analytical-QC decisions…")
    bundle.analytical_qc = evaluate_analytical_qc(selected)
    update(78, "Re-estimating dataset-adaptive thresholds and bootstrap intervals…")
    bundle.thresholds = recommend_thresholds(selected, float(bundle.parameters["mass_tolerance"]))
    bundle.warnings = [
        warning
        for warning in bundle.warnings
        if not warning.startswith("Only ")
        and not warning.startswith("No pooled-QC")
        and not warning.startswith("No blank samples")
    ]
    mapped_qc = int((normalized["Sample role"] == "Pooled QC").sum())
    mapped_blanks = int(normalized["Sample role"].isin({"Process blank", "Extraction blank", "Solvent blank"}).sum())
    if mapped_qc < 3:
        bundle.warnings.append(f"Only {mapped_qc} pooled-QC injections are mapped. QC CV and drift estimates require replicate pooled-QC measurements and should be interpreted cautiously.")
    if mapped_blanks == 0:
        bundle.warnings.append("No blank samples are mapped; blank contribution, carryover, and extraction-background filtering are not evaluable.")
    update(100, "Sample map applied and downstream evidence rebuilt")
    return bundle
