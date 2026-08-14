from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd


SAMPLE_ROLES = (
    "Biological",
    "Pooled QC",
    "Process blank",
    "Extraction blank",
    "Solvent blank",
    "Reference material",
    "System suitability",
    "Calibration",
    "Excluded",
    "Unknown",
)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return re.sub(r"\s+", " ", str(value)).strip()


def split_abundance_column(column: str) -> tuple[str, str, str]:
    parts = str(column).split("|", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    return parts[0] if parts else "", "", parts[-1] if parts else ""


def infer_sample_role(sample_name: str, condition: str, qc_pattern: str = "QC|Pool") -> tuple[str, str]:
    text = f"{condition} {sample_name}".casefold()
    rules = (
        (r"solvent[ _-]*blank|blank[ _-]*solvent", "Solvent blank", "solvent-blank name/condition"),
        (r"extract(?:ion)?[ _-]*blank|blank[ _-]*extract", "Extraction blank", "extraction-blank name/condition"),
        (r"process[ _-]*blank|blank[ _-]*process", "Process blank", "process-blank name/condition"),
        (r"\bblank\b|\bblk\b", "Process blank", "generic blank name/condition"),
        (r"system[ _-]*suit|\bsst\b", "System suitability", "system-suitability cue"),
        (r"reference|\bsrm\b|nist", "Reference material", "reference-material cue"),
        (r"calibr|cal[ _-]?\d", "Calibration", "calibration cue"),
    )
    for pattern, role, reason in rules:
        if re.search(pattern, text, re.IGNORECASE):
            return role, reason
    try:
        if re.search(qc_pattern or "QC|Pool", f"{condition} {sample_name}", re.IGNORECASE):
            return "Pooled QC", f"matched QC regex: {qc_pattern or 'QC|Pool'}"
    except re.error:
        if re.search(r"\bqc\b|pool", text, re.IGNORECASE):
            return "Pooled QC", "matched fallback QC rule (invalid user regex)"
    return "Biological", "default study sample"


def build_sample_metadata(
    normalized_columns: list[str],
    raw_columns: list[str],
    qc_pattern: str = "QC|Pool",
    defaults: dict[str, Any] | None = None,
) -> pd.DataFrame:
    defaults = defaults or {}
    raw_by_sample: dict[str, str] = {}
    for column in raw_columns:
        _kind, _condition, sample = split_abundance_column(column)
        raw_by_sample[sample] = column
    rows: list[dict[str, Any]] = []
    for order, column in enumerate(normalized_columns, 1):
        _kind, condition, sample = split_abundance_column(column)
        role, rule = infer_sample_role(sample, condition, qc_pattern)
        subject = re.sub(r"(?:[-_]R|[-_]rep(?:licate)?\d*)$", "", sample, flags=re.IGNORECASE)
        technical = subject != sample
        rows.append(
            {
                "Sample name": sample,
                "Condition / group": condition,
                "Sample role": role,
                "Include": role != "Excluded",
                "Subject / biological unit": subject,
                "Technical replicate": technical,
                "Batch": _clean(defaults.get("batch")) or "1",
                "Injection order": order,
                "Dilution factor": 1.0,
                "Sample matrix": _clean(defaults.get("sample_matrix")),
                "Biological system": _clean(defaults.get("biological_system")),
                "Extraction method": _clean(defaults.get("extraction_method")),
                "Analyzed phase": _clean(defaults.get("analyzed_phase")),
                "Normalized abundance column": column,
                "Raw abundance column": raw_by_sample.get(sample, ""),
                "Automatic role rule": rule,
                "Notes": "",
            }
        )
    return pd.DataFrame(rows)


def normalize_sample_metadata(metadata: pd.DataFrame, normalized_columns: list[str]) -> pd.DataFrame:
    required = build_sample_metadata(normalized_columns, [], qc_pattern=r"a^")
    output = metadata.copy()
    if "Normalized abundance column" not in output.columns and "Sample name" in output.columns:
        name_to_column = {split_abundance_column(column)[2]: column for column in normalized_columns}
        output["Normalized abundance column"] = output["Sample name"].map(name_to_column).fillna("")
    for column in required.columns:
        if column not in output.columns:
            output[column] = required[column] if len(required) == len(output) else ""
    output = output.loc[output["Normalized abundance column"].isin(normalized_columns)].copy()
    output = output.drop_duplicates("Normalized abundance column", keep="last")
    missing = [column for column in normalized_columns if column not in set(output["Normalized abundance column"])]
    if missing:
        output = pd.concat([output, build_sample_metadata(missing, [], qc_pattern=r"a^")], ignore_index=True)
    role_map = {role.casefold(): role for role in SAMPLE_ROLES}
    output["Sample role"] = output["Sample role"].map(lambda value: role_map.get(_clean(value).casefold(), "Unknown"))
    output["Include"] = output["Include"].map(
        lambda value: value if isinstance(value, bool) else _clean(value).casefold() not in {"", "0", "false", "no", "exclude", "excluded"}
    )
    output["Injection order"] = pd.to_numeric(output["Injection order"], errors="coerce")
    fallback_order = {column: index for index, column in enumerate(normalized_columns, 1)}
    output["Injection order"] = [
        value if pd.notna(value) else fallback_order.get(column, index + 1)
        for index, (value, column) in enumerate(zip(output["Injection order"], output["Normalized abundance column"]))
    ]
    output["Dilution factor"] = pd.to_numeric(output["Dilution factor"], errors="coerce").fillna(1.0)
    output["Technical replicate"] = output["Technical replicate"].map(
        lambda value: value if isinstance(value, bool) else _clean(value).casefold() in {"1", "true", "yes", "y"}
    )
    order_map = {column: index for index, column in enumerate(normalized_columns)}
    output["_sort"] = output["Normalized abundance column"].map(order_map)
    return output.sort_values("_sort", kind="mergesort").drop(columns="_sort").reset_index(drop=True)


def merge_sample_metadata(template: pd.DataFrame, imported: pd.DataFrame) -> pd.DataFrame:
    incoming = imported.copy()
    aliases = {
        "sample": "Sample name",
        "sample_id": "Sample name",
        "role": "Sample role",
        "group": "Condition / group",
        "condition": "Condition / group",
        "subject": "Subject / biological unit",
        "batch_id": "Batch",
        "order": "Injection order",
        "injection_order": "Injection order",
        "dilution": "Dilution factor",
        "matrix": "Sample matrix",
        "phase": "Analyzed phase",
    }
    incoming = incoming.rename(columns={column: aliases.get(str(column).strip().casefold(), column) for column in incoming.columns})
    if "Sample name" not in incoming.columns:
        raise ValueError("Imported metadata requires a 'Sample name' or 'sample_id' column.")
    output = template.set_index("Sample name", drop=False).copy()
    incoming = incoming.drop_duplicates("Sample name", keep="last").set_index("Sample name", drop=False)
    for sample in output.index.intersection(incoming.index):
        for column in incoming.columns:
            if column in output.columns and _clean(incoming.at[sample, column]) != "":
                output.at[sample, column] = incoming.at[sample, column]
    return output.reset_index(drop=True)


def metadata_requirement_table() -> pd.DataFrame:
    rows = [
        ("Required", "Sample name", "Must match the leaf sample name in the CM three-row header."),
        ("Required", "Sample role", "Biological, pooled QC, blank, reference, calibration, or excluded."),
        ("Required", "Condition / group", "Biological comparison group or QC/blank condition."),
        ("Recommended", "Subject / biological unit", "Groups technical replicates without treating them as independent samples."),
        ("Recommended", "Batch + injection order", "Enables pooled-QC drift diagnostics and batch-aware review."),
        ("Recommended", "Dilution factor", "Enables dilution-series response checks when applicable."),
        ("Required for interpretation", "Sample matrix + biological system", "Controls context priors; does not prove origin."),
        ("Required for phase model", "Extraction method + analyzed phase", "Controls partition/retention plausibility."),
        ("Strongly recommended", "Process/extraction/solvent blanks", "Required to assess blank contribution and carryover."),
    ]
    return pd.DataFrame(rows, columns=["Priority", "Metadata field", "Why it matters"])


def project_metadata_frame(settings: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(
        [(key.replace("_", " ").title(), value) for key, value in settings.items() if not isinstance(value, (dict, list, pd.DataFrame))],
        columns=["Parameter", "Value"],
    )


def suggested_results_name(project_name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", _clean(project_name)).strip("_") or "LCMS_Curation"
    return f"{safe}_Results"


def default_project_name(cm_path: str | Path) -> str:
    path = Path(cm_path)
    return path.parent.name or path.stem
