from __future__ import annotations

import math
import re
from typing import Any

import numpy as np
import pandas as pd


MONOISOTOPIC_MASS = {
    "H": 1.00782503223,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "F": 18.99840316273,
    "Cl": 34.968852682,
    "Br": 78.9183376,
    "I": 126.904468,
    "Na": 22.989769282,
    "K": 38.9637064864,
    "Si": 27.97692653465,
}


def parse_formula(formula: Any) -> tuple[dict[str, int], str]:
    text = re.sub(r"\s+", "", str(formula or ""))
    text = re.sub(r"^[+-]?\d+", "", text)
    text = re.sub(r"[+-]\d*$", "", text)
    if not text:
        return {}, "Missing"
    if any(symbol in text for symbol in "().·[]"):
        return {}, "Unsupported grouped formula"
    tokens = list(re.finditer(r"([A-Z][a-z]?)(\d*)", text))
    if not tokens or "".join(match.group(0) for match in tokens) != text:
        return {}, "Unparsed"
    counts: dict[str, int] = {}
    for match in tokens:
        element = match.group(1)
        count = int(match.group(2) or 1)
        counts[element] = counts.get(element, 0) + count
    unknown = sorted(set(counts) - set(MONOISOTOPIC_MASS))
    return counts, f"Unknown element(s): {', '.join(unknown)}" if unknown else "Parsed"


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def formula_properties(formula: Any) -> dict[str, Any]:
    counts, status = parse_formula(formula)
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    n = counts.get("N", 0)
    o = counts.get("O", 0)
    p = counts.get("P", 0)
    s = counts.get("S", 0)
    halogens = sum(counts.get(element, 0) for element in ("F", "Cl", "Br", "I"))
    exact_mass = sum(MONOISOTOPIC_MASS.get(element, 0.0) * count for element, count in counts.items()) if status == "Parsed" else np.nan
    dbe = 1 + c - (h + halogens) / 2 + n / 2 if counts else np.nan
    hetero = sum(count for element, count in counts.items() if element not in {"C", "H"})
    hc = h / c if c else np.nan
    oc = o / c if c else np.nan
    nc = n / c if c else np.nan
    hetero_c = hetero / c if c else np.nan
    polarity_score = 50.0
    if c:
        polarity_score = 100 * (1 - math.exp(-1.8 * (o + 0.9 * n + 1.4 * p + 0.7 * s) / max(c, 1)))
        polarity_score += 10 if counts.get("Na", 0) or counts.get("K", 0) else 0
        polarity_score = min(100.0, max(0.0, polarity_score))
    polarity_class = (
        "Unknown"
        if not counts
        else "High formula-polarity proxy"
        if polarity_score >= 65
        else "Intermediate / amphiphilic proxy"
        if polarity_score >= 35
        else "Low formula-polarity proxy"
    )
    lipid_like = bool(c >= 12 and _finite(hc) is not None and hc >= 1.2 and (_finite(oc) or 0) <= 0.45)
    kendrick_mass = exact_mass * 14.0 / 14.01565006446 if math.isfinite(exact_mass) else np.nan
    kendrick_nominal = round(kendrick_mass) if math.isfinite(kendrick_mass) else np.nan
    return {
        "Formula parse status": status,
        "Formula exact monoisotopic mass (Da)": exact_mass,
        "Formula DBE": dbe,
        "Formula C count": c,
        "Formula H count": h,
        "Formula N count": n,
        "Formula O count": o,
        "Formula P count": p,
        "Formula S count": s,
        "Formula halogen count": halogens,
        "Formula heteroatom count": hetero,
        "Formula H/C": hc,
        "Formula O/C": oc,
        "Formula N/C": nc,
        "Formula heteroatom/C": hetero_c,
        "Formula polarity proxy (0-100)": polarity_score if counts else np.nan,
        "Formula polarity class": polarity_class,
        "Lipid-like formula flag": lipid_like,
        "Kendrick mass (CH2)": kendrick_mass,
        "Kendrick mass defect (CH2)": kendrick_nominal - kendrick_mass if math.isfinite(kendrick_mass) else np.nan,
    }


def chemical_family(name: Any, description: Any = "") -> str:
    text = f"{name or ''} {description or ''}"
    rules = (
        (r"\b(PC|LPC|phosphatidylcholine)\b", "Glycerophosphocholine"),
        (r"\b(PE|LPE|phosphatidylethanolamine)\b", "Glycerophosphoethanolamine"),
        (r"\b(PG|PI|PS|PA)\b|phosphatidyl", "Other glycerophospholipid"),
        (r"\b(TG|triacylglycer|triglycer)\b", "Triacylglycerol"),
        (r"\b(DG|diacylglycer)\b", "Diacylglycerol"),
        (r"\b(Cer|ceramide|sphingo|ganglioside|SM\b)", "Sphingolipid"),
        (r"cholest|sterol|steroid|bile acid", "Sterol / steroid / bile acid"),
        (r"fatty acid|eicosanoid|prostagland|leukotriene", "Fatty acid / oxylipin"),
        (r"amino acid|peptide", "Amino acid / peptide"),
        (r"nucleoside|nucleotide|purine|pyrimidine", "Nucleotide / nucleoside"),
        (r"carbohydrate|saccharide|hexose|pentose|sugar", "Carbohydrate"),
        (r"organic acid|carboxylic acid", "Organic acid"),
        (r"drug|pharmaceutical|antibiotic|antifungal", "Drug / pharmaceutical"),
        (r"pesticide|herbicide|pollutant|plasticizer|industrial", "Environmental chemical"),
    )
    for pattern, label in rules:
        if re.search(pattern, text, re.IGNORECASE):
            return label
    return "Unclassified small molecule"


def adduct_polarity(adducts: Any) -> str:
    text = str(adducts or "")
    positive = bool(re.search(r"M\+H|M\+Na|M\+K|M\+NH4|\]\+", text, re.IGNORECASE))
    negative = bool(re.search(r"M-H|M\+HCOO|M\+CH3COO|M-HCOO|\]-", text, re.IGNORECASE))
    if positive and negative:
        return "Mixed / ambiguous"
    if positive:
        return "Positive-compatible"
    if negative:
        return "Negative-compatible"
    return "Unknown"


def enrich_chemical_properties(frame: pd.DataFrame, ion_mode: str = "mixed") -> pd.DataFrame:
    output = frame.copy()
    properties = pd.DataFrame([formula_properties(value) for value in output.get("Formula", pd.Series([""] * len(output)))], index=output.index)
    for column in properties.columns:
        output[column] = properties[column]
    output["Chemical family (text rule)"] = [
        chemical_family(name, description)
        for name, description in zip(output.get("Canonical name", ""), output.get("Accepted Description", ""))
    ]
    output["Adduct polarity class"] = output.get("Adducts", pd.Series([""] * len(output))).map(adduct_polarity)
    mode = str(ion_mode or "mixed").casefold()
    output["Ion-mode/adduct consistency"] = output["Adduct polarity class"].map(
        lambda value: "Consistent"
        if mode == "mixed" or value == "Unknown" or mode in value.casefold()
        else "Review"
    )
    neutral = pd.to_numeric(output.get("Neutral mass (Da)"), errors="coerce")
    formula_mass = pd.to_numeric(output["Formula exact monoisotopic mass (Da)"], errors="coerce")
    output["Formula-neutral mass difference (Da)"] = neutral - formula_mass
    output["Formula-neutral mass difference (ppm)"] = (neutral - formula_mass) / formula_mass * 1_000_000
    output["Chemical-property limitation"] = (
        "Formula-derived descriptors do not encode isomerism, pKa, measured logP, or measured recovery; polarity is a review proxy."
    )
    return output
