from __future__ import annotations

import math
import re
from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd

from .engine import clean

SOURCE_CLASSES = (
    "Human endogenous",
    "Microbiome-derived",
    "Host–microbe co-metabolic",
    "Food-derived",
    "Drug-derived",
    "Environmental-derived",
)
SOURCE_KEYS = ("human", "microbe", "co", "food", "drug", "environment")


def _softmax(scores: dict[str, float], temperature: float = 1.35) -> dict[str, float]:
    maximum = max(scores.values())
    exponentials = {key: math.exp((value - maximum) / temperature) for key, value in scores.items()}
    total = sum(exponentials.values())
    return {key: value / total * 100 for key, value in exponentials.items()}


def _descriptor(description: str, pattern: str) -> float | None:
    match = re.search(pattern, description, re.IGNORECASE)
    if not match:
        return None
    try:
        value = float(match.group(1))
        return value if math.isfinite(value) else None
    except ValueError:
        return None


def _normalized_entropy(probabilities: dict[str, float]) -> float:
    values = np.asarray([value / 100 for value in probabilities.values() if value > 0], dtype=float)
    return float(-(values * np.log(values)).sum() / math.log(len(SOURCE_KEYS))) if len(values) > 1 else 0.0


def _source_priors(settings: dict[str, Any]) -> tuple[dict[str, float], list[str]]:
    system = clean(settings.get("biological_system") or "human").casefold()
    matrix = clean(settings.get("sample_matrix")).casefold()
    assay = clean(settings.get("study_type") or settings.get("assay_type") or "metabolomics").casefold()
    if re.search(r"bacter|microbial culture|fung", system):
        priors = {"human": 0.08, "microbe": 0.62, "co": 0.08, "food": 0.06, "drug": 0.07, "environment": 0.09}
        notes = ["bacterial/microbial-system prior"]
    elif re.search(r"human.*microbiome|mixed host|fecal|stool|gut", f"{system} {matrix}"):
        priors = {"human": 0.22, "microbe": 0.31, "co": 0.25, "food": 0.10, "drug": 0.05, "environment": 0.07}
        notes = ["host–microbiome matrix prior"]
    elif re.search(r"cell line|cell culture|virus|viral", system):
        priors = {"human": 0.52, "microbe": 0.06, "co": 0.05, "food": 0.06, "drug": 0.18, "environment": 0.13}
        notes = ["cell/viral-exposure system prior"]
    elif re.search(r"environment|water|soil|plant", system):
        priors = {"human": 0.10, "microbe": 0.23, "co": 0.08, "food": 0.10, "drug": 0.12, "environment": 0.37}
        notes = ["environmental-system prior"]
    else:
        priors = {"human": 0.40, "microbe": 0.12, "co": 0.15, "food": 0.13, "drug": 0.09, "environment": 0.11}
        notes = ["human-sample prior"]
    if re.search(r"urine|saliva|feces|stool|gut", matrix):
        priors["microbe"] *= 1.25
        priors["co"] *= 1.25
        notes.append("microbiome-accessible matrix modifier")
    if re.search(r"serum|plasma|tissue|placenta", matrix):
        priors["human"] *= 1.15
        notes.append("human biofluid/tissue modifier")
    if assay == "lipidomics":
        priors["human"] *= 1.12
        notes.append("lipidomics assay modifier")
    total = sum(priors.values())
    return {key: value / total for key, value in priors.items()}, notes


def _phase_model(row: pd.Series, settings: dict[str, Any], mol_logp: float | None, tpsa: float | None) -> dict[str, Any]:
    description = clean(row.get("Accepted Description"))
    family = clean(row.get("Chemical family (text rule)"))
    polarity = row.get("Formula polarity proxy (0-100)")
    polarity = float(polarity) if pd.notna(polarity) else None
    lipid_like = bool(row.get("Lipid-like formula flag"))
    method = clean(settings.get("extraction_method") or "unknown").casefold()
    analyzed = clean(settings.get("analyzed_phase") or "unknown").casefold()
    explicit = ""
    logits = {"organic": 0.0, "aqueous": 0.0, "retained": -0.5, "flowthrough": -0.5, "pellet": -0.6, "uncertain": 0.3}

    if re.search(r"MTBE_DOMINANT|UPPER_MTBE|ORGANIC_DOMINANT", description, re.IGNORECASE):
        logits["organic"] += 4.3
        explicit = "description explicitly assigns organic/MTBE-dominant behavior"
    elif re.search(r"AQUEOUS_DOMINANT|LOWER_AQUEOUS", description, re.IGNORECASE):
        logits["aqueous"] += 4.3
        explicit = "description explicitly assigns aqueous-dominant behavior"
    elif re.search(r"CROSS_PHASE_OR_INTERPHASE|INTERPHASE_OR_SPLIT", description, re.IGNORECASE):
        logits["organic"] += 1.2
        logits["aqueous"] += 1.2
        logits["pellet"] += 1.2
        explicit = "description flags cross-phase/interphase behavior"
    elif re.search(r"PROTEIN_PELLET|PELLET_DOMINANT", description, re.IGNORECASE):
        logits["pellet"] += 4.3
        explicit = "description explicitly flags pellet association"

    property_notes: list[str] = []
    if mol_logp is not None:
        logits["organic"] += 0.68 * (mol_logp - 1.5)
        logits["aqueous"] -= 0.55 * (mol_logp - 1.0)
        property_notes.append(f"parsed MolLogP {mol_logp:g}")
    if tpsa is not None:
        logits["aqueous"] += 0.011 * (tpsa - 35)
        logits["organic"] -= 0.008 * max(tpsa - 35, 0)
        property_notes.append(f"parsed TPSA {tpsa:g} Å²")
    if polarity is not None:
        logits["aqueous"] += (polarity - 45) / 28
        logits["organic"] += (45 - polarity) / 30
        property_notes.append(f"formula-polarity proxy {polarity:.1f}/100")
    if lipid_like or re.search(r"lipid|glycer|sphingo|sterol|fatty", family, re.IGNORECASE):
        logits["organic"] += 1.35
        property_notes.append("lipid-class/formula cue")
    if clean(row.get("Adduct polarity class")) == "Mixed / ambiguous":
        logits["uncertain"] += 0.4

    if method in {"mtbe", "mtbe-biphasic", "folch", "bligh-dyer", "ethyl-acetate", "butanol", "generic-lle", "liquid-liquid"}:
        logits["retained"] -= 1.0
        logits["flowthrough"] -= 1.0
        if method.startswith("mtbe"):
            method_note = "MTBE: organic-rich layer is upper"
            organic_label, aqueous_label = "Upper organic", "Lower aqueous"
        elif method in {"folch", "bligh-dyer"}:
            method_note = "chloroform-based LLE: organic-rich layer is lower"
            organic_label, aqueous_label = "Lower organic", "Upper aqueous"
        else:
            method_note = "declared liquid–liquid extraction"
            organic_label, aqueous_label = "Organic-rich phase", "Aqueous-rich phase"
    elif method in {"spe-hlb", "spe-c18", "spe-reversed-phase"}:
        logits["retained"] += 1.5 + logits["organic"] * 0.5
        logits["flowthrough"] += 0.5 + logits["aqueous"] * 0.45
        method_note = "reversed-phase/HLB SPE; retention depends on sorbent, loading solvent, pH, and wash/elution"
        organic_label, aqueous_label = "Hydrophobic tendency", "Polar tendency"
    elif method in {"spe-hilic", "spe-normal-phase"}:
        logits["retained"] += 1.3 + logits["aqueous"] * 0.45
        logits["flowthrough"] += 0.4 + logits["organic"] * 0.35
        method_note = "HILIC/normal-phase SPE; polar retention is method-dependent"
        organic_label, aqueous_label = "Low-polarity tendency", "Polar tendency"
    elif method in {"spe-cation", "spe-anion", "spe-mixed-mode"}:
        logits["retained"] += 0.9
        logits["uncertain"] += 1.2
        method_note = "ion-exchange/mixed-mode SPE; pKa and loading/elution pH are not encoded"
        organic_label, aqueous_label = "Hydrophobic tendency", "Polar/ionic tendency"
    elif method in {"protein-precipitation", "ppt", "methanol-precipitation", "acetonitrile-precipitation"}:
        logits["aqueous"] += 0.5
        logits["pellet"] -= 0.3
        method_note = "protein precipitation; most small molecules remain in supernatant, with compound-specific losses"
        organic_label, aqueous_label = "Organic-soluble tendency", "Supernatant/polar tendency"
    elif method in {"monophasic-polar", "polar", "methanol-water", "acetonitrile-water"}:
        logits["aqueous"] += 1.0
        method_note = "monophasic polar extraction"
        organic_label, aqueous_label = "Low-polarity tendency", "Polar extract"
    else:
        logits["uncertain"] += 1.0
        method_note = "extraction method unspecified"
        organic_label, aqueous_label = "Organic-rich tendency", "Aqueous-rich tendency"

    probabilities = _softmax(logits, temperature=1.15)
    labels = {
        "organic": organic_label,
        "aqueous": aqueous_label,
        "retained": "SPE retained / eluate",
        "flowthrough": "SPE flow-through",
        "pellet": "Interphase / protein pellet",
        "uncertain": "Uncertain",
    }
    dominant_key = max(probabilities, key=probabilities.get)
    phase_map = {
        "upper-organic": "organic" if method.startswith("mtbe") or method in {"ethyl-acetate", "butanol", "generic-lle", "liquid-liquid"} else None,
        "lower-organic": "organic",
        "organic": "organic",
        "lower-aqueous": "aqueous",
        "upper-aqueous": "aqueous",
        "aqueous": "aqueous",
        "spe-eluate": "retained",
        "eluate": "retained",
        "spe-flowthrough": "flowthrough",
        "flow-through": "flowthrough",
        "interphase-pellet": "pellet",
        "pellet": "pellet",
        "supernatant": "aqueous" if method in {"protein-precipitation", "ppt", "methanol-precipitation", "acetonitrile-precipitation"} else None,
    }
    analyzed_key = phase_map.get(analyzed)
    analyzed_likelihood = probabilities.get(analyzed_key, np.nan) if analyzed_key else np.nan
    if analyzed == "both":
        analyzed_likelihood = min(100.0, probabilities["organic"] + probabilities["aqueous"])
    descriptor_count = int(mol_logp is not None) + int(tpsa is not None) + int(polarity is not None) + int(lipid_like)
    confidence = "Moderate" if explicit or descriptor_count >= 3 else "Low" if descriptor_count else "Very low"
    evidence = "; ".join(filter(None, [explicit, ", ".join(property_notes), method_note]))
    return {
        "Organic-rich extract likelihood (%)": probabilities["organic"],
        "Aqueous/polar extract likelihood (%)": probabilities["aqueous"],
        "SPE retained/eluate likelihood (%)": probabilities["retained"],
        "SPE flow-through likelihood (%)": probabilities["flowthrough"],
        "Interphase/pellet likelihood (%)": probabilities["pellet"],
        "Phase uncertainty (%)": probabilities["uncertain"],
        "Upper-organic phase likelihood (%)": probabilities["organic"] if method.startswith("mtbe") else np.nan,
        "Lower-aqueous phase likelihood (%)": probabilities["aqueous"] if method.startswith("mtbe") else np.nan,
        "Dominant predicted phase": labels[dominant_key],
        "Analyzed-phase likelihood (%)": analyzed_likelihood,
        "Phase confidence": confidence,
        "Phase evidence": evidence,
        "Phase-model coverage": f"{descriptor_count} chemical descriptor/cue dimensions; {'explicit behavior cue' if explicit else 'no explicit behavior cue'}",
    }


def score_provenance(
    frame: pd.DataFrame,
    settings: dict[str, Any],
    progress: Callable[[int, str], None] | None = None,
) -> pd.DataFrame:
    update = progress or (lambda _percent, _message: None)
    update(3, "Initializing assay, biological-system, matrix, and exposure priors…")
    priors, prior_notes = _source_priors(settings)
    study_type = clean(settings.get("study_type") or settings.get("assay_type") or "metabolomics").casefold()
    sample_matrix = clean(settings.get("sample_matrix"))
    biological_system = clean(settings.get("biological_system") or "human")
    exposure = clean(settings.get("exposure_context") or "none")
    output_rows: list[dict[str, Any]] = []

    total = len(frame)
    update(8, f"Scoring source and extraction evidence for {total:,} compounds…")
    for position, (_, row) in enumerate(frame.iterrows(), 1):
        description = clean(row.get("Accepted Description"))
        canonical = clean(row.get("Canonical name"))
        text = f"{canonical} {description}"
        logits = {key: math.log(max(value, 1e-9)) for key, value in priors.items()}
        evidence: list[str] = []
        evidence_strength = 0.0
        direct_evidence = 0

        def add(key: str, weight: float, note: str, direct: bool = False) -> None:
            nonlocal evidence_strength, direct_evidence
            logits[key] += weight
            evidence.append(note)
            evidence_strength += abs(weight)
            direct_evidence += int(direct)

        hmdb_endogenous = bool(re.search(r"HMDB carries an endogenous ontology assertion", description, re.IGNORECASE))
        observed = bool(re.search(r"status '(detected|quantified)'", description, re.IGNORECASE)) and not bool(re.search(r"not observed", description, re.IGNORECASE))
        if hmdb_endogenous:
            add("human", 3.6 if observed else 2.2, "HMDB endogenous ontology assertion" + (" with detected/quantified status" if observed else ""), observed)
        if re.search(r"human-only metabolite|endogenous human", description, re.IGNORECASE):
            add("human", 4.5, "explicit human-only/endogenous cue", True)
        if re.search(r"microbial-only|only be produced by microbial", description, re.IGNORECASE):
            add("microbe", 5.0, "explicit microbial-only cue", True)
        if re.search(r"MiMeDB carries a microbiome relation|microbiome relation", description, re.IGNORECASE):
            add("microbe", 2.5, "MiMeDB microbiome relation")
            logits["co"] += 0.7
        if re.search(r"co[- ]?metabol|host[- ]microb|microbial/host", description, re.IGNORECASE):
            add("co", 4.7, "host–microbe co-metabolic cue", True)
        if re.search(r"HMDB_ALSO_FOOD_SOURCE", description, re.IGNORECASE):
            add("food", 1.0, "secondary HMDB food-source annotation")
        if re.search(r"(?<!ALSO_)\bfood source\b|dietary|nutritional|food-derived", description, re.IGNORECASE):
            add("food", 4.0, "explicit food/dietary source cue", True)
        if re.search(r"\bdrug\b|pharmaceutical|medication|therapeutic|xenobiotic metabolite|antibiotic|antifungal", description, re.IGNORECASE):
            add("drug", 4.1, "drug/pharmaceutical cue", True)
        if re.search(r"environmental|pollutant|pesticide|herbicide|industrial|plasticizer|aquatic toxicity|embryotoxic", description, re.IGNORECASE):
            add("environment", 4.3, "environmental/pollutant cue", True)

        association_only = bool(re.search(r"associated_not_observed|NOT_OBSERVED|EXPECTED_OR_PREDICTED_NOT_OBSERVED|not proof of occurrence", description, re.IGNORECASE))
        if association_only:
            evidence.append("association/expected status—not direct observation")
            evidence_strength *= 0.65
            direct_evidence = 0
        if study_type == "lipidomics" and re.search(r"\b(PC|PE|TG|DG|Cer|SM|phosphatid|sphingo|cholest)\b", text, re.IGNORECASE):
            logits["human"] += 0.55
            evidence.append("lipid-class assay-context modifier")
        if re.search(r"drug|treatment|dose", exposure, re.IGNORECASE):
            logits["drug"] += 0.45
        if re.search(r"environment|chemical|pollut", exposure, re.IGNORECASE):
            logits["environment"] += 0.45

        probabilities = _softmax(logits)
        classes = list(zip(SOURCE_CLASSES, [probabilities[key] for key in SOURCE_KEYS]))
        classes.sort(key=lambda item: item[1], reverse=True)
        margin = classes[0][1] - classes[1][1]
        entropy = _normalized_entropy(probabilities)
        confidence = (
            "High"
            if direct_evidence and classes[0][1] >= 70 and margin >= 30
            else "Moderate"
            if evidence_strength >= 2.5 and classes[0][1] >= 50 and margin >= 15 and not association_only
            else "Low"
        )
        if not evidence:
            confidence = "Low"

        mol_logp = _descriptor(description, r"(?:RDKit\s+)?MolLogP\s+(-?\d+(?:\.\d+)?)")
        tpsa = _descriptor(description, r"TPSA\s+(-?\d+(?:\.\d+)?)")
        phase = _phase_model(row, settings, mol_logp, tpsa)
        result = row.to_dict()
        result.update(
            {
                "Study type": "Lipidomics" if study_type == "lipidomics" else "Metabolomics",
                "Biological system": biological_system,
                "Sample matrix": sample_matrix,
                "Exposure context": exposure,
                "Extraction method": settings.get("extraction_method", "unknown"),
                "Analyzed phase": settings.get("analyzed_phase", "unknown"),
                "Human endogenous likelihood (%)": probabilities["human"],
                "Microbiome-derived likelihood (%)": probabilities["microbe"],
                "Host–microbe co-metabolic likelihood (%)": probabilities["co"],
                "Food-derived likelihood (%)": probabilities["food"],
                "Drug-derived likelihood (%)": probabilities["drug"],
                "Environmental-derived likelihood (%)": probabilities["environment"],
                "Primary source class": classes[0][0],
                "Source confidence": confidence,
                "Source class margin (percentage points)": margin,
                "Source model entropy (0-1)": entropy,
                "Direct source evidence count": direct_evidence,
                "Source evidence strength (rule units)": evidence_strength,
                "Source prior context": "; ".join(prior_notes),
                "Source evidence": "; ".join(evidence) if evidence else "No explicit ontology/source cue; result is prior-driven.",
                "Source review priority": "High" if confidence == "Low" or entropy >= 0.75 else "Medium" if confidence == "Moderate" else "Routine",
                "Parsed MolLogP": mol_logp if mol_logp is not None else np.nan,
                "Parsed TPSA (angstrom^2)": tpsa if tpsa is not None else np.nan,
                **phase,
                "Interpretation limitation": "Evidence-weighted heuristic, not a calibrated causal-origin probability or measured extraction recovery. Confirm priority compounds with orthogonal standards, curated database review, blanks, and method-matched recovery experiments.",
            }
        )
        output_rows.append(result)
        interval = max(1, total // 40)
        if position == total or position % interval == 0:
            update(8 + int(position / max(total, 1) * 90), f"Contextually scored {position:,} / {total:,} compounds…")
    update(100, "Contextual source and extraction-location scoring complete")
    return pd.DataFrame(output_rows)
