"""
risk_classifier.py
Turns a matched compound-protein interaction (from Phase 3b) into a
transparent risk score, a concern level, and plain-English reasoning.

The score (0.0-1.0) is a weighted blend of three signals:
    active ratio   x 0.40   how often the compound was active
    potency        x 0.35   how low the effective concentration is (log scale)
    assay diversity x 0.25  how many distinct assay types confirm it

Every component is computed in its own small function and the breakdown
is returned alongside the score, so the formula is fully inspectable.
No network calls.
"""

import math

# Component weights (must sum to 1.0).
WEIGHT_ACTIVE = 0.40
WEIGHT_POTENCY = 0.35
WEIGHT_DIVERSITY = 0.25

# Potency scale (in uM): at/below "most potent" scores 1.0; at/above
# "least potent" scores 0.0; values in between scale logarithmically.
POTENCY_MOST_POTENT_UM = 0.1
POTENCY_LEAST_POTENT_UM = 100.0

# Number of distinct assay types at which diversity maxes out.
DIVERSITY_SATURATION = 3


def _active_component(active_ratio: float) -> float:
    """The active ratio is already 0-1, so it maps straight through."""
    return max(0.0, min(active_ratio, 1.0))


def _potency_component(min_potency_um: float | None) -> float:
    """
    Lower potency value (more potent) -> higher score, on a log scale.
    Returns 0.0 when no potency value is available.
    """
    if min_potency_um is None:
        return 0.0
    if min_potency_um <= POTENCY_MOST_POTENT_UM:
        return 1.0
    if min_potency_um >= POTENCY_LEAST_POTENT_UM:
        return 0.0
    span = math.log10(POTENCY_LEAST_POTENT_UM) - math.log10(POTENCY_MOST_POTENT_UM)
    return (math.log10(POTENCY_LEAST_POTENT_UM) - math.log10(min_potency_um)) / span


def _diversity_component(assay_types: list) -> float:
    """More distinct assay types -> higher score, capped at saturation."""
    return min(len(assay_types) / DIVERSITY_SATURATION, 1.0)


def _concern_level(score: float) -> str:
    """Map a 0-1 score onto the four named concern bands."""
    if score < 0.3:
        return "LOW"
    if score < 0.6:
        return "MODERATE"
    if score < 0.8:
        return "HIGH"
    return "CRITICAL"


def _confidence(total_assays: int, has_potency: bool) -> str:
    """How much we trust the score, based on the amount of evidence."""
    if total_assays >= 10 and has_potency:
        return "high"
    if total_assays >= 3:
        return "medium"
    return "low"


def _build_reasoning(match_result: dict, min_potency_um: float | None) -> str:
    """Compose a plain-English explanation of the score."""
    total = match_result["total_assays"]
    active = match_result["active_assays"]
    percent = round(match_result["active_ratio"] * 100)
    n_types = len(match_result["assay_types"])

    if min_potency_um is None:
        potency_phrase = "no potency values reported"
    else:
        if min_potency_um <= 1:
            quality = "very potent"
        elif min_potency_um <= 10:
            quality = "moderately potent"
        else:
            quality = "weakly potent"
        potency_phrase = f"most potent activity {min_potency_um} uM ({quality})"

    type_word = "type" if n_types == 1 else "types"
    return (
        f"{active} of {total} assays active ({percent}%); "
        f"{potency_phrase}; {n_types} assay {type_word}."
    )


def classify_risk(match_result: dict) -> dict:
    """
    Score a single protein's matched interaction. Returns the score,
    concern level, reasoning, data source, confidence, and a transparent
    breakdown of how each weighted component contributed.
    """
    # No experimental data means we can't assess risk from a lookup.
    if not match_result.get("match_found"):
        return {
            "score": 0.0,
            "concern_level": "UNKNOWN",
            "reasoning": (
                "No experimental BioAssay data for this compound-protein "
                "pair; risk cannot be assessed from database lookup."
            ),
            "data_source": "none",
            "confidence": "none",
        }

    active_ratio = match_result.get("active_ratio", 0.0)
    potency_range = match_result.get("potency_range_um")
    min_potency = potency_range["min"] if potency_range else None
    assay_types = match_result.get("assay_types", [])

    active_c = _active_component(active_ratio)
    potency_c = _potency_component(min_potency)
    diversity_c = _diversity_component(assay_types)

    score = (
        WEIGHT_ACTIVE * active_c
        + WEIGHT_POTENCY * potency_c
        + WEIGHT_DIVERSITY * diversity_c
    )
    score = round(score, 2)

    return {
        "score": score,
        "concern_level": _concern_level(score),
        "reasoning": _build_reasoning(match_result, min_potency),
        "data_source": "database_lookup",
        "confidence": _confidence(match_result.get("total_assays", 0), min_potency is not None),
        "score_breakdown": {
            "active_ratio_component": round(active_c, 3),
            "potency_component": round(potency_c, 3),
            "diversity_component": round(diversity_c, 3),
            "weights": {
                "active_ratio": WEIGHT_ACTIVE,
                "potency": WEIGHT_POTENCY,
                "diversity": WEIGHT_DIVERSITY,
            },
        },
    }


# Runs only when executed directly: live caffeine -> A2A receptor risk,
# plus the no-data case (insulin).
if __name__ == "__main__":
    try:
        from toxreadout import pubchem_client, matcher
    except ModuleNotFoundError:
        import pubchem_client
        import matcher

    raw = pubchem_client.get_bioassays(2519)  # caffeine
    for match in matcher.match_proteins(raw, ["P29274", "P01308"]):
        risk = classify_risk(match)
        print(f"Protein {match['protein_accession']}:")
        print(f"  Risk score:    {risk['score']}  ({risk['concern_level']})")
        print(f"  Confidence:    {risk['confidence']}")
        print(f"  Reasoning:     {risk['reasoning']}")
        if "score_breakdown" in risk:
            b = risk["score_breakdown"]
            print(f"  Breakdown:     active={b['active_ratio_component']} x{b['weights']['active_ratio']}, "
                  f"potency={b['potency_component']} x{b['weights']['potency']}, "
                  f"diversity={b['diversity_component']} x{b['weights']['diversity']}")
        print()
