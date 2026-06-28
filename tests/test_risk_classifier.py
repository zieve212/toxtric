"""
test_risk_classifier.py
Tests for the risk classifier. Pure math, no network. Verifies each
scoring component, the concern bands, the no-data path, and a full score.
Run with: pytest
"""

import pytest

from toxreadout import risk_classifier as rc


def test_potency_component_log_scale():
    # At/below the most-potent bound -> 1.0; at/above least-potent -> 0.0.
    assert rc._potency_component(0.1) == 1.0
    assert rc._potency_component(0.01) == 1.0
    assert rc._potency_component(100.0) == 0.0
    assert rc._potency_component(1000.0) == 0.0
    # 1 uM sits two log-decades above 0.1 across a three-decade span.
    assert rc._potency_component(1.0) == pytest.approx(2 / 3, abs=0.01)
    assert rc._potency_component(10.0) == pytest.approx(1 / 3, abs=0.01)
    # No value -> 0.0.
    assert rc._potency_component(None) == 0.0


def test_diversity_component_saturates():
    assert rc._diversity_component([]) == 0.0
    assert rc._diversity_component(["a"]) == pytest.approx(1 / 3, abs=0.01)
    assert rc._diversity_component(["a", "b", "c"]) == 1.0
    assert rc._diversity_component(["a", "b", "c", "d"]) == 1.0  # capped


def test_concern_level_bands():
    assert rc._concern_level(0.0) == "LOW"
    assert rc._concern_level(0.29) == "LOW"
    assert rc._concern_level(0.3) == "MODERATE"
    assert rc._concern_level(0.59) == "MODERATE"
    assert rc._concern_level(0.6) == "HIGH"
    assert rc._concern_level(0.79) == "HIGH"
    assert rc._concern_level(0.8) == "CRITICAL"
    assert rc._concern_level(1.0) == "CRITICAL"


def test_classify_no_match_is_unknown():
    result = rc.classify_risk({"protein_accession": "P01308", "match_found": False})
    assert result["concern_level"] == "UNKNOWN"
    assert result["score"] == 0.0
    assert result["confidence"] == "none"
    assert result["data_source"] == "none"


def test_classify_critical_case():
    """All signals maxed should produce a score of 1.0 / CRITICAL."""
    match = {
        "protein_accession": "X",
        "match_found": True,
        "total_assays": 20,
        "active_assays": 20,
        "active_ratio": 1.0,
        "potency_range_um": {"min": 0.1, "max": 0.5},
        "assay_types": ["binding", "functional", "cell"],
    }
    result = rc.classify_risk(match)
    assert result["score"] == 1.0
    assert result["concern_level"] == "CRITICAL"
    assert result["confidence"] == "high"
    assert result["data_source"] == "database_lookup"


def test_classify_realistic_caffeine_like_case():
    """A caffeine/A2A-like input should land in the MODERATE band."""
    match = {
        "protein_accession": "P29274",
        "match_found": True,
        "total_assays": 39,
        "active_assays": 6,
        "active_ratio": 0.15,
        "potency_range_um": {"min": 2.48, "max": 9.6},
        "assay_types": ["Confirmatory", "Other"],
    }
    result = rc.classify_risk(match)
    # 0.4*0.15 + 0.35*~0.535 + 0.25*0.667 ~= 0.41
    assert result["score"] == pytest.approx(0.41, abs=0.02)
    assert result["concern_level"] == "MODERATE"
    assert result["confidence"] == "high"
    assert "6 of 39 assays active" in result["reasoning"]
    # Breakdown must be present and weights exposed (transparency).
    assert result["score_breakdown"]["weights"]["active_ratio"] == 0.40


def test_score_always_within_bounds():
    match = {
        "match_found": True,
        "total_assays": 1,
        "active_assays": 1,
        "active_ratio": 1.0,
        "potency_range_um": {"min": 0.0001, "max": 0.0001},
        "assay_types": ["a", "b", "c", "d", "e"],
    }
    result = rc.classify_risk(match)
    assert 0.0 <= result["score"] <= 1.0
