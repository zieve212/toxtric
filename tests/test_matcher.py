"""
test_matcher.py
Tests for the compound-protein matcher. No network — sample assay rows
are matched against given accessions and the results checked.
Run with: pytest
"""

from toxreadout import matcher


# Raw rows shaped like pubchem_client.get_bioassays output. Note the
# version suffix on one P29274 row, to test normalization.
SAMPLE_RAW = [
    {"aid": "1", "activity_outcome": "Active", "target_accession": "P29274",
     "activity_value_um": 1.2, "activity_name": "IC50", "assay_type": "Confirmatory"},
    {"aid": "2", "activity_outcome": "Active", "target_accession": "P29274.2",
     "activity_value_um": 48.5, "activity_name": "EC50", "assay_type": "Screening"},
    {"aid": "3", "activity_outcome": "Inactive", "target_accession": "P29274",
     "activity_value_um": None, "activity_name": None, "assay_type": "Confirmatory"},
    {"aid": "4", "activity_outcome": "Active", "target_accession": "Q9Y5N1",
     "activity_value_um": 5.0, "activity_name": "Ki", "assay_type": "Other"},
]


def test_normalize_accession():
    assert matcher.normalize_accession("P29274.1") == "P29274"
    assert matcher.normalize_accession("  p29274 ") == "P29274"
    assert matcher.normalize_accession(None) == ""
    assert matcher.normalize_accession("") == ""


def test_match_protein_found():
    parsed = matcher.bp.parse_results(SAMPLE_RAW)
    result = matcher.match_protein(parsed, "P29274")
    assert result["match_found"] is True
    # All three P29274 rows match (including the .2 version-suffixed one).
    assert result["total_assays"] == 3
    assert result["active_assays"] == 2
    assert result["active_ratio"] == 0.67
    assert result["potency_range_um"] == {"min": 1.2, "max": 48.5}


def test_match_protein_not_found():
    parsed = matcher.bp.parse_results(SAMPLE_RAW)
    result = matcher.match_protein(parsed, "P99999")
    assert result["match_found"] is False
    assert result["total_assays"] == 0
    assert result["active_ratio"] == 0.0
    assert result["potency_range_um"] is None


def test_match_is_case_and_version_insensitive():
    parsed = matcher.bp.parse_results(SAMPLE_RAW)
    # Lowercase, version-suffixed user input should still match.
    result = matcher.match_protein(parsed, "p29274.9")
    assert result["match_found"] is True
    assert result["total_assays"] == 3


def test_match_proteins_multiple_targets():
    results = matcher.match_proteins(SAMPLE_RAW, ["P29274", "Q9Y5N1", "P00000"])
    assert len(results) == 3
    assert results[0]["match_found"] is True and results[0]["total_assays"] == 3
    assert results[1]["match_found"] is True and results[1]["total_assays"] == 1
    assert results[2]["match_found"] is False
