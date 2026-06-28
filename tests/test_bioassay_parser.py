"""
test_bioassay_parser.py
Tests for the BioAssay parser. No network involved — these feed sample
assay rows straight in and check the parsing, grouping, and summary math.
Run with: pytest
"""

from toxreadout import bioassay_parser as bp


# A small set of raw rows shaped like pubchem_client.get_bioassays output.
SAMPLE_RAW = [
    {"aid": "1", "activity_outcome": "Active", "target_accession": "P29274",
     "activity_value_um": 1.2, "activity_name": "IC50", "assay_type": "Confirmatory"},
    {"aid": "2", "activity_outcome": "Active", "target_accession": "P29274",
     "activity_value_um": 48.5, "activity_name": "EC50", "assay_type": "Screening"},
    {"aid": "3", "activity_outcome": "Inactive", "target_accession": "P29274",
     "activity_value_um": None, "activity_name": None, "assay_type": "Confirmatory"},
    {"aid": "4", "activity_outcome": "Inconclusive", "target_accession": "Q9Y5N1",
     "activity_value_um": None, "activity_name": None, "assay_type": "Other"},
    {"aid": "5", "activity_outcome": "Inactive", "target_accession": "",
     "activity_value_um": None, "activity_name": None, "assay_type": "Other"},
]


def test_normalize_outcome():
    assert bp.normalize_outcome("Active") == "ACTIVE"
    assert bp.normalize_outcome("  inactive ") == "INACTIVE"
    assert bp.normalize_outcome("Probe") == "PROBE"
    assert bp.normalize_outcome("weird value") == "UNSPECIFIED"
    assert bp.normalize_outcome(None) == "UNSPECIFIED"


def test_parse_results_normalizes_and_keeps_fields():
    parsed = bp.parse_results(SAMPLE_RAW)
    assert len(parsed) == 5
    assert parsed[0]["outcome"] == "ACTIVE"
    assert parsed[0]["target_accession"] == "P29274"
    assert parsed[2]["outcome"] == "INACTIVE"


def test_group_by_target_skips_untargeted():
    parsed = bp.parse_results(SAMPLE_RAW)
    groups = bp.group_by_target(parsed)
    # Two real targets; the blank-target row (aid 5) is dropped.
    assert set(groups.keys()) == {"P29274", "Q9Y5N1"}
    assert len(groups["P29274"]) == 3
    assert len(groups["Q9Y5N1"]) == 1


def test_activity_summary_math():
    parsed = bp.parse_results(SAMPLE_RAW)
    summary = bp.get_activity_summary(parsed)
    assert summary["total_assays"] == 5
    assert summary["active"] == 2
    assert summary["inactive"] == 2
    assert summary["inconclusive"] == 1
    assert summary["active_ratio"] == 0.4  # 2 of 5
    # Potency range spans only active assays with a value.
    assert summary["potency_range_um"] == {"min": 1.2, "max": 48.5}
    assert summary["assay_types"] == ["Confirmatory", "Other", "Screening"]
    assert summary["potency_measures"] == ["EC50", "IC50"]


def test_summary_with_no_potency_values():
    rows = [
        {"aid": "1", "activity_outcome": "Active", "target_accession": "X",
         "activity_value_um": None, "activity_name": None, "assay_type": "Other"},
    ]
    summary = bp.get_activity_summary(bp.parse_results(rows))
    assert summary["potency_range_um"] is None
    assert summary["active_ratio"] == 1.0


def test_summarize_by_target_end_to_end():
    by_target = bp.summarize_by_target(SAMPLE_RAW)
    assert set(by_target.keys()) == {"P29274", "Q9Y5N1"}
    p = by_target["P29274"]
    assert p["total_assays"] == 3
    assert p["active"] == 2
    assert p["active_ratio"] == 0.67  # 2 of 3, rounded
