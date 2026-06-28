"""
test_readout.py
Tests for the readout assembler. No network — sample upstream pieces are
assembled and the resulting structure / renderings are checked.
Run with: pytest
"""

import json

from toxreadout import readout as ro


COMPOUND = {
    "cid": 2519,
    "common_name": "Caffeine",
    "iupac_name": "1,3,7-trimethylpurine-2,6-dione",
    "molecular_formula": "C8H10N4O2",
    "molecular_weight": 194.19,
    "bioassay_results": [],
}

PROTEIN = {
    "accession": "P29274",
    "protein_name": "Adenosine receptor A2a",
    "gene": "ADORA2A",
    "organism_scientific": "Homo sapiens",
}

INTERACTION = {
    "match_found": True,
    "total_assays": 39,
    "active_assays": 6,
    "active_ratio": 0.15,
    "potency_range_um": {"min": 2.48, "max": 9.6},
    "assay_types": ["Confirmatory", "Other"],
}

RISK = {
    "score": 0.41,
    "concern_level": "MODERATE",
    "reasoning": "6 of 39 assays active (15%); most potent activity 2.48 uM.",
    "confidence": "high",
    "data_source": "database_lookup",
}

LITERATURE = {
    "query": "caffeine AND adenosine receptor A2a",
    "total_results": 31,
    "top_references": [{"pmid": "1", "title": "A paper", "url": "http://x/1/"}],
}


def _sample_readout(literature=LITERATURE):
    target = ro.build_target_entry(PROTEIN, INTERACTION, RISK)
    return ro.assemble_readout(
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        COMPOUND,
        [target],
        literature=literature,
        timestamp="2026-06-28T12:00:00Z",
    )


def test_build_target_entry_shape():
    entry = ro.build_target_entry(PROTEIN, INTERACTION, RISK)
    assert entry["protein"]["accession"] == "P29274"
    assert entry["protein"]["name"] == "Adenosine receptor A2a"
    assert entry["interaction"]["active_assays"] == 6
    assert entry["risk_assessment"]["concern_level"] == "MODERATE"


def test_assemble_readout_structure():
    readout = _sample_readout()
    data = readout.to_dict()
    assert data["timestamp"] == "2026-06-28T12:00:00Z"
    assert data["compound"]["cid"] == 2519
    assert data["compound"]["input_smiles"].startswith("CN1")
    assert len(data["targets"]) == 1
    assert data["literature"]["pubmed_hits"] == 31


def test_to_json_round_trips():
    readout = _sample_readout()
    parsed = json.loads(readout.to_json())
    assert parsed["compound"]["common_name"] == "Caffeine"
    assert parsed["targets"][0]["risk_assessment"]["score"] == 0.41


def test_to_pretty_string_contains_key_facts():
    text = _sample_readout().to_pretty_string()
    assert "TOXICOLOGY READOUT" in text
    assert "Caffeine" in text
    assert "P29274" in text
    assert "MODERATE" in text
    assert "PubMed hits: 31" in text


def test_multiple_targets_and_no_match():
    no_match_interaction = {
        "match_found": False,
        "total_assays": 0,
        "active_assays": 0,
        "active_ratio": 0.0,
        "potency_range_um": None,
        "assay_types": [],
    }
    no_match_risk = {
        "score": 0.0,
        "concern_level": "UNKNOWN",
        "reasoning": "No data.",
        "confidence": "none",
        "data_source": "none",
    }
    insulin = {"accession": "P01308", "protein_name": "Insulin",
               "gene": "INS", "organism_scientific": "Homo sapiens"}

    targets = [
        ro.build_target_entry(PROTEIN, INTERACTION, RISK),
        ro.build_target_entry(insulin, no_match_interaction, no_match_risk),
    ]
    readout = ro.assemble_readout("smiles", COMPOUND, targets, timestamp="t")
    assert len(readout.to_dict()["targets"]) == 2
    text = readout.to_pretty_string()
    assert "no BioAssay data for this pair" in text


def test_literature_optional():
    readout = _sample_readout(literature=None)
    assert readout.to_dict()["literature"] is None
    # Pretty string should simply omit the literature section.
    assert "LITERATURE" not in readout.to_pretty_string()
