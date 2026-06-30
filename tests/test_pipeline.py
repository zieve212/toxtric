"""
test_pipeline.py
Integration tests for the orchestrator. The network-touching calls
(PubChem, NCBI, PubMed) are mocked; everything else — SMILES/FASTA
validation, matching, scoring, assembly — runs for real.
Run with: pytest
"""

import pytest

from toxreadout import pipeline


CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
FASTA_A2A = ">sp|P29274|AA2AR_HUMAN Adenosine receptor A2a\nMPIMGSSVYITVELAIA"

FAKE_COMPOUND = {
    "cid": 2519,
    "common_name": "Caffeine",
    "iupac_name": "1,3,7-trimethylpurine-2,6-dione",
    "molecular_formula": "C8H10N4O2",
    "molecular_weight": 194.19,
    "bioassay_count": 2,
    "bioassay_results": [
        {"aid": "1", "activity_outcome": "Active", "target_accession": "P29274",
         "activity_value_um": 2.48, "activity_name": "IC50", "assay_type": "Confirmatory"},
        {"aid": "2", "activity_outcome": "Inactive", "target_accession": "P29274",
         "activity_value_um": None, "activity_name": None, "assay_type": "Other"},
    ],
}

FAKE_PROTEIN = {
    "accession": "P29274",
    "protein_name": "Adenosine receptor A2a",
    "gene": "ADORA2A",
    "organism_scientific": "Homo sapiens",
}

FAKE_LIT = {
    "query": "Caffeine AND Adenosine receptor A2a",
    "total_results": 31,
    "top_references": [],
}


@pytest.fixture
def mocked_network(monkeypatch):
    monkeypatch.setattr(pipeline.pubchem_client, "fetch_compound",
                        lambda *_a, **_k: FAKE_COMPOUND)
    monkeypatch.setattr(pipeline.ncbi_client, "fetch_protein",
                        lambda *_a, **_k: FAKE_PROTEIN)
    monkeypatch.setattr(pipeline.literature_client, "literature_lookup",
                        lambda *_a, **_k: FAKE_LIT)


def test_full_run_assembles_readout(mocked_network):
    result = pipeline.run_pipeline(CAFFEINE, [FASTA_A2A], timestamp="t")
    data = result.to_dict()

    assert data["compound"]["cid"] == 2519
    assert len(data["targets"]) == 1
    target = data["targets"][0]
    assert target["protein"]["accession"] == "P29274"
    assert target["interaction"]["match_found"] is True
    assert target["interaction"]["total_assays"] == 2
    assert target["interaction"]["active_assays"] == 1
    assert target["risk_assessment"]["concern_level"] in {
        "LOW", "MODERATE", "HIGH", "CRITICAL"
    }
    assert data["literature"]["pubmed_hits"] == 31


def test_progress_callback_receives_messages(mocked_network):
    messages = []
    pipeline.run_pipeline(CAFFEINE, [FASTA_A2A], on_progress=messages.append, timestamp="t")
    joined = " ".join(messages)
    assert "Validating SMILES" in joined
    assert "Assembling readout" in joined
    assert "Done." in messages[-1]


def test_invalid_smiles_raises_pipeline_error(mocked_network):
    with pytest.raises(pipeline.PipelineError) as exc:
        pipeline.run_pipeline("not_a_smiles", [FASTA_A2A])
    assert "SMILES validation failed" in str(exc.value)


def test_fasta_without_accession_raises(mocked_network):
    # A bare sequence has no accession, so there's nothing to look up.
    with pytest.raises(pipeline.PipelineError) as exc:
        pipeline.run_pipeline(CAFFEINE, ["MPIMGSSVYITVELAIA"])
    assert "No protein accession" in str(exc.value)


def test_protein_metadata_failure_is_non_fatal(monkeypatch):
    """If one protein's metadata lookup fails, the run still completes."""
    monkeypatch.setattr(pipeline.pubchem_client, "fetch_compound",
                        lambda *_a, **_k: FAKE_COMPOUND)
    monkeypatch.setattr(pipeline.literature_client, "literature_lookup",
                        lambda *_a, **_k: FAKE_LIT)

    def boom(*_a, **_k):
        raise RuntimeError("NCBI is down")
    monkeypatch.setattr(pipeline.ncbi_client, "fetch_protein", boom)

    result = pipeline.run_pipeline(CAFFEINE, [FASTA_A2A], timestamp="t")
    target = result.to_dict()["targets"][0]
    # Falls back to accession-only metadata, but interaction/risk still there.
    assert target["protein"]["accession"] == "P29274"
    assert target["protein"]["name"] is None
    assert target["interaction"]["match_found"] is True


def test_literature_can_be_skipped(mocked_network):
    result = pipeline.run_pipeline(
        CAFFEINE, [FASTA_A2A], include_literature=False, timestamp="t"
    )
    assert result.to_dict()["literature"] is None
