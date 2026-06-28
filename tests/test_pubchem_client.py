"""
test_pubchem_client.py
Tests for the PubChem client.

These tests never touch the real network. We replace requests.get and
requests.post with fakes that return canned responses, so we're testing
our own parsing and retry logic — fast, offline, and deterministic.
Run with: pytest
"""

import pytest
import requests

from toxreadout import pubchem_client


class FakeResponse:
    """A stand-in for a requests.Response object."""

    def __init__(self, status_code=200, json_data=None):
        self.status_code = status_code
        self._json = json_data

    def json(self):
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    """Skip all real waiting so tests run instantly."""
    monkeypatch.setattr(pubchem_client.time, "sleep", lambda *_: None)


def test_get_cid_parses_response(monkeypatch):
    """A normal CID response should return the first CID as an int."""
    fake = FakeResponse(200, {"IdentifierList": {"CID": [2519]}})
    monkeypatch.setattr(pubchem_client.requests, "post", lambda *a, **k: fake)

    assert pubchem_client.get_cid("CN1C=NC2=C1C(=O)N(C(=O)N2C)C") == 2519


def test_get_cid_not_found_raises(monkeypatch):
    """A 404 from PubChem should raise a clear ValueError."""
    fake = FakeResponse(404, None)
    monkeypatch.setattr(pubchem_client.requests, "post", lambda *a, **k: fake)

    with pytest.raises(ValueError):
        pubchem_client.get_cid("not_a_real_smiles")


def test_get_properties_parses_and_converts(monkeypatch):
    """Weight should come back as a float, other fields as given."""
    payload = {
        "PropertyTable": {
            "Properties": [
                {
                    "MolecularFormula": "C8H10N4O2",
                    "MolecularWeight": "194.19",
                    "IUPACName": "1,3,7-trimethylpurine-2,6-dione",
                }
            ]
        }
    }
    monkeypatch.setattr(
        pubchem_client.requests, "get", lambda *a, **k: FakeResponse(200, payload)
    )

    props = pubchem_client.get_properties(2519)
    assert props["molecular_formula"] == "C8H10N4O2"
    assert props["molecular_weight"] == 194.19
    assert isinstance(props["molecular_weight"], float)


def test_get_common_name_parses_title(monkeypatch):
    """The record Title should be returned as the common name."""
    payload = {
        "InformationList": {
            "Information": [
                {"CID": 2519, "Title": "Caffeine"},
                {"CID": 2519, "Description": "A methylxanthine alkaloid..."},
            ]
        }
    }
    monkeypatch.setattr(
        pubchem_client.requests, "get", lambda *a, **k: FakeResponse(200, payload)
    )

    assert pubchem_client.get_common_name(2519) == "Caffeine"


def test_get_bioassays_flattens_table(monkeypatch):
    """Assay rows should be flattened into per-assay dicts."""
    payload = {
        "Table": {
            "Columns": {
                "Column": [
                    "AID",
                    "Activity Outcome",
                    "Target Accession",
                    "Target GeneID",
                    "Activity Value [uM]",
                    "Activity Name",
                    "Assay Type",
                    "Assay Name",
                ]
            },
            "Row": [
                {"Cell": ["1234", "Active", "P29274", "135", "12.5", "IC50",
                          "Confirmatory", "A2A binding assay"]},
                {"Cell": ["1235", "Inactive", "", "", "", "", "Other",
                          "Yeast phenotypic screen"]},
            ],
        }
    }
    monkeypatch.setattr(
        pubchem_client.requests, "get", lambda *a, **k: FakeResponse(200, payload)
    )

    assays = pubchem_client.get_bioassays(2519)
    assert len(assays) == 2
    assert assays[0]["aid"] == "1234"
    assert assays[0]["activity_outcome"] == "Active"
    assert assays[0]["target_accession"] == "P29274"
    assert assays[0]["activity_value_um"] == 12.5
    assert assays[0]["activity_name"] == "IC50"
    assert assays[0]["assay_type"] == "Confirmatory"
    # Blank cells become None so the parser can skip untargeted assays.
    assert assays[1]["target_accession"] is None


def test_get_bioassays_no_data_returns_empty(monkeypatch):
    """A 404 (compound has no assays) should return an empty list."""
    monkeypatch.setattr(
        pubchem_client.requests, "get", lambda *a, **k: FakeResponse(404, None)
    )

    assert pubchem_client.get_bioassays(2519) == []


def test_request_retries_on_server_error(monkeypatch):
    """A 503 then a 200 should retry and ultimately succeed."""
    responses = [FakeResponse(503), FakeResponse(503), FakeResponse(200, {"ok": True})]
    calls = {"n": 0}

    def flaky_get(*a, **k):
        resp = responses[calls["n"]]
        calls["n"] += 1
        return resp

    monkeypatch.setattr(pubchem_client.requests, "get", flaky_get)

    result = pubchem_client._request("http://example/test")
    assert result.status_code == 200
    assert calls["n"] == 3  # two failures + one success
