"""
test_ncbi_client.py
Tests for the NCBI/UniProt protein client.

Like the PubChem tests, these never touch the network: requests.get is
replaced with fakes returning canned UniProt and E-utilities payloads.
Run with: pytest
"""

import pytest
import requests

from toxreadout import ncbi_client


class FakeResponse:
    def __init__(self, status_code=200, json_data=None):
        self.status_code = status_code
        self._json = json_data

    def json(self):
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


# A trimmed-down but realistic UniProt record for insulin (P01308).
UNIPROT_INSULIN = {
    "primaryAccession": "P01308",
    "proteinDescription": {
        "recommendedName": {"fullName": {"value": "Insulin"}}
    },
    "genes": [{"geneName": {"value": "INS"}}],
    "organism": {
        "scientificName": "Homo sapiens",
        "commonName": "Human",
        "taxonId": 9606,
        "lineage": ["Eukaryota", "Metazoa", "Chordata", "Mammalia", "Primates"],
    },
    "sequence": {"length": 110, "molWeight": 11981},
    "comments": [
        {
            "commentType": "FUNCTION",
            "texts": [{"value": "Insulin decreases blood glucose concentration."}],
        }
    ],
}


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    """Skip real waiting so tests run instantly."""
    monkeypatch.setattr(ncbi_client.time, "sleep", lambda *_: None)


def test_parse_uniprot_extracts_all_fields():
    """The parser should pull every field we care about out of the record."""
    parsed = ncbi_client._parse_uniprot(UNIPROT_INSULIN)
    assert parsed["accession"] == "P01308"
    assert parsed["protein_name"] == "Insulin"
    assert parsed["gene"] == "INS"
    assert parsed["organism_scientific"] == "Homo sapiens"
    assert parsed["organism_common"] == "Human"
    assert parsed["taxon_id"] == 9606
    assert parsed["length"] == 110
    assert parsed["molecular_weight"] == 11981
    assert "glucose" in parsed["function"]


def test_get_uniprot_record_not_found_raises(monkeypatch):
    """A 404 from UniProt should raise a clear ValueError."""
    monkeypatch.setattr(
        ncbi_client.requests, "get", lambda *a, **k: FakeResponse(404, None)
    )
    with pytest.raises(ValueError):
        ncbi_client.get_uniprot_record("NOT_AN_ACCESSION")


def test_get_gi_number_parses_idlist(monkeypatch):
    """The first UID from esearch is the GI number."""
    payload = {"esearchresult": {"idlist": ["124617"]}}
    monkeypatch.setattr(
        ncbi_client.requests, "get", lambda *a, **k: FakeResponse(200, payload)
    )
    assert ncbi_client.get_gi_number("P01308") == "124617"


def test_get_gi_number_no_match_returns_none(monkeypatch):
    """An empty idlist should yield None rather than an error."""
    payload = {"esearchresult": {"idlist": []}}
    monkeypatch.setattr(
        ncbi_client.requests, "get", lambda *a, **k: FakeResponse(200, payload)
    )
    assert ncbi_client.get_gi_number("P01308") is None


def test_fetch_protein_combines_sources(monkeypatch):
    """fetch_protein should merge UniProt metadata with the NCBI GI number."""
    def fake_get(url, *a, **k):
        if "uniprot" in url:
            return FakeResponse(200, UNIPROT_INSULIN)
        return FakeResponse(200, {"esearchresult": {"idlist": ["124617"]}})

    monkeypatch.setattr(ncbi_client.requests, "get", fake_get)

    record = ncbi_client.fetch_protein("P01308")
    assert record["protein_name"] == "Insulin"
    assert record["gi_number"] == "124617"
