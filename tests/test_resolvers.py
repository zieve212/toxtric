"""
test_resolvers.py
Tests for the name resolvers. No network: the request helpers in
pubchem_client and ncbi_client are replaced with fakes that route canned
responses by URL.
Run with: pytest
"""

import pytest
import requests

from toxreadout import resolvers


class FakeResponse:
    def __init__(self, status_code=200, json_data=None):
        self.status_code = status_code
        self._json = json_data

    def json(self):
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


def test_resolve_compound(monkeypatch):
    def fake_pubchem_request(url, *a, **k):
        if "/compound/name/" in url:
            return FakeResponse(200, {"IdentifierList": {"CID": [2519]}})
        if "/property/SMILES" in url:
            return FakeResponse(
                200,
                {"PropertyTable": {"Properties": [{"SMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}]}},
            )
        if "/property/MolecularFormula" in url:
            return FakeResponse(
                200,
                {"PropertyTable": {"Properties": [{
                    "MolecularFormula": "C8H10N4O2",
                    "MolecularWeight": "194.19",
                    "IUPACName": "1,3,7-trimethylpurine-2,6-dione",
                }]}},
            )
        if "/description/" in url:
            return FakeResponse(
                200, {"InformationList": {"Information": [{"Title": "Caffeine"}]}}
            )
        raise AssertionError(f"unexpected url: {url}")

    monkeypatch.setattr(resolvers.pubchem_client, "_request", fake_pubchem_request)

    result = resolvers.resolve_compound("caffeine")
    assert result["cid"] == 2519
    assert result["common_name"] == "Caffeine"
    assert result["smiles"] == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    assert result["molecular_formula"] == "C8H10N4O2"


def test_resolve_compound_not_found(monkeypatch):
    monkeypatch.setattr(
        resolvers.pubchem_client, "_request", lambda *a, **k: FakeResponse(404, None)
    )
    with pytest.raises(ValueError):
        resolvers.resolve_compound("notacompoundxyz")


def test_resolve_protein(monkeypatch):
    payload = {
        "results": [
            {
                "primaryAccession": "P29274",
                "uniProtkbId": "AA2AR_HUMAN",
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Adenosine receptor A2a"}}
                },
                "organism": {"scientificName": "Homo sapiens"},
                "genes": [{"geneName": {"value": "ADORA2A"}}],
                "sequence": {"value": "MPIMGSSVYITVELAIA"},
            }
        ]
    }
    monkeypatch.setattr(
        resolvers.ncbi_client, "_request", lambda *a, **k: FakeResponse(200, payload)
    )

    result = resolvers.resolve_protein("adenosine A2A receptor")
    assert result["accession"] == "P29274"
    assert result["protein_name"] == "Adenosine receptor A2a"
    assert result["gene"] == "ADORA2A"
    assert result["organism"] == "Homo sapiens"
    # FASTA should be built from the header and sequence.
    assert result["fasta"].startswith(">sp|P29274|AA2AR_HUMAN")
    assert "MPIMGSSVYITVELAIA" in result["fasta"]


def test_resolve_protein_not_found(monkeypatch):
    monkeypatch.setattr(
        resolvers.ncbi_client, "_request", lambda *a, **k: FakeResponse(200, {"results": []})
    )
    with pytest.raises(ValueError):
        resolvers.resolve_protein("notaproteinxyz")
