"""
test_cli.py
Tests for the command-line interface. The pipeline is mocked so these
run offline; we verify argument handling, output formats, and file output.
Run with: pytest
"""

import json

import pytest
from click.testing import CliRunner

from toxreadout import cli, readout


SAMPLE_DATA = {
    "timestamp": "2026-06-28T12:00:00Z",
    "compound": {
        "input_smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "cid": 2519,
        "common_name": "Caffeine",
        "iupac_name": "1,3,7-trimethylpurine-2,6-dione",
        "molecular_formula": "C8H10N4O2",
        "molecular_weight": 194.19,
    },
    "targets": [
        {
            "protein": {"accession": "P29274", "name": "Adenosine receptor A2a",
                        "gene": "ADORA2A", "organism": "Homo sapiens"},
            "interaction": {"match_found": True, "total_assays": 39,
                            "active_assays": 6, "active_ratio": 0.15,
                            "potency_range_um": {"min": 2.48, "max": 9.6},
                            "assay_types": ["Confirmatory", "Other"]},
            "risk_assessment": {"score": 0.41, "concern_level": "MODERATE",
                                "reasoning": "6 of 39 assays active.",
                                "confidence": "high",
                                "data_source": "database_lookup"},
        }
    ],
    "literature": {"query": "Caffeine AND Adenosine receptor A2a",
                   "pubmed_hits": 31, "top_references": []},
}


@pytest.fixture
def mock_pipeline(monkeypatch):
    """Replace run_pipeline with one that returns a fixed Readout."""
    def fake_run(*_args, **kwargs):
        # Honor the on_progress callback if provided, like the real one.
        cb = kwargs.get("on_progress")
        if cb:
            cb("Validating SMILES...")
        return readout.Readout(SAMPLE_DATA)

    monkeypatch.setattr(cli.pipeline, "run_pipeline", fake_run)


def test_missing_required_args_fails():
    result = CliRunner().invoke(cli.main, [])
    assert result.exit_code != 0  # Click reports the missing -s/-f


def test_json_output_to_stdout(mock_pipeline):
    result = CliRunner().invoke(cli.main, ["-s", "CN1...", "-f", "P29274", "--format", "json"])
    assert result.exit_code == 0
    assert '"cid": 2519' in result.output or '"cid":2519' in result.output
    assert "Caffeine" in result.output


def test_pretty_output_to_stdout(mock_pipeline):
    # Give Rich a wide "terminal" so table cells aren't truncated.
    result = CliRunner().invoke(
        cli.main, ["-s", "CN1...", "-f", "P29274"], env={"COLUMNS": "200"}
    )
    assert result.exit_code == 0
    assert "Compound" in result.output
    assert "P29274" in result.output
    assert "MODERATE" in result.output


def test_output_to_file(mock_pipeline, tmp_path):
    out_file = tmp_path / "result.json"
    result = CliRunner().invoke(
        cli.main, ["-s", "CN1...", "-f", "P29274", "-o", str(out_file)]
    )
    assert result.exit_code == 0
    assert out_file.exists()
    # File defaults to JSON and should parse back cleanly.
    data = json.loads(out_file.read_text(encoding="utf-8"))
    assert data["compound"]["cid"] == 2519


def test_fasta_file_is_read(monkeypatch, tmp_path):
    """A -f file path should be read and its contents passed to the pipeline."""
    fasta_file = tmp_path / "p.fasta"
    fasta_file.write_text(">sp|P29274|AA2AR_HUMAN\nMPIMGSSVYA", encoding="utf-8")

    captured = {}

    def fake_run(smiles, fasta_inputs, **kwargs):
        captured["fasta_inputs"] = fasta_inputs
        return readout.Readout(SAMPLE_DATA)

    monkeypatch.setattr(cli.pipeline, "run_pipeline", fake_run)

    result = CliRunner().invoke(cli.main, ["-s", "CN1...", "-f", str(fasta_file)])
    assert result.exit_code == 0
    assert captured["fasta_inputs"][0].startswith(">sp|P29274|")


def test_pipeline_error_exits_nonzero(monkeypatch):
    def boom(*_a, **_k):
        raise cli.pipeline.PipelineError("bad SMILES")
    monkeypatch.setattr(cli.pipeline, "run_pipeline", boom)

    result = CliRunner().invoke(cli.main, ["-s", "junk", "-f", "P29274"])
    assert result.exit_code == 1
