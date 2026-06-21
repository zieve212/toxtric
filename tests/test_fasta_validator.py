"""
test_fasta_validator.py
Tests for the FASTA validator module.
Run with: pytest
"""

import pytest
from toxreadout.fasta_validator import validate_fasta


def test_valid_single_protein():
    """A single valid FASTA record should parse with the right fields."""
    fasta = ">sp|P29274|AA2AR_HUMAN Adenosine receptor A2a\nMPIMGSSVYITVELAIA"
    result = validate_fasta(fasta)
    assert len(result) == 1
    assert result[0]["valid"] is True
    assert result[0]["accession"] == "P29274"
    assert result[0]["sequence"] == "MPIMGSSVYITVELAIA"
    assert result[0]["length"] == 17


def test_valid_multi_protein():
    """Two sequences in one input should both come back (protein complex)."""
    fasta = (
        ">sp|P29274|AA2AR_HUMAN\nMPIMGSSVYA\n"
        ">gi|12345 some other protein\nGNVLVCWAVW"
    )
    result = validate_fasta(fasta)
    assert len(result) == 2
    assert result[0]["accession"] == "P29274"
    assert result[1]["accession"] == "12345"


def test_bare_sequence_without_header():
    """A raw sequence with no >header line should still validate."""
    result = validate_fasta("MPIMGSSVYA")
    assert len(result) == 1
    assert result[0]["valid"] is True
    assert result[0]["accession"] is None


def test_invalid_sequence_with_numbers():
    """A sequence containing digits should raise a ValueError."""
    with pytest.raises(ValueError):
        validate_fasta(">bad\nMPIM1234GS")


def test_empty_input():
    """An empty string should raise a ValueError."""
    with pytest.raises(ValueError):
        validate_fasta("")
