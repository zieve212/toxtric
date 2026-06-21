"""
test_smiles_validator.py
Tests for the SMILES validator module.
Run with: pytest
"""

import pytest
from toxreadout.smiles_validator import validate_smiles


def test_valid_caffeine():
    """A valid SMILES (caffeine) should return valid=True with correct weight."""
    result = validate_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert result["valid"] is True
    assert result["molecular_formula"] == "C8H10N4O2"
    assert round(result["molecular_weight"]) == 194


def test_valid_red40():
    """A second valid SMILES (Red 40 dye) should parse without error."""
    red40 = "OC1=CC2=CC(=CC=C2C(=C1)N=NC1=CC=C(C=C1S(O)(=O)=O)S(O)(=O)=O)S(O)(=O)=O"
    result = validate_smiles(red40)
    assert result["valid"] is True


def test_invalid_string():
    """Gibberish text should raise a ValueError."""
    with pytest.raises(ValueError):
        validate_smiles("not_a_smiles")


def test_empty_string():
    """An empty string should raise a ValueError."""
    with pytest.raises(ValueError):
        validate_smiles("")
        