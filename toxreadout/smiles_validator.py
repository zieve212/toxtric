"""
smiles_validator.py
Validates a SMILES string (a text representation of a chemical compound)
and returns its canonical form plus basic molecular properties.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors


def validate_smiles(smiles: str) -> dict:
    """
    Takes a SMILES string and checks whether it describes a real molecule.
    Returns a dictionary with validity, canonical SMILES, formula, and weight.
    Raises ValueError if the SMILES is empty or cannot be parsed.
    """
    # Guard against empty or whitespace-only input
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string is empty.")

    # Ask RDKit to parse the string into a molecule object
    mol = Chem.MolFromSmiles(smiles)

    # RDKit returns None when the string isn't a valid molecule
    if mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    # Build the result dictionary
    return {
        "valid": True,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
    }


# This block runs only when you execute this file directly,
# letting you test it without the rest of the pipeline.
if __name__ == "__main__":
    caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    result = validate_smiles(caffeine)
    print(result)
    