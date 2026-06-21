"""
fasta_validator.py
Validates one or more protein sequences supplied in FASTA format and
returns clean, structured information about each one.
"""

import io

from Bio import SeqIO

# The 20 standard amino-acid single-letter codes, plus X for "unknown residue".
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYX")


def _extract_accession(header: str) -> str | None:
    """
    Pull an accession ID out of a FASTA header if one is present.

    Recognizes two common formats:
      UniProt:  sp|P29274|AA2AR_HUMAN   -> "P29274"
      NCBI:     gi|12345                -> "12345"
    Returns None if neither pattern is found.
    """
    if not header:
        return None

    # The accession lives in the first whitespace-delimited token.
    first_token = header.split()[0]
    parts = first_token.split("|")

    # UniProt style: db|ACCESSION|ENTRY_NAME -> middle field
    if parts[0] in ("sp", "tr") and len(parts) >= 2:
        return parts[1]

    # NCBI style: gi|ACCESSION -> field right after "gi"
    if parts[0] == "gi" and len(parts) >= 2:
        return parts[1]

    return None


def validate_fasta(fasta: str) -> list[dict]:
    """
    Takes a FASTA string (one or more sequences) and validates each one.

    Returns a list of dictionaries, one per protein, each containing:
        valid, header, sequence, length, accession

    Raises ValueError if the input is empty or contains no parseable
    sequence, or if any sequence contains invalid amino-acid characters.
    """
    # Guard against empty or whitespace-only input
    if not fasta or not fasta.strip():
        raise ValueError("FASTA input is empty.")

    text = fasta.strip()

    # The guide allows a bare sequence with no ">header" line.
    # BioPython needs a header to parse, so add a placeholder one.
    if not text.startswith(">"):
        text = ">query\n" + text

    # SeqIO reads from a file-like object, so wrap the string in StringIO.
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))

    if not records:
        raise ValueError("No valid FASTA sequence found in input.")

    results = []
    for record in records:
        sequence = str(record.seq).upper()

        # Reject any residue that isn't a recognized amino-acid code.
        invalid_chars = set(sequence) - VALID_AMINO_ACIDS
        if invalid_chars:
            bad = ", ".join(sorted(invalid_chars))
            raise ValueError(
                f"Invalid amino-acid character(s) in sequence: {bad}"
            )

        header = record.description
        results.append(
            {
                "valid": True,
                "header": header,
                "sequence": sequence,
                "length": len(sequence),
                "accession": _extract_accession(header),
            }
        )

    return results


# Runs only when you execute this file directly, so you can see it work
# without the rest of the pipeline.
if __name__ == "__main__":
    sample = ">sp|P29274|AA2AR_HUMAN Adenosine receptor A2a\nMPIMGSSVYITVELAIAVLAILGNVLVCWAVWLNSNLQNVTNYFVVSLAAADIAVGVLAIP"
    for protein in validate_fasta(sample):
        print(protein)
