"""
resolvers.py
Convenience lookups that turn plain-English names into the precise inputs
the pipeline needs, so a user can type "caffeine" and "adenosine A2A
receptor" instead of a SMILES string and a FASTA sequence.

    resolve_compound(name) -> SMILES + CID + properties (via PubChem)
    resolve_protein(name)  -> FASTA + accession + metadata (via UniProt)

Each resolver also reports the English name and the resolved code, so the
caller can show the user exactly what was picked (and let them copy/edit
the SMILES or FASTA for precise re-runs). Name lookup is convenient but
fuzzy: a vague name may match the wrong entry, so always surface the
choice. Reuses the rate-limited request helpers from the API clients.
"""

from urllib.parse import quote

# Work whether imported as a package, run with -m, or run directly.
try:
    from toxreadout import pubchem_client, ncbi_client
except ModuleNotFoundError:
    import pubchem_client
    import ncbi_client

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def _name_to_cid(name: str) -> int:
    """Resolve a compound name to its PubChem CID."""
    url = f"{pubchem_client.BASE_URL}/compound/name/{quote(name)}/cids/JSON"
    response = pubchem_client._request(url)
    if response.status_code == 404:
        raise ValueError(f"No PubChem compound found for name: '{name}'")
    response.raise_for_status()
    return response.json()["IdentifierList"]["CID"][0]


def _cid_to_smiles(cid: int) -> str:
    """Fetch a canonical SMILES string for a CID."""
    url = f"{pubchem_client.BASE_URL}/compound/cid/{cid}/property/SMILES/JSON"
    response = pubchem_client._request(url)
    response.raise_for_status()
    props = response.json()["PropertyTable"]["Properties"][0]
    # PubChem returns the structure under "SMILES" (older "CanonicalSMILES"
    # now comes back as "ConnectivitySMILES").
    return props.get("SMILES") or props.get("ConnectivitySMILES")


def resolve_compound(name: str) -> dict:
    """
    Turn a compound name into the SMILES and identity info the pipeline
    needs. Raises ValueError if the name can't be resolved.
    """
    cid = _name_to_cid(name)
    smiles = _cid_to_smiles(cid)
    properties = pubchem_client.get_properties(cid)
    common_name = pubchem_client.get_common_name(cid)

    return {
        "query": name,
        "cid": cid,
        "common_name": common_name,
        "smiles": smiles,
        "iupac_name": properties.get("iupac_name"),
        "molecular_formula": properties.get("molecular_formula"),
        "molecular_weight": properties.get("molecular_weight"),
    }


def resolve_protein(name: str, organism_id: int = 9606, reviewed: bool = True) -> dict:
    """
    Turn a protein name into a FASTA and identity info. Defaults to the
    reviewed human entry; pass organism_id=None to search all organisms.
    Raises ValueError if nothing matches.
    """
    query_parts = [name]
    if organism_id:
        query_parts.append(f"organism_id:{organism_id}")
    if reviewed:
        query_parts.append("reviewed:true")

    params = {
        "query": " AND ".join(query_parts),
        "format": "json",
        "size": "1",
        "fields": "accession,id,protein_name,organism_name,gene_primary,sequence",
    }
    response = ncbi_client._request(UNIPROT_SEARCH_URL, params=params)
    response.raise_for_status()

    results = response.json().get("results", [])
    if not results:
        raise ValueError(f"No UniProt protein found for name: '{name}'")

    record = results[0]
    accession = record.get("primaryAccession")
    entry_name = record.get("uniProtkbId")
    protein_name = (
        record.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
    )
    organism = record.get("organism", {}).get("scientificName")
    genes = record.get("genes") or []
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    sequence = record.get("sequence", {}).get("value", "")

    fasta = f">sp|{accession}|{entry_name} {protein_name}\n{sequence}"

    return {
        "query": name,
        "accession": accession,
        "entry_name": entry_name,
        "protein_name": protein_name,
        "organism": organism,
        "gene": gene,
        "sequence": sequence,
        "fasta": fasta,
    }


# Runs only when executed directly: resolve a name pair and show the codes.
if __name__ == "__main__":
    compound = resolve_compound("caffeine")
    print(f"Compound: {compound['common_name']} (CID {compound['cid']})")
    print(f"  SMILES: {compound['smiles']}")
    print(f"  Formula: {compound['molecular_formula']} "
          f"({compound['molecular_weight']} g/mol)\n")

    protein = resolve_protein("adenosine A2A receptor")
    print(f"Protein: {protein['protein_name']} "
          f"({protein['accession']}, {protein['organism']})")
    print(f"  Gene: {protein['gene']}")
    print(f"  FASTA header: >sp|{protein['accession']}|{protein['entry_name']}")
    print(f"  Sequence: {protein['sequence'][:40]}... ({len(protein['sequence'])} aa)")
