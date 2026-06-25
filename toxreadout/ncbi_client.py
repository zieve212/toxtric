"""
ncbi_client.py
Looks up rich protein metadata for a given accession number.

Most fields (official name, gene symbol, organism, sequence length,
molecular weight, and a plain-English functional description) come from
the UniProt REST API, which carries this data cleanly. The NCBI GI
number is fetched separately from NCBI E-utilities, since UniProt no
longer reports it and Phase 3 needs the GI to cross-reference BioAssay
targets.

Like the PubChem client, this makes live calls, so it rate-limits and
retries transient server errors.
"""

import time

import requests

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI allows 3 requests/second without an API key; stay just under that.
MIN_INTERVAL = 0.34  # seconds between requests
MAX_RETRIES = 3
RETRY_STATUS = (500, 502, 503)

_last_request_time = 0.0


def _throttle() -> None:
    """Pause if needed so we stay under the rate limit."""
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_request_time = time.time()


def _request(url: str, params: dict | None = None):
    """Rate-limited GET with retry/backoff on transient server errors."""
    last_response = None
    for attempt in range(MAX_RETRIES):
        _throttle()
        last_response = requests.get(url, params=params, timeout=30)
        if last_response.status_code in RETRY_STATUS:
            time.sleep(2 ** attempt)
            continue
        return last_response
    return last_response


def _parse_uniprot(data: dict) -> dict:
    """Pull the fields we care about out of a UniProtKB JSON record."""
    # Official protein name (fall back to a submission name if needed).
    description = data.get("proteinDescription", {})
    recommended = description.get("recommendedName") or {}
    name = recommended.get("fullName", {}).get("value")
    if not name:
        submissions = description.get("submissionNames") or []
        if submissions:
            name = submissions[0].get("fullName", {}).get("value")

    # Gene symbol (e.g. "INS").
    genes = data.get("genes") or []
    gene = genes[0].get("geneName", {}).get("value") if genes else None

    # Organism: scientific name, English (common) name, taxon id, lineage.
    organism = data.get("organism", {})

    # Sequence length and molecular weight (UniProt reports weight in Daltons).
    sequence = data.get("sequence", {})

    # First FUNCTION comment is the plain-English description of what it does.
    function = None
    for comment in data.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts") or []
            if texts:
                function = texts[0].get("value")
            break

    return {
        "accession": data.get("primaryAccession"),
        "protein_name": name,
        "gene": gene,
        "organism_scientific": organism.get("scientificName"),
        "organism_common": organism.get("commonName"),
        "taxon_id": organism.get("taxonId"),
        "lineage": organism.get("lineage", []),
        "length": sequence.get("length"),
        "molecular_weight": sequence.get("molWeight"),
        "function": function,
    }


def get_uniprot_record(accession: str) -> dict:
    """Fetch and parse a UniProtKB record for an accession."""
    url = f"{UNIPROT_URL}/{accession}.json"
    response = _request(url)

    if response.status_code in (400, 404):
        raise ValueError(f"No UniProt record found for accession: '{accession}'")
    response.raise_for_status()

    return _parse_uniprot(response.json())


def get_gi_number(accession: str) -> str | None:
    """
    Look up the NCBI GI number for an accession via E-utilities esearch.
    In the protein database, the returned UID is the GI number.
    Returns None if no match is found.
    """
    url = f"{NCBI_EUTILS}/esearch.fcgi"
    params = {"db": "protein", "term": accession, "retmode": "json"}
    response = _request(url, params=params)
    response.raise_for_status()

    id_list = response.json().get("esearchresult", {}).get("idlist", [])
    return id_list[0] if id_list else None


def fetch_protein(accession: str) -> dict:
    """
    Full protein lookup: UniProt metadata plus the NCBI GI number.
    Returns a single record describing the protein target.
    """
    record = get_uniprot_record(accession)
    record["gi_number"] = get_gi_number(accession)
    return record


# Runs only when you execute this file directly: a live lookup of a small
# gene (insulin, gene INS, accession P01308).
if __name__ == "__main__":
    record = fetch_protein("P01308")

    weight = record["molecular_weight"]
    weight_str = f"{weight:,} Da ({weight / 1000:.2f} kDa)" if weight else "n/a"
    lineage = " > ".join(record["lineage"][:5]) if record["lineage"] else "n/a"

    print(f"Accession:         {record['accession']}")
    print(f"GI number:         {record['gi_number']}")
    print(f"Protein name:      {record['protein_name']}")
    print(f"Gene symbol:       {record['gene']}")
    print(f"Organism:          {record['organism_scientific']} "
          f"({record['organism_common']}), taxon {record['taxon_id']}")
    print(f"Taxonomy:          {lineage} ...")
    print(f"Sequence length:   {record['length']} residues")
    print(f"Molecular weight:  {weight_str}")
    print(f"Function:          {record['function']}")
