"""
pubchem_client.py
Talks to PubChem's public PUG REST API to turn a SMILES string into a
compound record: its CID, basic properties, and BioAssay results.

This is the first module that makes live calls over the internet, so it
includes polite rate limiting (max 5 requests/second) and automatic
retries with exponential backoff when PubChem's servers hiccup.
"""

import time

import requests

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# PubChem asks callers to stay under 5 requests per second.
MIN_INTERVAL = 0.2  # seconds between requests -> at most 5/sec
MAX_RETRIES = 3
RETRY_STATUS = (500, 503)  # server-side errors worth retrying

# Tracks when we last hit the API, so we can space requests out.
_last_request_time = 0.0


def _throttle() -> None:
    """Pause if needed so we never exceed the rate limit."""
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_request_time = time.time()


def _request(url: str, method: str = "GET", data: dict | None = None):
    """
    Make a rate-limited HTTP request with retry/backoff on server errors.
    Returns the requests.Response object.
    """
    last_response = None
    for attempt in range(MAX_RETRIES):
        _throttle()
        if method == "POST":
            last_response = requests.post(url, data=data, timeout=30)
        else:
            last_response = requests.get(url, timeout=30)

        # Retry only on transient server errors; wait 1s, 2s, 4s...
        if last_response.status_code in RETRY_STATUS:
            time.sleep(2 ** attempt)
            continue
        return last_response

    return last_response


def get_cid(smiles: str) -> int:
    """
    Resolve a SMILES string to its PubChem Compound ID (CID).

    Uses POST so special characters in the SMILES don't break the URL.
    Raises ValueError if PubChem has no compound for this SMILES.
    """
    url = f"{BASE_URL}/compound/smiles/cids/JSON"
    response = _request(url, method="POST", data={"smiles": smiles})

    if response.status_code == 404:
        raise ValueError(f"No PubChem compound found for SMILES: '{smiles}'")
    response.raise_for_status()

    cids = response.json()["IdentifierList"]["CID"]
    return cids[0]


def get_properties(cid: int) -> dict:
    """Fetch molecular formula, weight, and IUPAC name for a CID."""
    props = "MolecularFormula,MolecularWeight,IUPACName"
    url = f"{BASE_URL}/compound/cid/{cid}/property/{props}/JSON"
    response = _request(url)
    response.raise_for_status()

    row = response.json()["PropertyTable"]["Properties"][0]
    return {
        "molecular_formula": row.get("MolecularFormula"),
        # PubChem returns weight as a string; convert to a number.
        "molecular_weight": float(row["MolecularWeight"]),
        "iupac_name": row.get("IUPACName"),
    }


def get_bioassays(cid: int) -> list[dict]:
    """
    Fetch all BioAssay summary rows for a CID.

    Returns a list of assay records, each with the fields the build guide
    needs. Returns an empty list if the compound has no assay data (404).
    """
    url = f"{BASE_URL}/compound/cid/{cid}/assaysummary/JSON"
    response = _request(url)

    if response.status_code == 404:
        return []
    response.raise_for_status()

    table = response.json().get("Table", {})
    columns = table.get("Columns", {}).get("Column", [])
    rows = table.get("Row", [])

    results = []
    for row in rows:
        # Pair each column name with its cell value for easy lookup.
        record = dict(zip(columns, row.get("Cell", [])))
        results.append(
            {
                "aid": record.get("AID"),
                "activity_outcome": record.get("Activity Outcome"),
                "target_gi": record.get("Target GI"),
                "target_name": record.get("Target Name"),
                "target_accession": record.get("Target Accession"),
            }
        )
    return results


def fetch_compound(smiles: str) -> dict:
    """
    Full lookup for one compound: SMILES -> CID -> properties + assays.
    Returns the record shaped exactly like the build guide's schema.
    """
    cid = get_cid(smiles)
    properties = get_properties(cid)
    bioassays = get_bioassays(cid)

    return {
        "cid": cid,
        "iupac_name": properties["iupac_name"],
        "molecular_formula": properties["molecular_formula"],
        "molecular_weight": properties["molecular_weight"],
        "bioassay_count": len(bioassays),
        "bioassay_results": bioassays,
    }


# Runs only when you execute this file directly: a live caffeine lookup.
if __name__ == "__main__":
    caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    record = fetch_compound(caffeine)
    print(f"CID:               {record['cid']}")
    print(f"IUPAC name:        {record['iupac_name']}")
    print(f"Molecular formula: {record['molecular_formula']}")
    print(f"Molecular weight:  {record['molecular_weight']}")
    print(f"BioAssay rows:     {record['bioassay_count']}")
