"""
matcher.py
The compound-protein matcher: given a compound's BioAssay results and the
protein target(s) the user actually cares about, it filters down to just
the assays that tested the compound against each target and summarizes
the interaction.

This is the step that makes the tool answer a specific question:
"how does THIS compound act on THIS protein?" — instead of dumping every
assay the compound ever appeared in. No network calls; it processes data
from Phase 2a / 3a.
"""

# Work whether imported as a package, run with -m, or run directly.
try:
    from toxreadout import bioassay_parser as bp
except ModuleNotFoundError:
    import bioassay_parser as bp


def normalize_accession(accession: str | None) -> str:
    """
    Put an accession into a comparable form: drop any version suffix
    (P29274.1 -> P29274), trim whitespace, and uppercase it. Returns ""
    for missing values.
    """
    if not accession:
        return ""
    return accession.strip().split(".")[0].upper()


def match_protein(parsed_results: list[dict], accession: str) -> dict:
    """
    Find the assays in `parsed_results` that tested against `accession`
    and summarize them. Returns a per-protein result dict; if the compound
    was never tested against this protein, match_found is False.
    """
    target = normalize_accession(accession)

    matched = [
        row
        for row in parsed_results
        if normalize_accession(row.get("target_accession")) == target
    ]

    if not matched:
        return {
            "protein_accession": accession,
            "match_found": False,
            "total_assays": 0,
            "active_assays": 0,
            "active_ratio": 0.0,
            "potency_range_um": None,
            "assay_types": [],
        }

    summary = bp.get_activity_summary(matched)
    return {
        "protein_accession": accession,
        "match_found": True,
        "total_assays": summary["total_assays"],
        "active_assays": summary["active"],
        "active_ratio": summary["active_ratio"],
        "potency_range_um": summary["potency_range_um"],
        "assay_types": summary["assay_types"],
    }


def match_proteins(raw_results: list[dict], accessions: list[str]) -> list[dict]:
    """
    Match a compound's raw assay rows (from Phase 2a) against one or more
    user-provided protein accessions. Returns one result per protein.
    """
    parsed = bp.parse_results(raw_results)
    return [match_protein(parsed, accession) for accession in accessions]


# Runs only when executed directly: live caffeine matched against the
# adenosine A2A receptor (expected hit) and insulin (expected no data).
if __name__ == "__main__":
    try:
        from toxreadout import pubchem_client
    except ModuleNotFoundError:
        import pubchem_client

    raw = pubchem_client.get_bioassays(2519)  # caffeine
    targets = ["P29274", "P01308"]  # adenosine A2A receptor, insulin

    for result in match_proteins(raw, targets):
        print(f"Protein {result['protein_accession']}:")
        if not result["match_found"]:
            print("  No BioAssay data - caffeine was never tested on this target.\n")
            continue
        potency = result["potency_range_um"]
        potency_str = (
            f"{potency['min']}-{potency['max']} uM" if potency else "no potency values"
        )
        print(f"  Match found: {result['total_assays']} assays, "
              f"{result['active_assays']} active "
              f"(ratio {result['active_ratio']})")
        print(f"  Potency range: {potency_str}")
        print(f"  Assay types:   {', '.join(result['assay_types'])}\n")
