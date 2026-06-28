"""
bioassay_parser.py
Turns the raw BioAssay rows from the PubChem client into clean, grouped,
summarized activity data. This module makes no network calls — it only
processes data already fetched in Phase 2a.

The key jobs:
  - normalize messy outcome labels into a fixed set of categories
  - group assay results by the protein target they tested against
  - summarize activity (totals, active ratio, potency range, assay types)
"""

# The fixed vocabulary we normalize every outcome into.
VALID_OUTCOMES = {"ACTIVE", "INACTIVE", "INCONCLUSIVE", "PROBE", "UNSPECIFIED"}


def normalize_outcome(raw_outcome: str | None) -> str:
    """
    Map PubChem's outcome text onto our fixed enum. Anything unrecognized
    (including blanks) becomes "UNSPECIFIED".
    """
    if not raw_outcome:
        return "UNSPECIFIED"
    key = raw_outcome.strip().upper()
    return key if key in VALID_OUTCOMES else "UNSPECIFIED"


def parse_results(raw_results: list[dict]) -> list[dict]:
    """
    Clean each raw assay row: normalize its outcome and keep just the
    fields the cross-referencing engine needs.
    """
    parsed = []
    for row in raw_results:
        parsed.append(
            {
                "aid": row.get("aid"),
                "outcome": normalize_outcome(row.get("activity_outcome")),
                "target_accession": row.get("target_accession"),
                "activity_value_um": row.get("activity_value_um"),
                "activity_name": row.get("activity_name"),
                "assay_type": row.get("assay_type"),
            }
        )
    return parsed


def group_by_target(parsed_results: list[dict]) -> dict[str, list[dict]]:
    """
    Group parsed results by protein target accession. Assays with no
    protein target (cell/phenotypic screens) are skipped.
    """
    groups: dict[str, list[dict]] = {}
    for row in parsed_results:
        target = row.get("target_accession")
        if not target:
            continue
        groups.setdefault(target, []).append(row)
    return groups


def get_activity_summary(parsed_results: list[dict]) -> dict:
    """
    Summarize a list of parsed assay results: outcome counts, the active
    ratio, the potency range (in uM, over active assays that report a
    value), and the distinct assay types and potency measures seen.
    """
    total = len(parsed_results)

    # Tally each outcome category.
    counts = {outcome: 0 for outcome in VALID_OUTCOMES}
    for row in parsed_results:
        counts[row["outcome"]] += 1

    active = counts["ACTIVE"]
    active_ratio = round(active / total, 2) if total else 0.0

    # Potency range looks only at active assays that reported a value.
    potencies = [
        row["activity_value_um"]
        for row in parsed_results
        if row["outcome"] == "ACTIVE" and row["activity_value_um"] is not None
    ]
    potency_range = (
        {"min": min(potencies), "max": max(potencies)} if potencies else None
    )

    assay_types = sorted({r["assay_type"] for r in parsed_results if r["assay_type"]})
    measures = sorted({r["activity_name"] for r in parsed_results if r["activity_name"]})

    return {
        "total_assays": total,
        "active": active,
        "inactive": counts["INACTIVE"],
        "inconclusive": counts["INCONCLUSIVE"],
        "probe": counts["PROBE"],
        "unspecified": counts["UNSPECIFIED"],
        "active_ratio": active_ratio,
        "potency_range_um": potency_range,
        "assay_types": assay_types,
        "potency_measures": measures,
    }


def summarize_by_target(raw_results: list[dict]) -> dict[str, dict]:
    """
    Convenience pipeline: parse raw rows, group by target, and produce an
    activity summary for each protein target. Returns {accession: summary}.
    """
    parsed = parse_results(raw_results)
    return {
        target: get_activity_summary(rows)
        for target, rows in group_by_target(parsed).items()
    }


# Runs only when you execute this file directly: a live look at caffeine's
# real BioAssay data, grouped and summarized by protein target.
if __name__ == "__main__":
    # Work whether run directly (python toxreadout/bioassay_parser.py) or
    # as a module (python -m toxreadout.bioassay_parser).
    try:
        from toxreadout import pubchem_client
    except ModuleNotFoundError:
        import pubchem_client

    raw = pubchem_client.get_bioassays(2519)  # caffeine
    parsed = parse_results(raw)

    overall = get_activity_summary(parsed)
    print(f"Caffeine: {overall['total_assays']} total BioAssay rows")
    print(f"  Active: {overall['active']}   Inactive: {overall['inactive']}   "
          f"Inconclusive: {overall['inconclusive']}   "
          f"Unspecified: {overall['unspecified']}")

    by_target = summarize_by_target(raw)
    print(f"\nAssays mapped to {len(by_target)} distinct protein targets.")

    # Show the targets with the most assays.
    top = sorted(by_target.items(), key=lambda kv: kv[1]["total_assays"], reverse=True)
    print("\nTop protein targets by assay count:")
    for accession, summary in top[:5]:
        potency = summary["potency_range_um"]
        potency_str = (
            f"{potency['min']}-{potency['max']} uM" if potency else "no potency values"
        )
        print(f"  {accession}: {summary['total_assays']} assays, "
              f"{summary['active']} active (ratio {summary['active_ratio']}), "
              f"{potency_str}")
