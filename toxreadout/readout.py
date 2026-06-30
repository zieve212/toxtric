"""
readout.py
Assembles every upstream result into one final toxicology readout.

This module adds no analysis of its own — it just gathers the compound
info (Phase 2a), each protein's metadata (2b), the matched interaction
(3b), the risk assessment (4a), and the literature hits (2c) into a
single structured object that can be emitted as JSON or a human-readable
report. The `targets` array holds one entry per protein, so multi-FASTA
inputs are supported.
"""

import json
from datetime import datetime, timezone


# Plain-English interpretation of each concern level. These describe what
# the experimental evidence implies — not any specific molecular mechanism.
_CONCERN_INTERPRETATION = {
    "LOW": "a limited or weak interaction; unlikely to be a primary target "
           "at typical exposures.",
    "MODERATE": "a measurable interaction that may be biologically relevant "
                "and warrants further review.",
    "HIGH": "a strong interaction that is likely biologically significant.",
    "CRITICAL": "a very strong interaction and a high priority for "
                "toxicological concern.",
    "UNKNOWN": "no experimental evidence, so the interaction cannot be "
               "characterized from database data.",
}


def explain_interaction(compound_name: str, target: dict) -> str:
    """
    Build a plain-English explanation of how the compound interacts with one
    protein, derived entirely from the assay evidence and risk score.
    """
    protein = target["protein"]
    interaction = target["interaction"]
    risk = target["risk_assessment"]
    protein_name = protein.get("name") or protein.get("accession")
    concern = risk.get("concern_level", "UNKNOWN")
    interpretation = _CONCERN_INTERPRETATION.get(concern, "")

    if not interaction.get("match_found"):
        return (
            f"{compound_name} has no experimental BioAssay records against "
            f"{protein_name}, so no interaction could be measured. This means "
            f"the databases hold no evidence either way — not that the compound "
            f"is necessarily safe or inactive."
        )

    potency = interaction.get("potency_range_um")
    potency_text = (
        f"with potency as strong as {potency['min']} uM"
        if potency
        else "though no potency value was reported"
    )
    percent = round(interaction.get("active_ratio", 0) * 100)
    return (
        f"{compound_name} was tested against {protein_name} in "
        f"{interaction['total_assays']} assays, of which "
        f"{interaction['active_assays']} showed activity ({percent}%), "
        f"{potency_text}. On the transparent 0-1 risk scale this scores "
        f"{risk['score']} ({concern}) — {interpretation}"
    )


def build_target_entry(protein: dict, interaction: dict, risk: dict) -> dict:
    """
    Combine one protein's metadata, its matched interaction, and its risk
    assessment into a single target entry for the readout.
    """
    return {
        "protein": {
            "accession": protein.get("accession"),
            "name": protein.get("protein_name"),
            "gene": protein.get("gene"),
            "organism": protein.get("organism_scientific"),
        },
        "interaction": {
            "match_found": interaction.get("match_found"),
            "total_assays": interaction.get("total_assays"),
            "active_assays": interaction.get("active_assays"),
            "active_ratio": interaction.get("active_ratio"),
            "potency_range_um": interaction.get("potency_range_um"),
            "assay_types": interaction.get("assay_types"),
        },
        "risk_assessment": {
            "score": risk.get("score"),
            "concern_level": risk.get("concern_level"),
            "reasoning": risk.get("reasoning"),
            "confidence": risk.get("confidence"),
            "data_source": risk.get("data_source"),
        },
    }


class Readout:
    """Holds the assembled readout data and renders it on demand."""

    def __init__(self, data: dict):
        self.data = data

    def to_dict(self) -> dict:
        return self.data

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.data, indent=indent)

    def to_pretty_string(self) -> str:
        """Render a readable, terminal-friendly report."""
        d = self.data
        c = d["compound"]
        lines = [
            "=" * 64,
            "TOXICOLOGY READOUT",
            "=" * 64,
            f"Generated: {d['timestamp']}",
            "",
            "COMPOUND",
            f"  Common name:  {c.get('common_name')}",
            f"  IUPAC name:   {c.get('iupac_name')}",
            f"  Formula:      {c.get('molecular_formula')} "
            f"({c.get('molecular_weight')} g/mol)",
            f"  PubChem CID:  {c.get('cid')}",
            f"  Input SMILES: {c.get('input_smiles')}",
            "",
            f"PROTEIN TARGETS ({len(d['targets'])})",
        ]

        for target in d["targets"]:
            protein = target["protein"]
            interaction = target["interaction"]
            risk = target["risk_assessment"]

            lines.append("-" * 64)
            lines.append(f"  {protein.get('accession')}  {protein.get('name')}")
            lines.append(
                f"    Gene: {protein.get('gene')}   "
                f"Organism: {protein.get('organism')}"
            )

            if not interaction.get("match_found"):
                lines.append("    Interaction: no BioAssay data for this pair.")
                lines.append(f"    Risk: {risk.get('concern_level')}")
                continue

            potency = interaction.get("potency_range_um")
            potency_str = (
                f"{potency['min']}-{potency['max']} uM" if potency else "n/a"
            )
            lines.append(
                f"    Interaction: {interaction['active_assays']}/"
                f"{interaction['total_assays']} active "
                f"(ratio {interaction['active_ratio']}), potency {potency_str}"
            )
            lines.append(
                f"    Risk: {risk['score']} {risk['concern_level']} "
                f"(confidence {risk['confidence']})"
            )
            lines.append(f"    Reasoning: {risk['reasoning']}")

        lines.append("-" * 64)

        # Plain-English interaction explanation per target.
        compound_name = c.get("common_name") or "The compound"
        lines.append("")
        lines.append("INTERACTION EXPLANATION")
        for target in d["targets"]:
            lines.append(f"  {target['protein'].get('accession')}:")
            lines.append(f"    {explain_interaction(compound_name, target)}")

        lit = d.get("literature")
        if lit:
            lines.append("")
            lines.append("LITERATURE")
            lines.append(
                f"  PubMed hits: {lit.get('pubmed_hits')} "
                f'for "{lit.get("query")}"'
            )

            discussion = lit.get("discussion")
            if discussion:
                lines.append("")
                lines.append("  Summary (from the most relevant paper):")
                lines.append(f"    {discussion}")
                source = lit.get("discussion_source")
                if source:
                    lines.append(f"    Source: {source.get('title')}")
                    lines.append(f"            {source.get('url')}")

            references = lit.get("top_references") or []
            if references:
                lines.append("")
                lines.append("  Citations:")
                for i, ref in enumerate(references, 1):
                    lines.append(
                        f"    [{i}] {ref.get('title')} "
                        f"({ref.get('journal')}, {ref.get('year')})"
                    )
                    lines.append(f"        {ref.get('url')}")

        lines.append("=" * 64)
        return "\n".join(lines)


def assemble_readout(
    input_smiles: str,
    compound: dict,
    targets: list[dict],
    literature: dict | None = None,
    timestamp: str | None = None,
) -> Readout:
    """
    Build the complete readout from upstream pieces.

    `targets` is a list of entries from build_target_entry(). `literature`
    is the optional result from literature_client.literature_lookup().
    `timestamp` can be supplied for reproducibility; otherwise the current
    UTC time is used.
    """
    if timestamp is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    data = {
        "timestamp": timestamp,
        "compound": {
            "input_smiles": input_smiles,
            "cid": compound.get("cid"),
            "common_name": compound.get("common_name"),
            "iupac_name": compound.get("iupac_name"),
            "molecular_formula": compound.get("molecular_formula"),
            "molecular_weight": compound.get("molecular_weight"),
        },
        "targets": targets,
        "literature": None,
    }

    if literature is not None:
        data["literature"] = {
            "query": literature.get("query"),
            "pubmed_hits": literature.get("total_results"),
            "discussion": literature.get("discussion"),
            "discussion_source": literature.get("discussion_source"),
            "top_references": literature.get("top_references", []),
        }

    return Readout(data)


# Runs only when executed directly: a full live readout for caffeine
# against the adenosine A2A receptor (a hit) and insulin (no data).
if __name__ == "__main__":
    import sys

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    try:
        from toxreadout import (
            pubchem_client,
            ncbi_client,
            matcher,
            risk_classifier,
            literature_client,
        )
    except ModuleNotFoundError:
        import pubchem_client
        import ncbi_client
        import matcher
        import risk_classifier
        import literature_client

    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # caffeine
    accessions = ["P29274", "P01308"]  # A2A receptor (hit), insulin (no data)

    compound = pubchem_client.fetch_compound(smiles)
    interactions = matcher.match_proteins(compound["bioassay_results"], accessions)

    targets = []
    for accession, interaction in zip(accessions, interactions):
        protein = ncbi_client.fetch_protein(accession)
        risk = risk_classifier.classify_risk(interaction)
        targets.append(build_target_entry(protein, interaction, risk))

    first_protein = ncbi_client.fetch_protein(accessions[0])
    literature = literature_client.literature_lookup(
        compound.get("common_name") or "caffeine",
        first_protein.get("protein_name") or accessions[0],
    )

    readout = assemble_readout(smiles, compound, targets, literature)
    print(readout.to_pretty_string())
