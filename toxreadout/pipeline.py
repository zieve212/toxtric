"""
pipeline.py
Orchestrates the whole ToxReadout pipeline end to end.

It wires the modules together in order:
    validate SMILES -> validate FASTA -> fetch compound (PubChem)
    -> match assays to each protein -> fetch protein metadata + score risk
    -> search literature -> assemble the final readout.

Progress is reported through an optional callback so callers (like the
CLI) can display it however they like. Fatal problems raise PipelineError
with a clear message; non-fatal hiccups (e.g. one protein's metadata
lookup) are logged and skipped so the rest of the run still completes.
"""

# Work whether imported as a package, run with -m, or run directly.
try:
    from toxreadout import (
        smiles_validator,
        fasta_validator,
        pubchem_client,
        ncbi_client,
        matcher,
        risk_classifier,
        literature_client,
        readout,
    )
except ModuleNotFoundError:
    import smiles_validator
    import fasta_validator
    import pubchem_client
    import ncbi_client
    import matcher
    import risk_classifier
    import literature_client
    import readout


class PipelineError(Exception):
    """Raised when a stage fails in a way that stops the whole run."""


def run_pipeline(
    smiles: str,
    fasta_inputs: list[str],
    include_literature: bool = True,
    on_progress=None,
    timestamp: str | None = None,
) -> readout.Readout:
    """
    Run the full compound-vs-protein analysis and return a Readout.

    Args:
        smiles: the compound's SMILES string.
        fasta_inputs: one or more FASTA strings (protein targets).
        include_literature: whether to run the PubMed search.
        on_progress: optional callable(str) for progress messages.
        timestamp: optional fixed timestamp for reproducible output.
    """
    log = on_progress or (lambda _message: None)

    # --- Stage 1: validate the SMILES (offline) ---------------------------
    log("Validating SMILES...")
    try:
        smiles_info = smiles_validator.validate_smiles(smiles)
    except ValueError as exc:
        raise PipelineError(f"SMILES validation failed: {exc}") from exc
    canonical_smiles = smiles_info["canonical_smiles"]

    # --- Stage 2: validate the FASTA input(s) (offline) -------------------
    log("Validating protein FASTA input...")
    proteins = []
    for fasta in fasta_inputs:
        try:
            proteins.extend(fasta_validator.validate_fasta(fasta))
        except ValueError as exc:
            raise PipelineError(f"FASTA validation failed: {exc}") from exc

    accessions = [p["accession"] for p in proteins if p.get("accession")]
    if not accessions:
        raise PipelineError(
            "No protein accession found in the FASTA input(s). Each header "
            "must include a UniProt (sp|XXXXX|NAME) or NCBI (gi|XXXXX) accession."
        )

    # --- Stage 3: fetch compound data from PubChem (live) -----------------
    log("Resolving PubChem CID and fetching compound data...")
    try:
        compound = pubchem_client.fetch_compound(canonical_smiles)
    except Exception as exc:  # network/HTTP/parse problems are all fatal here
        raise PipelineError(f"PubChem lookup failed: {exc}") from exc
    log(f"Fetched {compound['bioassay_count']} BioAssay rows (CID {compound['cid']}).")

    # --- Stage 4: match assays to targets, fetch metadata, score risk -----
    log("Matching compound activity to protein targets...")
    interactions = matcher.match_proteins(compound["bioassay_results"], accessions)

    targets = []
    for accession, interaction in zip(accessions, interactions):
        log(f"  {accession}: fetching protein metadata and scoring risk...")
        try:
            protein = ncbi_client.fetch_protein(accession)
        except Exception as exc:
            # One protein failing shouldn't sink the whole run.
            log(f"  Warning: metadata lookup for {accession} failed ({exc}); "
                "continuing with accession only.")
            protein = {"accession": accession}

        risk = risk_classifier.classify_risk(interaction)
        targets.append(readout.build_target_entry(protein, interaction, risk))

    # --- Stage 5: literature search (live, optional, non-fatal) -----------
    literature = None
    if include_literature:
        log("Searching PubMed literature...")
        compound_name = compound.get("common_name") or "compound"
        protein_name = next(
            (t["protein"]["name"] for t in targets if t["protein"].get("name")),
            accessions[0],
        )
        try:
            literature = literature_client.literature_lookup(compound_name, protein_name)
        except Exception as exc:
            log(f"  Warning: literature lookup failed ({exc}); skipping.")
            literature = None

    # --- Stage 6: assemble the final readout ------------------------------
    log("Assembling readout...")
    result = readout.assemble_readout(
        smiles, compound, targets, literature, timestamp=timestamp
    )
    log("Done.")
    return result


# Runs only when executed directly: a full live run on caffeine against
# two protein targets, with progress printed to the screen.
if __name__ == "__main__":
    import sys

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    fasta_a2a = ">sp|P29274|AA2AR_HUMAN Adenosine receptor A2a\nMPIMGSSVYITVELAIA"
    fasta_insulin = ">sp|P01308|INS_HUMAN Insulin\nMALWMRLLPLLALLALWGPD"

    result = run_pipeline(
        caffeine,
        [fasta_a2a, fasta_insulin],
        on_progress=lambda msg: print(f"[pipeline] {msg}"),
    )

    print()
    print(result.to_pretty_string())
