"""
cli.py
Command-line interface for ToxReadout.

Usage:
    toxreadout -s "<SMILES>" -f protein.fasta
    toxreadout -s "<SMILES>" -f a.fasta -f b.fasta -o results.json --format json

Arguments are parsed with Click; output is formatted with Rich. Progress
messages are written to stderr so that JSON sent to stdout stays clean
and pipeable.
"""

import os

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

# Work whether imported as a package, run with -m, or run directly.
try:
    from toxreadout import pipeline, resolvers, readout
except ModuleNotFoundError:
    import pipeline
    import resolvers
    import readout

# Results go to stdout; progress/status go to stderr.
console = Console()
err_console = Console(stderr=True)

# Colour per concern level for the Rich table.
CONCERN_COLORS = {
    "LOW": "green",
    "MODERATE": "yellow",
    "HIGH": "dark_orange",
    "CRITICAL": "bold red",
    "UNKNOWN": "dim",
}


def render_rich(data: dict) -> None:
    """Display the readout as a Rich panel + table on stdout."""
    c = data["compound"]
    console.print(
        Panel.fit(
            f"[bold]{c.get('common_name') or '(unnamed compound)'}[/]\n"
            f"IUPAC:   {c.get('iupac_name')}\n"
            f"Formula: {c.get('molecular_formula')}  "
            f"({c.get('molecular_weight')} g/mol)\n"
            f"CID:     {c.get('cid')}\n"
            f"SMILES:  {c.get('input_smiles')}",
            title="Compound",
            border_style="cyan",
        )
    )

    table = Table(title="Protein Targets", show_lines=True)
    table.add_column("Accession", style="bold")
    table.add_column("Protein")
    table.add_column("Gene")
    table.add_column("Assays", justify="right")
    table.add_column("Active", justify="right")
    table.add_column("Potency uM")
    table.add_column("Score", justify="right")
    table.add_column("Concern")

    for target in data["targets"]:
        p = target["protein"]
        inter = target["interaction"]
        risk = target["risk_assessment"]
        concern = risk.get("concern_level", "UNKNOWN")
        color = CONCERN_COLORS.get(concern, "white")

        if not inter.get("match_found"):
            table.add_row(
                p.get("accession") or "-", p.get("name") or "-",
                p.get("gene") or "-", "0", "-", "-", "-",
                f"[{color}]NO DATA[/]",
            )
            continue

        potency = inter.get("potency_range_um")
        potency_str = f"{potency['min']}-{potency['max']}" if potency else "-"
        table.add_row(
            p.get("accession") or "-",
            p.get("name") or "-",
            p.get("gene") or "-",
            str(inter.get("total_assays", 0)),
            str(inter.get("active_assays", 0)),
            potency_str,
            str(risk.get("score")),
            f"[{color}]{concern}[/]",
        )

    console.print(table)

    # Plain-English explanation of how the compound interacts with each target.
    compound_name = data["compound"].get("common_name") or "The compound"
    console.print()
    console.rule("[bold]How they interact[/]")
    for target in data["targets"]:
        console.print(
            f"[bold cyan]{target['protein'].get('accession')}[/]  "
            + readout.explain_interaction(compound_name, target)
        )

    lit = data.get("literature")
    if lit:
        console.print()
        console.rule("[bold]Literature[/]")
        console.print(
            f'[dim]{lit.get("pubmed_hits")} PubMed hits for[/] "{lit.get("query")}"'
        )

        discussion = lit.get("discussion")
        if discussion:
            console.print()
            console.print("[bold]Summary from the most relevant paper:[/]")
            console.print(discussion)
            source = lit.get("discussion_source")
            if source and source.get("url"):
                console.print(
                    f"[dim]Source:[/] [link={source['url']}]{source.get('title')}[/link]"
                )

        references = lit.get("top_references") or []
        if references:
            console.print()
            console.print("[bold]Citations:[/]")
            for i, ref in enumerate(references, 1):
                url = ref.get("url")
                meta = f"({ref.get('journal')}, {ref.get('year')})"
                console.print(f"  [cyan]{i}.[/] {ref.get('title')} [dim]{meta}[/]")
                if url:
                    console.print(f"     [link={url}]{url}[/link]")


def _load_fasta_inputs(fasta_args: tuple) -> list:
    """Each -f value may be a file path or a raw sequence/accession string."""
    inputs = []
    for item in fasta_args:
        if os.path.isfile(item):
            with open(item, encoding="utf-8") as handle:
                inputs.append(handle.read())
        else:
            inputs.append(item)
    return inputs


def _resolve_compound_name(name: str) -> str:
    """Resolve a compound name to a SMILES, showing the user what was picked."""
    with err_console.status(f"[bold cyan]Resolving compound '{name}'..."):
        try:
            compound = resolvers.resolve_compound(name)
        except Exception as exc:
            err_console.print(f"[bold red]Error resolving compound:[/] {exc}")
            raise SystemExit(1)
    err_console.print(
        f"[green]Compound:[/] {compound['common_name']} "
        f"(CID {compound['cid']}, {compound['molecular_formula']})"
    )
    err_console.print(f"  [cyan]SMILES:[/] {compound['smiles']}")
    return compound["smiles"]


def _resolve_protein_name(name: str) -> str:
    """Resolve a protein name to a FASTA, showing the user what was picked."""
    with err_console.status(f"[bold cyan]Resolving protein '{name}'..."):
        try:
            protein = resolvers.resolve_protein(name)
        except Exception as exc:
            err_console.print(f"[bold red]Error resolving protein:[/] {exc}")
            raise SystemExit(1)
    err_console.print(
        f"[green]Protein:[/] {protein['protein_name']} "
        f"({protein['accession']}, {protein['organism']})"
    )
    err_console.print(
        f"  [cyan]FASTA:[/] >sp|{protein['accession']}|{protein['entry_name']} "
        f"({len(protein['sequence'])} aa) - use -f with this accession to re-run precisely"
    )
    return protein["fasta"]


@click.command()
@click.option("--smiles", "-s", default=None,
              help="Compound SMILES string.")
@click.option("--compound-name", "compound_name", default=None,
              help="Compound English name (resolved to SMILES via PubChem).")
@click.option("--fasta", "-f", "fasta", multiple=True,
              help="Protein FASTA file path or raw sequence/accession. "
                   "Repeat -f for multiple targets.")
@click.option("--protein-name", "protein_names", multiple=True,
              help="Protein English name (resolved to FASTA via UniProt). "
                   "Repeat for multiple targets.")
@click.option("--output", "-o", type=click.Path(), default=None,
              help="Write the readout to this file instead of the screen.")
@click.option("--format", "output_format",
              type=click.Choice(["pretty", "json"]), default=None,
              help="Output format. Defaults to pretty on screen, json to file.")
@click.option("--no-literature", is_flag=True, default=False,
              help="Skip the PubMed literature search (faster).")
def main(smiles, compound_name, fasta, protein_names, output, output_format, no_literature):
    """Generate a toxicology readout for a compound against protein target(s).

    Provide the compound as a SMILES (-s) or an English name
    (--compound-name), and the protein(s) as FASTA (-f) or English name(s)
    (--protein-name). Name lookups print the resolved SMILES/FASTA so you
    can copy and tweak them.
    """
    # Scientific abstracts contain Unicode; keep the Windows console happy.
    import sys
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")

    # --- Resolve the compound (name takes the resolver path) ---
    if smiles and compound_name:
        err_console.print("[bold red]Error:[/] use either -s or --compound-name, not both.")
        raise SystemExit(2)
    if compound_name:
        smiles = _resolve_compound_name(compound_name)
    if not smiles:
        err_console.print("[bold red]Error:[/] provide a compound via -s or --compound-name.")
        raise SystemExit(2)

    # --- Resolve the protein target(s) from files, raw input, and names ---
    fasta_inputs = _load_fasta_inputs(fasta)
    for name in protein_names:
        fasta_inputs.append(_resolve_protein_name(name))
    if not fasta_inputs:
        err_console.print("[bold red]Error:[/] provide a protein via -f or --protein-name.")
        raise SystemExit(2)

    # Pretty for the screen, json for files, unless the user overrides.
    if output_format is None:
        output_format = "json" if output else "pretty"

    try:
        with err_console.status("[bold cyan]Starting...") as status:
            result = pipeline.run_pipeline(
                smiles,
                fasta_inputs,
                include_literature=not no_literature,
                on_progress=lambda msg: status.update(f"[bold cyan]{msg}"),
            )
    except pipeline.PipelineError as exc:
        err_console.print(f"[bold red]Error:[/] {exc}")
        raise SystemExit(1)

    # Write to a file...
    if output:
        content = result.to_json() if output_format == "json" else result.to_pretty_string()
        with open(output, "w", encoding="utf-8") as handle:
            handle.write(content)
        err_console.print(f"[green]Wrote {output_format} readout to[/] {output}")
        return

    # ...or render to the screen.
    if output_format == "json":
        console.print_json(result.to_json())
    else:
        render_rich(result.to_dict())


if __name__ == "__main__":
    main()
