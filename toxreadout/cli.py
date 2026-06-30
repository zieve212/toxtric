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
    from toxreadout import pipeline
except ModuleNotFoundError:
    import pipeline

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

    # Reasoning for each matched target.
    for target in data["targets"]:
        if target["interaction"].get("match_found"):
            risk = target["risk_assessment"]
            console.print(
                f"  [dim]{target['protein'].get('accession')}:[/] {risk.get('reasoning')}"
            )

    lit = data.get("literature")
    if lit:
        console.print(
            f'\n[dim]Literature:[/] {lit.get("pubmed_hits")} PubMed hits '
            f'for "{lit.get("query")}"'
        )


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


@click.command()
@click.option("--smiles", "-s", required=True,
              help="Compound SMILES string.")
@click.option("--fasta", "-f", "fasta", required=True, multiple=True,
              help="Protein FASTA file path or raw sequence/accession. "
                   "Repeat -f for multiple targets.")
@click.option("--output", "-o", type=click.Path(), default=None,
              help="Write the readout to this file instead of the screen.")
@click.option("--format", "output_format",
              type=click.Choice(["pretty", "json"]), default=None,
              help="Output format. Defaults to pretty on screen, json to file.")
@click.option("--no-literature", is_flag=True, default=False,
              help="Skip the PubMed literature search (faster).")
def main(smiles, fasta, output, output_format, no_literature):
    """Generate a toxicology readout for a compound against protein target(s)."""
    fasta_inputs = _load_fasta_inputs(fasta)

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
