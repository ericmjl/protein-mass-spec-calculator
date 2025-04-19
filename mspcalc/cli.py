"""Command line interface for the Mass Spectrometry Calculator for Proteins."""

import json
import sys
from enum import Enum
from pathlib import Path
from typing import List, Optional

import typer
from loguru import logger
from rich.console import Console
from rich.table import Table

from mspcalc.core import (
    COMMON_PTMS,
    MassSpectrum,
    Peptide,
    calculate_mz,
    calculate_protein_mass,
    digest_protein,
    generate_theoretical_spectrum,
    get_supported_enzymes,
    read_fasta,
    validate_sequence,
)

# Configure logger
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",  # noqa: E501
)

# Create Typer app
app = typer.Typer(
    name="mspcalc",
    help="Mass Spectrometry Calculator for Proteins",
    add_completion=False,
)

# Console for rich output
console = Console()


class OutputFormat(str, Enum):
    """Supported output formats."""

    CSV = "csv"
    TSV = "tsv"
    JSON = "json"


class PlotFormat(str, Enum):
    """Supported plot formats."""

    PNG = "png"
    SVG = "svg"
    PDF = "pdf"


@app.command(name="protein-mass")
def protein_mass(
    sequence: Optional[str] = typer.Option(
        None, "--sequence", "-s", help="Amino acid sequence of the protein"
    ),
    fasta: Optional[Path] = typer.Option(
        None, "--fasta", "-f", help="Path to a FASTA file containing protein sequences"
    ),
    ptms: List[str] = typer.Option(
        [],
        "--ptm",
        help="PTMs to apply, format: 'position:modification' "
        "(e.g., '42:phosphorylation')",
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.JSON, "--format", "-fmt", help="Output format"
    ),
    output_file: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Path to output file"
    ),
):
    """Calculate the mass of intact proteins."""
    if sequence is None and fasta is None:
        logger.error("Either --sequence or --fasta must be provided")
        raise typer.Exit(code=1)

    # Process PTMs
    ptm_dict = {}
    for ptm in ptms:
        try:
            pos, mod = ptm.split(":")
            pos = int(pos)
            if mod not in COMMON_PTMS:
                valid_ptms = ", ".join(COMMON_PTMS.keys())
                logger.error(f"Unknown PTM: {mod}. Valid options are: {valid_ptms}")
                raise typer.Exit(code=1)
            ptm_dict[pos] = COMMON_PTMS[mod]
        except ValueError:
            logger.error(f"Invalid PTM format: {ptm}. Use 'position:modification'")
            raise typer.Exit(code=1)

    results = []

    if sequence:
        try:
            sequence = validate_sequence(sequence)
            mass = calculate_protein_mass(sequence, ptm_dict)
            results.append({"sequence": sequence, "mass": mass})
        except ValueError as e:
            logger.error(f"Error processing sequence: {e}")
            raise typer.Exit(code=1)

    if fasta:
        try:
            sequences = read_fasta(fasta)
            for seq_id, seq in sequences.items():
                try:
                    seq = validate_sequence(seq)
                    mass = calculate_protein_mass(seq, ptm_dict)
                    results.append({"id": seq_id, "sequence": seq, "mass": mass})
                except ValueError as e:
                    logger.warning(f"Error processing sequence {seq_id}: {e}")
        except FileNotFoundError:
            logger.error(f"FASTA file not found: {fasta}")
            raise typer.Exit(code=1)

    # Output results
    if output_format == OutputFormat.JSON:
        output = json.dumps(results, indent=2)
    elif output_format == OutputFormat.CSV:
        headers = ["id", "sequence", "mass"]
        rows = []
        for r in results:
            row = [r.get("id", ""), r.get("sequence", ""), r.get("mass", "")]
            rows.append(",".join([str(x) for x in row]))
        output = "\n".join([",".join(headers)] + rows)
    elif output_format == OutputFormat.TSV:
        headers = ["id", "sequence", "mass"]
        rows = []
        for r in results:
            row = [r.get("id", ""), r.get("sequence", ""), r.get("mass", "")]
            rows.append("\t".join([str(x) for x in row]))
        output = "\n".join(["\t".join(headers)] + rows)

    if output_file:
        output_file.write_text(output)
        logger.info(f"Results written to {output_file}")
    else:
        console.print(output)


@app.command(name="digest")
def digest(
    sequence: Optional[str] = typer.Option(
        None, "--sequence", "-s", help="Amino acid sequence of the protein"
    ),
    fasta: Optional[Path] = typer.Option(
        None, "--fasta", "-f", help="Path to a FASTA file containing protein sequences"
    ),
    enzyme: str = typer.Option(
        "trypsin", "--enzyme", "-e", help="Enzyme to use for digestion"
    ),
    missed_cleavages: int = typer.Option(
        0, "--missed-cleavages", "-mc", help="Number of missed cleavages to allow"
    ),
    min_length: int = typer.Option(
        4, "--min-length", "-min", help="Minimum peptide length"
    ),
    max_length: Optional[int] = typer.Option(
        None, "--max-length", "-max", help="Maximum peptide length"
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.JSON, "--format", "-fmt", help="Output format"
    ),
    output_file: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Path to output file"
    ),
):
    """Perform enzymatic digestion of proteins."""
    if sequence is None and fasta is None:
        logger.error("Either --sequence or --fasta must be provided")
        raise typer.Exit(code=1)

    # Validate enzyme
    supported_enzymes = get_supported_enzymes()
    if enzyme not in supported_enzymes:
        logger.error(
            f"Unsupported enzyme: {enzyme}. "
            f"Valid options: {', '.join(supported_enzymes)}"
        )
        raise typer.Exit(code=1)

    results = []

    if sequence:
        try:
            sequence = validate_sequence(sequence)
            peptides = digest_protein(
                sequence, enzyme, missed_cleavages, min_length, max_length
            )
            peptide_data = [
                {
                    "sequence": p.sequence,
                    "start": p.start_pos,
                    "end": p.end_pos,
                    "length": len(p.sequence),
                    "mass": calculate_protein_mass(p.sequence),
                }
                for p in peptides
            ]
            results.append({"protein": sequence, "peptides": peptide_data})
        except ValueError as e:
            logger.error(f"Error processing sequence: {e}")
            raise typer.Exit(code=1)

    if fasta:
        try:
            sequences = read_fasta(fasta)
            for seq_id, seq in sequences.items():
                try:
                    seq = validate_sequence(seq)
                    peptides = digest_protein(
                        seq, enzyme, missed_cleavages, min_length, max_length
                    )
                    peptide_data = [
                        {
                            "sequence": p.sequence,
                            "start": p.start_pos,
                            "end": p.end_pos,
                            "length": len(p.sequence),
                            "mass": calculate_protein_mass(p.sequence),
                        }
                        for p in peptides
                    ]
                    results.append(
                        {"id": seq_id, "protein": seq, "peptides": peptide_data}
                    )
                except ValueError as e:
                    logger.warning(f"Error processing sequence {seq_id}: {e}")
        except FileNotFoundError:
            logger.error(f"FASTA file not found: {fasta}")
            raise typer.Exit(code=1)

    # Output results
    if output_format == OutputFormat.JSON:
        output = json.dumps(results, indent=2)
    elif output_format == OutputFormat.CSV:
        headers = ["protein_id", "peptide_sequence", "start", "end", "length", "mass"]
        rows = []
        for r in results:
            for p in r.get("peptides", []):
                row = [
                    r.get("id", ""),
                    p.get("sequence", ""),
                    p.get("start", ""),
                    p.get("end", ""),
                    p.get("length", ""),
                    p.get("mass", ""),
                ]
                rows.append(",".join([str(x) for x in row]))
        output = "\n".join([",".join(headers)] + rows)
    elif output_format == OutputFormat.TSV:
        headers = ["protein_id", "peptide_sequence", "start", "end", "length", "mass"]
        rows = []
        for r in results:
            for p in r.get("peptides", []):
                row = [
                    r.get("id", ""),
                    p.get("sequence", ""),
                    p.get("start", ""),
                    p.get("end", ""),
                    p.get("length", ""),
                    p.get("mass", ""),
                ]
                rows.append("\t".join([str(x) for x in row]))
        output = "\n".join(["\t".join(headers)] + rows)

    if output_file:
        output_file.write_text(output)
        logger.info(f"Results written to {output_file}")
    else:
        console.print(output)


@app.command(name="mz-calc")
def mz_calc(
    sequence: Optional[str] = typer.Option(
        None, "--sequence", "-s", help="Amino acid sequence"
    ),
    digest_results: Optional[Path] = typer.Option(
        None, "--digest-results", "-d", help="Path to digest results JSON file"
    ),
    charges: List[int] = typer.Option(
        [1, 2, 3], "--charge", "-c", help="Charge states to calculate"
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.JSON, "--format", "-fmt", help="Output format"
    ),
    output_file: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Path to output file"
    ),
):
    """Calculate m/z ratios for peptides."""
    if sequence is None and digest_results is None:
        logger.error("Either --sequence or --digest-results must be provided")
        raise typer.Exit(code=1)

    # Validate charges
    if not charges or any(c <= 0 for c in charges):
        logger.error("Charge states must be positive integers")
        raise typer.Exit(code=1)

    results = []

    if sequence:
        try:
            sequence = validate_sequence(sequence)
            mass = calculate_protein_mass(sequence)
            mz_values = []
            for charge in charges:
                mz = calculate_mz(mass, charge)
                mz_values.append({"charge": charge, "mz": mz})
            results.append({"sequence": sequence, "mass": mass, "mz_values": mz_values})
        except ValueError as e:
            logger.error(f"Error processing sequence: {e}")
            raise typer.Exit(code=1)

    if digest_results:
        try:
            digest_data = json.loads(digest_results.read_text())
            for protein_data in digest_data:
                protein_id = protein_data.get("id", "")
                protein_mz = []

                for peptide_data in protein_data.get("peptides", []):
                    peptide_seq = peptide_data.get("sequence", "")
                    mass = peptide_data.get("mass", calculate_protein_mass(peptide_seq))

                    mz_values = []
                    for charge in charges:
                        mz = calculate_mz(mass, charge)
                        mz_values.append({"charge": charge, "mz": mz})

                    protein_mz.append(
                        {
                            "sequence": peptide_seq,
                            "start": peptide_data.get("start", 0),
                            "end": peptide_data.get("end", 0),
                            "mass": mass,
                            "mz_values": mz_values,
                        }
                    )

                results.append({"id": protein_id, "peptides": protein_mz})
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading digest results: {e}")
            raise typer.Exit(code=1)

    # Output results
    if output_format == OutputFormat.JSON:
        output = json.dumps(results, indent=2)
    elif output_format == OutputFormat.CSV:
        headers = ["protein_id", "peptide_sequence", "mass", "charge", "mz"]
        rows = []
        for r in results:
            if "peptides" in r:  # From digest results
                for p in r.get("peptides", []):
                    for mz_value in p.get("mz_values", []):
                        row = [
                            r.get("id", ""),
                            p.get("sequence", ""),
                            p.get("mass", ""),
                            mz_value.get("charge", ""),
                            mz_value.get("mz", ""),
                        ]
                        rows.append(",".join([str(x) for x in row]))
            else:  # From single sequence
                for mz_value in r.get("mz_values", []):
                    row = [
                        "",
                        r.get("sequence", ""),
                        r.get("mass", ""),
                        mz_value.get("charge", ""),
                        mz_value.get("mz", ""),
                    ]
                    rows.append(",".join([str(x) for x in row]))
        output = "\n".join([",".join(headers)] + rows)
    elif output_format == OutputFormat.TSV:
        headers = ["protein_id", "peptide_sequence", "mass", "charge", "mz"]
        rows = []
        for r in results:
            if "peptides" in r:  # From digest results
                for p in r.get("peptides", []):
                    for mz_value in p.get("mz_values", []):
                        row = [
                            r.get("id", ""),
                            p.get("sequence", ""),
                            p.get("mass", ""),
                            mz_value.get("charge", ""),
                            mz_value.get("mz", ""),
                        ]
                        rows.append("\t".join([str(x) for x in row]))
            else:  # From single sequence
                for mz_value in r.get("mz_values", []):
                    row = [
                        "",
                        r.get("sequence", ""),
                        r.get("mass", ""),
                        mz_value.get("charge", ""),
                        mz_value.get("mz", ""),
                    ]
                    rows.append("\t".join([str(x) for x in row]))
        output = "\n".join(["\t".join(headers)] + rows)

    if output_file:
        output_file.write_text(output)
        logger.info(f"Results written to {output_file}")
    else:
        console.print(output)


@app.command(name="simulate")
def simulate(
    sequence: Optional[str] = typer.Option(
        None, "--sequence", "-s", help="Amino acid sequence"
    ),
    mz_results: Optional[Path] = typer.Option(
        None, "--mz-results", "-m", help="Path to m/z results JSON file"
    ),
    digest_results: Optional[Path] = typer.Option(
        None, "--digest-results", "-d", help="Path to digest results JSON file"
    ),
    enzyme: Optional[str] = typer.Option(
        None, "--enzyme", "-e", help="Enzyme to use for digestion"
    ),
    missed_cleavages: int = typer.Option(
        0, "--missed-cleavages", "-mc", help="Number of missed cleavages to allow"
    ),
    charges: List[int] = typer.Option(
        [1, 2, 3], "--charge", "-c", help="Charge states to simulate"
    ),
    output_file: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Path to output plot file"
    ),
    plot_format: PlotFormat = typer.Option(
        PlotFormat.PNG, "--plot-format", "-pf", help="Format for the output plot"
    ),
    title: Optional[str] = typer.Option(
        None, "--title", "-t", help="Title for the spectrum plot"
    ),
    no_display: bool = typer.Option(
        False, "--no-display", help="Don't display the plot (just save to file)"
    ),
):
    """Generate simulated mass spectra."""
    # Input validation
    input_sources = sum(
        1 for x in [sequence, mz_results, digest_results] if x is not None
    )
    if input_sources == 0:
        logger.error(
            "One of --sequence, --mz-results or --digest-results must be provided"
        )
        raise typer.Exit(code=1)

    if input_sources > 1:
        logger.warning(
            "Multiple input sources provided. Using the first available in order: "
            "sequence, digest-results, mz-results"
        )

    # Prepare output file path if provided
    if output_file:
        output_file = Path(output_file)
        if not output_file.suffix:
            output_file = output_file.with_suffix(f".{plot_format.value}")

    # Process input and generate spectrum
    spectrum = None

    # Case 1: Direct sequence input with optional digestion
    if sequence:
        try:
            sequence = validate_sequence(sequence)

            # If enzyme is provided, digest the sequence first
            if enzyme:
                supported_enzymes = get_supported_enzymes()
                if enzyme not in supported_enzymes:
                    logger.error(
                        f"Unsupported enzyme: {enzyme}. "
                        f"Valid options: {', '.join(supported_enzymes)}"
                    )
                    raise typer.Exit(code=1)

                peptides = digest_protein(sequence, enzyme, missed_cleavages)
                spectrum_title = title or f"Theoretical Spectrum: {enzyme} digest"
                spectrum = generate_theoretical_spectrum(
                    peptides, charges, title=spectrum_title
                )
            else:
                # Just calculate m/z for the intact protein
                mass = calculate_protein_mass(sequence)
                spectrum = MassSpectrum(title=title or "Theoretical Spectrum")

                for charge in charges:
                    mz = calculate_mz(mass, charge)
                    intensity = (
                        100.0 / charge
                    )  # Simple model: higher charges have lower intensity
                    spectrum.add_peak(mz, intensity, f"{sequence} ({charge}+)")
        except ValueError as e:
            logger.error(f"Error processing sequence: {e}")
            raise typer.Exit(code=1)

    # Case 2: Digest results input
    elif digest_results:
        try:
            digest_data = json.loads(digest_results.read_text())
            all_peptides = []

            for protein_data in digest_data:
                for peptide_data in protein_data.get("peptides", []):
                    peptide_seq = peptide_data.get("sequence", "")
                    start = peptide_data.get("start", 0)
                    end = peptide_data.get("end", 0)

                    all_peptides.append(Peptide(peptide_seq, start, end))

            if all_peptides:
                spectrum_title = title or "Theoretical Spectrum from Digest Results"
                spectrum = generate_theoretical_spectrum(
                    all_peptides, charges, title=spectrum_title
                )
            else:
                logger.error("No peptides found in digest results")
                raise typer.Exit(code=1)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading digest results: {e}")
            raise typer.Exit(code=1)

    # Case 3: m/z results input
    elif mz_results:
        try:
            mz_data = json.loads(mz_results.read_text())
            spectrum = MassSpectrum(
                title=title or "Theoretical Spectrum from m/z Results"
            )

            for item in mz_data:
                if "peptides" in item:  # From digest results
                    for peptide in item.get("peptides", []):
                        peptide_seq = peptide.get("sequence", "")

                        for mz_value in peptide.get("mz_values", []):
                            mz = mz_value.get("mz", 0)
                            charge = mz_value.get("charge", 1)
                            intensity = 100.0 / (len(peptide_seq) * charge) * 10
                            spectrum.add_peak(
                                mz, intensity, f"{peptide_seq} ({charge}+)"
                            )
                else:  # From single sequence
                    seq = item.get("sequence", "")

                    for mz_value in item.get("mz_values", []):
                        mz = mz_value.get("mz", 0)
                        charge = mz_value.get("charge", 1)
                        intensity = 100.0 / charge
                        spectrum.add_peak(mz, intensity, f"{seq} ({charge}+)")
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading m/z results: {e}")
            raise typer.Exit(code=1)

    # Plot the spectrum
    if spectrum:
        if not spectrum.peaks:
            logger.error("No peaks to display")
            raise typer.Exit(code=1)

        try:
            spectrum.plot(output_file, show=not no_display)
            if output_file:
                logger.info(f"Plot saved to {output_file}")
        except Exception as e:
            logger.error(f"Error generating plot: {e}")
            raise typer.Exit(code=1)
    else:
        logger.error("Failed to generate spectrum")
        raise typer.Exit(code=1)


@app.command(name="enzymes")
def list_enzymes():
    """List all supported enzymes for digestion."""
    enzymes = get_supported_enzymes()

    table = Table(title="Supported Enzymes")
    table.add_column("Enzyme", style="cyan")
    table.add_column("Description", style="green")

    descriptions = {
        "trypsin": "Cleaves after K and R, except before P",
        "chymotrypsin": "Cleaves after F, Y, W",
        "pepsin": "Cleaves after F, L, W, Y, A, E, Q",
        "asp-n": "Cleaves before D",
        "glu-c": "Cleaves after E, sometimes D",
        "lys-c": "Cleaves after K",
        "arg-c": "Cleaves after R",
    }

    for enzyme in sorted(enzymes):
        table.add_row(enzyme, descriptions.get(enzyme, ""))

    console.print(table)


if __name__ == "__main__":
    app()
