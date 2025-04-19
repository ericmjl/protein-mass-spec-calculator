"""Protein mass calculation functionality."""

from pathlib import Path
from typing import Dict, Optional, Union

from .amino_acids import (
    AMINO_ACID_WEIGHTS,
    NON_STANDARD_AMINO_ACID_WEIGHTS,
    PROTON_MASS,
    WATER_MASS,
)


def validate_sequence(sequence: str) -> str:
    """Validate an amino acid sequence.

    :param sequence: The amino acid sequence to validate
    :return: The validated sequence
    :raises ValueError: If sequence contains invalid amino acid characters
    """
    sequence = sequence.strip().upper()
    valid_amino_acids = set(AMINO_ACID_WEIGHTS.keys()) | set(
        NON_STANDARD_AMINO_ACID_WEIGHTS.keys()
    )

    invalid_chars = [
        char for char in sequence if char not in valid_amino_acids and char != "*"
    ]
    if invalid_chars:
        raise ValueError(
            f"Invalid amino acid(s) in sequence: {', '.join(set(invalid_chars))}"
        )

    return sequence


def read_fasta(fasta_path: Union[str, Path]) -> Dict[str, str]:
    """Read a fasta file and return a dictionary of sequence identifiers to sequences.

    :param fasta_path: Path to the fasta file
    :return: Dictionary mapping sequence identifiers to their sequences
    :raises FileNotFoundError: If the fasta file does not exist
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    sequences = {}
    current_id = None
    current_sequence = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous sequence if it exists
                if current_id is not None:
                    sequences[current_id] = "".join(current_sequence)

                # Start a new sequence
                current_id = line[1:].split()[0]  # Take the first word after '>' as ID
                current_sequence = []
            else:
                if current_id is not None:
                    current_sequence.append(line)

    # Add the last sequence
    if current_id is not None:
        sequences[current_id] = "".join(current_sequence)

    return sequences


def calculate_protein_mass(
    sequence: str, ptms: Optional[Dict[int, float]] = None, monoisotopic: bool = True
) -> float:
    """Calculate the mass of a protein from its amino acid sequence.

    :param sequence: The amino acid sequence
    :param ptms: Dictionary mapping positions (0-indexed) to mass modifications
    :param monoisotopic: Whether to use monoisotopic masses (True)
        or average masses (False)
    :return: The calculated protein mass in Daltons
    """
    sequence = validate_sequence(sequence)

    # Calculate the base mass from amino acids
    mass = sum(
        AMINO_ACID_WEIGHTS.get(aa, NON_STANDARD_AMINO_ACID_WEIGHTS.get(aa, 0))
        for aa in sequence
    )

    # Add the mass of a water molecule (peptide bond formation removes water)
    mass += WATER_MASS

    # Add PTMs if provided
    if ptms:
        for _, ptm_mass in ptms.items():
            mass += ptm_mass

    return mass


def calculate_mz(mass: float, charge: int) -> float:
    """Calculate the mass-to-charge ratio (m/z) for a given mass and charge.

    :param mass: The mass in Daltons
    :param charge: The charge state (integer)
    :return: The calculated m/z value
    :raises ValueError: If charge is less than or equal to 0
    """
    if charge <= 0:
        raise ValueError("Charge must be greater than 0")

    # m/z = (mass + (charge * proton mass)) / charge
    return (mass + (charge * PROTON_MASS)) / charge
