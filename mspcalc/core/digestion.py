"""Enzymatic digestion functionality for protein mass spec analysis."""

from typing import List, Optional

from .amino_acids import ENZYME_CLEAVAGE_RULES
from .protein import validate_sequence


class Peptide:
    """Class representing a peptide fragment.

    :param sequence: The peptide sequence
    :param start_pos: The starting position of the peptide
    :param end_pos: The ending position of the peptide
    """

    def __init__(self, sequence: str, start_pos: int, end_pos: int):
        self.sequence = sequence
        self.start_pos = start_pos
        self.end_pos = end_pos

    def __str__(self) -> str:
        """String representation of the peptide.

        :return: String representation of the peptide
        """
        return f"{self.sequence} ({self.start_pos}-{self.end_pos})"

    def __repr__(self) -> str:
        """Representation of the peptide.

        :return: Representation of the peptide
        """
        return f"Peptide({repr(self.sequence)}, {self.start_pos}, {self.end_pos})"


def digest_protein(
    sequence: str,
    enzyme: str = "trypsin",
    missed_cleavages: int = 0,
    min_length: int = 4,
    max_length: Optional[int] = None,
) -> List[Peptide]:
    """Digest a protein sequence using the specified enzyme.

    :param sequence: The protein sequence to digest
    :param enzyme: The enzyme to use for digestion
    :param missed_cleavages: Number of missed cleavages to allow
    :param min_length: Minimum peptide length to include in results
    :param max_length: Maximum peptide length to include in results
    :return: List of Peptide objects representing the fragments
    :raises ValueError: If enzyme is not supported
    """
    sequence = validate_sequence(sequence)

    if enzyme not in ENZYME_CLEAVAGE_RULES:
        valid_enzymes = ", ".join(ENZYME_CLEAVAGE_RULES.keys())
        raise ValueError(
            f"Unsupported enzyme: {enzyme}. Valid options are: {valid_enzymes}"
        )

    # Get cleavage rules for the enzyme
    cleavage_sites, blocked_by = ENZYME_CLEAVAGE_RULES[enzyme]

    # Find all potential cleavage sites
    sites = []

    # Special handling for Asp-N (cleaves before D)
    if enzyme == "asp-n":
        for i in range(1, len(sequence)):
            if sequence[i] in blocked_by:
                sites.append(i)
    else:
        for i in range(len(sequence) - 1):
            if sequence[i] in cleavage_sites and sequence[i + 1] not in blocked_by:
                sites.append(i + 1)

    # Always include the beginning and end of the sequence
    all_sites = [0] + sites + [len(sequence)]

    # Generate peptides with the requested number of missed cleavages
    peptides = []

    for i in range(len(all_sites) - 1):
        for mc in range(min(missed_cleavages + 1, len(all_sites) - i - 1)):
            start = all_sites[i]
            end = all_sites[i + mc + 1]
            peptide_seq = sequence[start:end]

            if min_length <= len(peptide_seq) and (
                max_length is None or len(peptide_seq) <= max_length
            ):
                peptides.append(Peptide(peptide_seq, start, end - 1))

    return peptides


def get_supported_enzymes() -> List[str]:
    """Get a list of supported enzymes for digestion.

    :return: List of supported enzyme names
    """
    return list(ENZYME_CLEAVAGE_RULES.keys())
