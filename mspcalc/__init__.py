"""Mass Spectrometry Calculator for Proteins.

A command-line interface (CLI) and library for protein mass spectrometry analysis.
"""

__version__ = "0.1.0"

from mspcalc.core import (
    # Constants
    AMINO_ACID_WEIGHTS,
    COMMON_PTMS,
    ENZYME_CLEAVAGE_RULES,
    NON_STANDARD_AMINO_ACID_WEIGHTS,
    PROTON_MASS,
    WATER_MASS,
    MassSpectrum,
    # Spectrum analysis
    Peak,
    # Digestion
    Peptide,
    calculate_mz,
    calculate_peptide_mz_values,
    # Protein calculation
    calculate_protein_mass,
    digest_protein,
    generate_theoretical_spectrum,
    get_supported_enzymes,
    read_fasta,
    validate_sequence,
)

__all__ = [
    "calculate_protein_mass",
    "calculate_mz",
    "read_fasta",
    "validate_sequence",
    "Peptide",
    "digest_protein",
    "get_supported_enzymes",
    "Peak",
    "MassSpectrum",
    "calculate_peptide_mz_values",
    "generate_theoretical_spectrum",
    "AMINO_ACID_WEIGHTS",
    "NON_STANDARD_AMINO_ACID_WEIGHTS",
    "ENZYME_CLEAVAGE_RULES",
    "COMMON_PTMS",
    "WATER_MASS",
    "PROTON_MASS",
]
