"""Amino acid properties for mass spectrometry calculations."""

from typing import Dict

# Molecular weights of amino acids in Daltons (Da)
# These are monoisotopic masses (most abundant isotope)
AMINO_ACID_WEIGHTS: Dict[str, float] = {
    "A": 71.03711,  # Alanine
    "R": 156.10111,  # Arginine
    "N": 114.04293,  # Asparagine
    "D": 115.02694,  # Aspartic acid
    "C": 103.00919,  # Cysteine
    "E": 129.04259,  # Glutamic acid
    "Q": 128.05858,  # Glutamine
    "G": 57.02146,  # Glycine
    "H": 137.05891,  # Histidine
    "I": 113.08406,  # Isoleucine
    "L": 113.08406,  # Leucine
    "K": 128.09496,  # Lysine
    "M": 131.04049,  # Methionine
    "F": 147.06841,  # Phenylalanine
    "P": 97.05276,  # Proline
    "S": 87.03203,  # Serine
    "T": 101.04768,  # Threonine
    "W": 186.07931,  # Tryptophan
    "Y": 163.06333,  # Tyrosine
    "V": 99.06841,  # Valine
}

# Add non-standard amino acids
NON_STANDARD_AMINO_ACID_WEIGHTS: Dict[str, float] = {
    "U": 150.95364,  # Selenocysteine
    "O": 237.14773,  # Pyrrolysine
}

# Terminal modifications (water molecule)
WATER_MASS = 18.01056  # H2O

# Masses for different elements and common modifications
PROTON_MASS = 1.00728  # Mass of a proton (for charge calculation)

# Enzyme cleavage sites
# Format: Dictionary mapping enzyme name to a tuple of (cleavage_sites, blocked_by)
ENZYME_CLEAVAGE_RULES: Dict[str, tuple] = {
    "trypsin": ({"R", "K"}, {"P"}),  # Cleaves after R and K, unless followed by P
    "chymotrypsin": ({"F", "Y", "W"}, set()),  # Cleaves after F, Y, W
    "pepsin": (
        {"F", "L", "W", "Y", "A", "E", "Q"},
        set(),
    ),  # Cleaves after F, L, W, Y, A, E, Q
    # Additional enzymes as per design doc:
    "asp-n": (set(), {"D"}),  # Cleaves before D (special case)
    "glu-c": ({"E"}, set()),  # Cleaves after E (sometimes also D)
    "lys-c": ({"K"}, set()),  # Cleaves after K
    "arg-c": ({"R"}, set()),  # Cleaves after R
}

# Common PTMs (Post-Translational Modifications)
COMMON_PTMS: Dict[str, float] = {
    "phosphorylation": 79.96633,  # Addition of phosphate group
    "acetylation": 42.01056,  # Addition of acetyl group
    "methylation": 14.01565,  # Addition of methyl group
    "oxidation": 15.99491,  # Oxidation (e.g., of Methionine)
    "carbamidomethylation": 57.02146,  # Common Cysteine modification
}
