"""Tests for mspcalc package."""

import pytest

from mspcalc.core import (
    AMINO_ACID_WEIGHTS,
    calculate_mz,
    calculate_protein_mass,
    digest_protein,
    get_supported_enzymes,
)


def test_amino_acid_weights():
    """Test that amino acid weights are defined."""
    assert len(AMINO_ACID_WEIGHTS) == 20
    assert "A" in AMINO_ACID_WEIGHTS
    assert "R" in AMINO_ACID_WEIGHTS
    assert AMINO_ACID_WEIGHTS["A"] == pytest.approx(71.03711)


def test_calculate_protein_mass():
    """Test protein mass calculation."""
    # Alanine (A) has a mass of 71.03711 Da
    mass = calculate_protein_mass("A")
    # Mass should be A + water (18.01056)
    assert mass == pytest.approx(71.03711 + 18.01056)

    # Test a longer sequence
    mass = calculate_protein_mass("AAA")
    # Mass should be 3*A + water
    assert mass == pytest.approx(3 * 71.03711 + 18.01056)

    # Test with invalid amino acid
    with pytest.raises(ValueError):
        calculate_protein_mass("AXA")


def test_calculate_mz():
    """Test m/z calculation."""
    # Mass of 1000 Da with charge +1
    mz = calculate_mz(1000.0, 1)
    # m/z = (mass + charge*proton_mass) / charge
    assert mz == pytest.approx(1001.00728)

    # Mass of 1000 Da with charge +2
    mz = calculate_mz(1000.0, 2)
    # m/z = (1000 + 2*1.00728) / 2
    assert mz == pytest.approx(501.00364)

    # Invalid charge
    with pytest.raises(ValueError):
        calculate_mz(1000.0, 0)


def test_digest_protein():
    """Test protein digestion."""
    # Trypsin cleaves after K and R, except before P
    peptides = digest_protein("AKRGP", enzyme="trypsin")
    assert len(peptides) == 2
    assert peptides[0].sequence == "A"
    assert peptides[1].sequence == "KRGP"

    # Test with missed cleavages
    peptides = digest_protein("AKRGP", enzyme="trypsin", missed_cleavages=1)
    assert len(peptides) == 3
    # Original peptides plus one with a missed cleavage
    sequences = [p.sequence for p in peptides]
    assert "A" in sequences
    assert "KRGP" in sequences
    assert "AKRGP" in sequences

    # Invalid enzyme
    with pytest.raises(ValueError):
        digest_protein("AKRP", enzyme="invalid_enzyme")


def test_get_supported_enzymes():
    """Test getting supported enzymes."""
    enzymes = get_supported_enzymes()
    assert len(enzymes) >= 7  # At least the 7 enzymes from the design doc
    assert "trypsin" in enzymes
    assert "chymotrypsin" in enzymes
    assert "pepsin" in enzymes
