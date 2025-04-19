"""Mass spectrometry spectrum simulation functionality."""

from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np

from .digestion import Peptide
from .protein import calculate_mz, calculate_protein_mass


class Peak:
    """Class representing a mass spec peak.

    :param mz: The mass-to-charge ratio (m/z)
    :param intensity: The intensity value
    :param label: Optional label for the peak
    """

    def __init__(self, mz: float, intensity: float, label: Optional[str] = None):
        self.mz = mz
        self.intensity = intensity
        self.label = label

    def __str__(self) -> str:
        """String representation of the peak.

        :return: String representation of the peak
        """
        if self.label:
            return f"{self.label}: m/z={self.mz:.4f}, intensity={self.intensity:.2f}"
        return f"m/z={self.mz:.4f}, intensity={self.intensity:.2f}"

    def __repr__(self) -> str:
        """Representation of the peak.

        :return: Representation of the peak
        """
        if self.label:
            return f"Peak({self.mz:.4f}, {self.intensity:.2f}, {repr(self.label)})"
        return f"Peak({self.mz:.4f}, {self.intensity:.2f})"


class MassSpectrum:
    """Class representing a mass spectrum.

    :param peaks: List of Peak objects
    :param title: Title for the spectrum
    """

    def __init__(
        self, peaks: Optional[List[Peak]] = None, title: str = "Mass Spectrum"
    ):
        self.peaks = peaks or []
        self.title = title

    def add_peak(
        self, mz: float, intensity: float, label: Optional[str] = None
    ) -> None:
        """Add a peak to the spectrum.

        :param mz: The mass-to-charge ratio (m/z)
        :param intensity: The intensity value
        :param label: Optional label for the peak
        """
        self.peaks.append(Peak(mz, intensity, label))

    def get_mz_array(self) -> np.ndarray:
        """Get an array of m/z values.

        :return: NumPy array of m/z values
        """
        return np.array([peak.mz for peak in self.peaks])

    def get_intensity_array(self) -> np.ndarray:
        """Get an array of intensity values.

        :return: NumPy array of intensity values
        """
        return np.array([peak.intensity for peak in self.peaks])

    def plot(
        self,
        output_file: Optional[Union[str, Path]] = None,
        show: bool = True,
        figsize: Tuple[int, int] = (10, 6),
    ) -> None:
        """Plot the mass spectrum.

        :param output_file: Path to save the plot (optional)
        :param show: Whether to display the plot
        :param figsize: Figure size (width, height) in inches
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            sns.set_style("whitegrid")

            mz_array = self.get_mz_array()
            intensity_array = self.get_intensity_array()

            plt.figure(figsize=figsize)
            plt.stem(mz_array, intensity_array, markerfmt=" ", basefmt=" ")
            plt.plot(mz_array, intensity_array, "b-", linewidth=0.8, alpha=0.7)

            # Add labels for the top 10 peaks by intensity
            if self.peaks:
                top_peaks = sorted(self.peaks, key=lambda x: x.intensity, reverse=True)[
                    :10
                ]
                for peak in top_peaks:
                    plt.annotate(
                        f"{peak.mz:.2f}" if not peak.label else peak.label,
                        xy=(peak.mz, peak.intensity),
                        xytext=(0, 5),
                        textcoords="offset points",
                        ha="center",
                        fontsize=8,
                    )

            plt.xlabel("m/z")
            plt.ylabel("Intensity")
            plt.title(self.title)
            plt.grid(True, alpha=0.3)

            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches="tight")

            if show:
                plt.show()
            else:
                plt.close()

        except ImportError:
            print("Plotting requires matplotlib and seaborn to be installed.")


def calculate_peptide_mz_values(
    peptides: List[Peptide], charges: List[int] = [1, 2, 3]
) -> List[Tuple[float, Peptide, int]]:
    """Calculate m/z values for a list of peptides with different charge states.

    :param peptides: List of Peptide objects
    :param charges: List of charge states to calculate
    :return: List of tuples containing (m/z value, peptide, charge)
    """
    results = []

    for peptide in peptides:
        mass = calculate_protein_mass(peptide.sequence)

        for charge in charges:
            mz = calculate_mz(mass, charge)
            results.append((mz, peptide, charge))

    # Sort by m/z value
    return sorted(results, key=lambda x: x[0])


def generate_theoretical_spectrum(
    peptides: List[Peptide],
    charges: List[int] = [1, 2, 3],
    base_intensity: float = 100.0,
    title: str = "Theoretical Mass Spectrum",
) -> MassSpectrum:
    """Generate a theoretical mass spectrum from peptide fragments.

    :param peptides: List of Peptide objects
    :param charges: List of charge states to include
    :param base_intensity: Base intensity value for peaks
    :param title: Title for the spectrum
    :return: MassSpectrum object
    """
    spectrum = MassSpectrum(title=title)

    # Calculate m/z values for all peptides and charge states
    mz_values = calculate_peptide_mz_values(peptides, charges)

    # Create a simple model for intensity based on peptide length and charge
    for mz, peptide, charge in mz_values:
        # Simple intensity model:
        # longer peptides and higher charges have lower intensity
        intensity = base_intensity * (1.0 / len(peptide.sequence)) * (1.0 / charge) * 10

        # Label format: sequence (charge)
        label = f"{peptide.sequence} ({charge}+)"

        spectrum.add_peak(mz, intensity, label)

    return spectrum
