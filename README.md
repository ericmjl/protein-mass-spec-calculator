# Protein Mass Spectrometry Calculator

A powerful command-line tool and library for protein mass spectrometry analysis.

Made with ❤️ by Eric J. Ma (@ericmjl).

## Features

- Calculate molecular weights of proteins from amino acid sequences
- Simulate enzymatic digestion (trypsin, chymotrypsin, pepsin, and more)
- Calculate m/z ratios for protein fragments
- Generate theoretical mass spectra
- Support for post-translational modifications (PTMs)
- FASTA file support for batch processing
- Multiple output formats (JSON, CSV, TSV)
- Visualization of mass spectra

## Installation

To get started with development:

```bash
git clone git@github.com:ericmjl/protein-mass-spec-calculator
cd protein-mass-spec-calculator
pixi install
```

## CLI Usage

The `mspcalc` command-line tool provides several subcommands for protein mass spectrometry analysis:

### Calculate Protein Mass

Calculate the molecular weight of a protein from its amino acid sequence:

```bash
mspcalc protein-mass --sequence MKWVTFISLLLLFSSAYS
```

Or calculate masses for proteins in a FASTA file:

```bash
mspcalc protein-mass --fasta proteins.fasta --output results.json
```

With post-translational modifications:

```bash
mspcalc protein-mass --sequence MKWVTFISLLLLFSSAYS --ptm 5:phosphorylation --ptm 10:methylation
```

### Enzymatic Digestion

Digest a protein with a specific enzyme:

```bash
mspcalc digest --sequence MKWVTFISLLLLFSSAYS --enzyme trypsin --missed-cleavages 1
```

List all supported enzymes:

```bash
mspcalc enzymes
```

### Calculate m/z Ratios

Calculate m/z ratios for a protein or peptide with different charge states:

```bash
mspcalc mz-calc --sequence MKWVTFISLLLLFSSAYS --charge 1 --charge 2 --charge 3
```

Or process the results of a digestion:

```bash
mspcalc digest --sequence MKWVTFISLLLLFSSAYS --enzyme trypsin --output digest.json
mspcalc mz-calc --digest-results digest.json --charge 1 --charge 2
```

### Simulate Mass Spectrum

Generate a theoretical mass spectrum:

```bash
mspcalc simulate --sequence MKWVTFISLLLLFSSAYS --enzyme trypsin --output spectrum.png
```

Process m/z calculation results to create a spectrum:

```bash
mspcalc mz-calc --sequence MKWVTFISLLLLFSSAYS --output mz.json
mspcalc simulate --mz-results mz.json --output spectrum.png
```

## Output Formats

The tool supports multiple output formats:

```bash
mspcalc protein-mass --sequence MKWVTFISLLLLFSSAYS --format json
mspcalc protein-mass --sequence MKWVTFISLLLLFSSAYS --format csv
mspcalc protein-mass --sequence MKWVTFISLLLLFSSAYS --format tsv
```

For full documentation on each command, use the `--help` option:

```bash
mspcalc --help
mspcalc protein-mass --help
mspcalc digest --help
mspcalc mz-calc --help
mspcalc simulate --help
```

## Python API

You can also use the library in your Python code:

```python
from mspcalc import calculate_protein_mass, digest_protein, generate_theoretical_spectrum

# Calculate protein mass
mass = calculate_protein_mass("MKWVTFISLLLLFSSAYS")
print(f"Protein mass: {mass} Da")

# Digest a protein
peptides = digest_protein("MKWVTFISLLLLFSSAYS", enzyme="trypsin", missed_cleavages=1)
for peptide in peptides:
    print(f"Peptide: {peptide}")

# Generate a spectrum
spectrum = generate_theoretical_spectrum(peptides)
spectrum.plot("spectrum.png")
```

## License

[MIT License]
