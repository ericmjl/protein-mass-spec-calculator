# Mass Spectrometry Calculator for Proteins: Design Document

## Project Overview

This tool will provide both a command-line interface (CLI) and web application for protein mass spectrometry analysis, allowing users to:

1. Calculate the mass of intact proteins
2. Predict peptide fragments from enzymatic digests (trypsin and others)
3. Calculate expected mass-to-charge (m/z) ratios for these fragments
4. Simulate mass spectrometry spectra based on these calculations

## Core Functionality Requirements

1. **Protein Mass Calculation**
   - Calculate molecular weight from amino acid sequence
   - Account for post-translational modifications (PTMs)
   - Handle different charge states

2. **Enzymatic Digestion Simulation**
   - Trypsin (cleaves after K and R, except before P)
   - Chymotrypsin (cleaves after F, Y, W)
   - Pepsin (cleaves after F, L, W, Y, A, E, Q)
   <!-- NOTE: THE FOLLOWING ARE NOT GOING TO BE IMPLEMENTED IN THE FIRST VERSION! NEED TO CHECK WITH SME FIRST! -->
   - Asp-N (cleaves before D)
   - Glu-C (cleaves after E, sometimes D)
   - Lys-C (cleaves after K)
   - Arg-C (cleaves after R)
   - Custom cleavage rules

3. **Mass-to-Charge (m/z) Calculation**
   - Calculate m/z for each peptide fragment
   - Account for different charge states (+1, +2, +3, etc.)
   - Consider isotopic distributions

4. **Spectrum Simulation**
   - Generate theoretical peak patterns
   - Account for relative intensities
   - Visualize expected spectra

## Command Line Interface (CLI) Design

1. **Basic Structure**
   - Main command: `mspcalc`
   - Subcommands for different functionalities

2. **Subcommands**
   - `protein-mass`: Calculate mass of intact protein
   - `digest`: Perform enzymatic digestion
   - `mz-calc`: Calculate m/z ratios
   - `simulate`: Generate simulated spectrum

3. **Input Options**
   - Sequence input (direct or from FASTA file)
   - Enzyme selection
   - Charge state specification
   - Missed cleavages allowance
   - PTM definitions

4. **Output Options**
   - Format (CSV, TSV, JSON)
   - Visualization (PNG, SVG, PDF)
   - Verbosity levels

5. **Example Usage Patterns**
   - `mspcalc protein-mass --sequence MKWVTFISLLLLFSSAYS`
   - `mspcalc digest --fasta protein.fasta --enzyme trypsin --missed-cleavages 2`
   - `mspcalc mz-calc --digest-results digest.json --charges 1,2,3`
   - `mspcalc simulate --mz-results mz.json --output spectrum.png`

## Web Application Design

1. **Overall Layout**
   - Clean, modern interface with responsive design
   - Tabbed interface for different functionalities
   - Results display area with visualization capabilities

2. **Input Section**
   - Text area for sequence input
   - File upload for FASTA
   - Dropdown menus for enzyme selection
   - Numeric inputs for charge states, missed cleavages
   - Advanced options collapsible panel

3. **Results Section**
   - Interactive table of results (sortable, filterable)
   - Mass spectrum visualization with zoom/pan capabilities
   - Export options (CSV, PNG, PDF)

4. **Visualization Features**
   - Interactive spectrum plot
   - Ability to overlay multiple spectra
   - Peak labeling
   - Customizable display options (colors, scales)

5. **User Experience Considerations**
   - Real-time calculations for small proteins
   - Progress indicators for larger calculations
   - Persistent storage of recent calculations
   - Shareable result links

## Technical Architecture

1. **Core Library**
   - Python package with core calculation functions
   - Modular design for extensibility
   - Comprehensive test suite

2. **CLI Implementation**
   - Built with Typer
   - Rich for terminal output formatting
   - Matplotlib/Seaborn for visualization

3. **Web Application**
   - FastAPI backend
   - HTMX for interactive frontend
   - SQLite for lightweight persistent storage
   - Matplotlib/Seaborn for visualization

4. **Data Flow**
   - Input validation → Calculation → Results formatting → Visualization

5. **Deployment Considerations**
   - Containerization with Docker
   - CI/CD pipeline
   - Performance optimization for larger proteins

## Development Phases

1. **Phase 1: Core Functionality**
   - Implement basic protein mass calculation
   - Implement trypsin digestion
   - Basic CLI interface

2. **Phase 2: Extended Functionality**
   - Add additional enzymes
   - Implement m/z calculations
   - Enhance CLI capabilities

3. **Phase 3: Visualization**
   - Implement spectrum simulation
   - Add visualization to CLI output

4. **Phase 4: Web Application**
   - Develop FastAPI backend
   - Create HTMX frontend
   - Implement interactive visualizations

5. **Phase 5: Refinement**
   - Performance optimization
   - User experience improvements
   - Additional features based on feedback

## Future Expansion Possibilities

1. Batch processing capabilities
2. Integration with experimental data for comparison
3. Machine learning for spectrum prediction refinement
4. Database of common proteins and their expected spectra
5. Plugin system for custom calculations
