
# Walsh-Floquet Theory of Periodic Kick Drives

This Zenodo record contains the complete codebase, data, and reproducibility materials for the manuscript "Walsh-Floquet Theory of Periodic Kick Drives" by James Walkling and Marin Bukov. The repository provides all source code, analysis notebooks, generated data, and publication figures necessary to reproduce the results presented in the paper.

## Citation

### Plain Text Citation
J. Walkling and M. Bukov, “Walsh-Floquet Theory of Periodic Kick Drives,” arXiv:2505.11071 (2025).

### BibTeX Citation
```bibtex
@article{Walkling2025,
  author       = {James Walkling and Marin Bukov},
  title        = {Walsh-Floquet Theory of Periodic Kick Drives},
  journal      = {arXiv preprint},
  eprint       = {2505.11071},
  archivePrefix= {arXiv},
  primaryClass = {quant-ph},
  year         = {2025},
  month        = may,
}
```

## Links
- **arXiv**: https://doi.org/10.48550/arXiv.2505.11071
- **Zenodo DOI**: TK 

## Table of Contents
- [Repository Overview](#repository-overview)
- [Installation and Setup](#installation-and-setup)
- [Directory Structure](#directory-structure)
- [Data Formats and Organization](#data-formats-and-organization)
- [Reproducing Results](#reproducing-results)
- [Figure Generation](#figure-generation)
- [License Information](#license-information)


## Repository Overview

This repository provides the associated code framework to implement the Floquet extended Hilbert space formalism in the Walsh basis.

## Installation and Setup

### Prerequisites
- Python 3.8 or higher
- Git for repository management
- At least 8GB RAM (16GB+ recommended for larger systems)

### Environment Setup
Create a dedicated Python environment and install all required packages.


## Directory Structure

```
WalshFloquet/
├── README.md                   # This file
├── requirements.txt            # Python package dependencies
├── environment.yml             # yaml file to create python environment
├── LICENSE                     # BSD 3-Clause license for code
├── LICENSE-DATA                # CC-BY 4.0 license for data
├── src_code/                   # Source code implementation
│   ├── README.md               # Source code documentation
│   ├── src/                    # walsh_floquet library
│   └── scripts/                # Standalone execution scripts
├── data/                       # Generated and processed data
│   ├── README.md              # Data organization documentation
│   └── processed/              # Processed data for plotting
├── visual_elements/            # Publication figures and assets
│   ├── README.md               # Figure documentation
│   ├── figs/                   # Manuscript figures
└── src_latex/                  # LaTeX source for manuscript
```

## Data Formats and Organization

### Data Formats
Data are given in **pkl** format for easy i/o in python.


## Reproducing Results
For the relative directory navigation to the walsh_floquet library to work, scripts need to be run from within their folder.

### Complete Reproduction Pipeline

1. **Environment Setup**: Follow installation instructions above
2. **Data Generation**: Run computational scripts to generate raw data. 
3. **Visualization**: Generate figures using plotting scripts


## Figure Generation

Figures can be regenerated using Python scripts or svg files. Further interactive exploration of results is also available.
Some figures were created using Inkscape, for which the **svg** file is provided.

### Figures
All figures together with scripts are placed at `visual_elements/figs/`.

- **Figure 1** `loc_schematic.svg` (SVG) -> `loc_schematic.pdf`
  **Schematic representation of localisation of the response to a periodic drive in terms of modes in the Fourier and Walsh bases**

- **Figure 2** `walsh_examples.svg` (SVG) -> `walsh_examples.pdf`
  **Example plots of the Walsh functions with corresponding Hadamard matrix and tip of a sawtooth wave expanded in Fourier vs. Walsh basis**

- **Figure 3** `rlscheme.svg` (SVG) -> `rlscheme.pdf`
  **Schematic Reinforcement Learning (RL) feedback loop**

- **Figure 4** `greedy_vs_rl.svg` (SVG) -> `greedy_vs_rl.pdf`
  **Greedy vs. learned active-feedback control strategies**

- **Figure 5** `rl_p.py` (python) -> `rl_p.pdf`
  **Partially informed learned strategies**

- **Figure 6** `example_frame.py` (python) -> `example_frame.pdf`
  **Representative snapshot from supplementary video *learned 128 050.mp4***

- **Figure 7** `bipartite_graph.tex` (TikZ/TeX)
  **Schematic of the composite bipartite system**

- **Figure 8** `transient.py` (python) -> `transient_greedy.pdf` `transient_random.pdf`
  **Ensemble-averaged entanglement Stot(t) dynamics under random/greedy strategy**

- **Figure 9** `full_extended.py` (python) -> `full_extended.pdf`
  **Normalized ensemble-averaged steady-state entanglement ⟨Stot⟩ versus bias parameter *p***

- **Figure 10** `greedy_binder.py` (python) -> `greedy_binder.pdf`
  **Binder cumulant analysis for the greedy strategy**

- **Figure 11** `pyramids_human.py` (python) -> `pyramids_human.pdf`
  **Total entanglement Stot(τ) dynamics under learned and pyramid human-designed strategy**

- **Figure 12** `bell.py` (python) -> `bell.pdf`
  **Bell pair targeting for learned model**

- **Figure 13** `encoder.svg` (SVG) -> `encoder.pdf`
  **Encoder architecture schematic**

- **Figure 14** `stabilizer_diagrams.tex` (TikZ/TeX)
  **All 21 possible local stabilizer structures in a diagrammatic representation**
  

## License Information

### Source Code
- **License**: [BSD 3-Clause License](./LICENSE)
- **Applies to**: All files in `src_code/` directory
- **Permissions**: Commercial use, modification, distribution
- **Requirements**: License and copyright notice

### Data
- **License**: [Creative Commons Attribution 4.0 International (CC-BY 4.0)](./LICENSE-DATA)
- **Applies to**: All files in `data/` directory
- **Permissions**: Share, adapt, commercial use
- **Requirements**: Attribution to original authors


---

**Repository Version**: v1.0
**Last Updated**: TK
**Zenodo DOI**: TK
