# Walsh-Floquet Functions

This repository provides a modular and efficient collection of Walsh–Floquet utilities for analysing periodically driven quantum systems using the Walsh function basis. The codebase provides efficient implementations of Walsh–Hadamard transforms and sequency-ordered Walsh functions, together with tools for constructing Walsh expansions of time-dependent Hamiltonians. It also includes utilities for computing Floquet operators, effective Hamiltonians, and kicked dynamics, as well as helper routines for visualisation, benchmarking, and fully reproducible numerical experiments.

---

## 📑 Table of Contents

- [Overview](#overview)
- [Code Structure](#code-structure)
- [Source Code Overview](#source-code-overview)
  - [Key Functions](#key-functions)
- [Setup and Execution](#setup-and-execution)

---

## Overview

This code is designed to allow for easy computation of comparing exact diagonalisation of the quasienergy matrix in extended space with the exact numerical solutions. There is also additional functionality to calculate localisation on the frequency lattice.

---

## Code Structure

The repository is organised as below

```plaintext
.
├── scripts/       # Scripts for data generation of figures.
├── src/            # Source code for the extended Floquet Walsh operators.
└── CODE_STRUCTURE.txt      # High-level structure of the codebase.
```

### Directories and Main Files

1. **[scripts/](scripts/)**:
   - Scripts for generating the data for the figures.
     - `spin_kick.py` and `mb_spin_kick.py`: Calculates quasienergy spectrum and error due to truncation for kicked spins.
     - `localisation_vs_error.py`: Scanning parameters to find error vs. localisation for a kicked spin.
     - `error_scaling.py` and `error_scaling_strong_kick.py`: Scaling of error with truncation order for square and kick drive.
     - `quasienergy_matrices.py` and `discrete_generator`: matrix elements for the quasienergy operator in different bases.
     - `walsh_series.py`: scaling of analytic Walsh expansion vs. numerics with driving frequency.
     - `alias_sqwave_res_polariton.py`: square wave error vs. localisation, resonance effects and Walsh polariton

2. **[src/](src/)**:
   - Core source directory for the Walsh-Floquet utilities.

---

## Source Code Overview

The following key functions are important for implementing the Walsh-Floquet formalism:

### Key Functions

--- **walsh_generator**: generator of time-translations in the Walsh basis
--- **Q_walsh**: quasienergy operator in the Walsh basis
--- **compute_PPE**: calculate the photon participation entropy (measure of localisation traced over physical degrees of freedom).

## Setup and Execution


### 1. Install Dependencies
Create a virtual environment (recommended) and install dependencies.

```bash
conda env create -f environment.yml
```
### 2. Run the scripts from inside the folder
The code is set-up with relative pathing. 

---

