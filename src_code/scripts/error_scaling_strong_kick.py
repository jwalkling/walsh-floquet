"""
Strong-Kick Floquet Error Scaling: Walsh vs. Fourier Truncation
=================================================================

This script computes the scaling of quasi-energy errors for a *strong*
up–down kick drive in both the **single-particle** and **many-body (L=3)**
cases. The comparison is performed between:

    • Fourier-truncated extended space (F)
    • Walsh-truncated extended space (W)

For each truncation size `n`, the script constructs the extended Hamiltonians
(Q_W and Q_F), extracts the physical photon block of the quasi-energy spectrum,
and compares the result to the *exact (untruncated)* Floquet spectrum.

Two cases are computed:

1. Single-particle, π/2 kick
   --------------------------------------------------
   Hamiltonian:
       H0 = Z
       H1 = (π/2) X
   Drive:
       Strong up–down kick with frequency ω = 10

   Output arrays:
       • spstrongkick_F   – Fourier truncation errors vs n
       • spstrongkick_W   – Walsh truncation errors vs n

2. Many-body (L=3), π/4 × 1.1 transverse-field kick
   --------------------------------------------------
   Hamiltonian:
       H0 = J Σ ZZ + hz Σ Z
       H1 = hx Σ X      with hx = 1.1 × (π/4)
   Drive:
       Strong up–down kick with frequency ω = 10

   Output arrays:
       • mbstrongkick_F   – Fourier truncation errors vs n
       • mbstrongkick_W   – Walsh truncation errors vs n


Definitions
-----------
For each value of n:

    N = 2^n                     # number of Walsh modes
    maxmode = (N-2)/2           # Fourier cutoff
    d = 2^L                     # physical Hilbert-space dimension

The physical photon block of the extended Hamiltonian spectrum is extracted as:

    • Fourier:  eigenvalues[d*maxmode : d*maxmode + d]
    • Walsh:    eigenvalues[(N//2 - 1)*d : (N//2)*d]

All errors computed are *mean absolute quasi-energy differences*:

    error_F = mean(|E_F(n) - E_exact|)
    error_W = mean(|E_W(n) - E_exact|)


Output
------
The script saves one pickle file:

    ../../data/processed/strong_kick_error_scaling.pkl

containing:

    {
        "nvals":            array of single-particle truncation sizes
        "nvalsmb":          array of many-body truncation sizes
        "spstrongkick_F":   single-particle Fourier errors
        "spstrongkick_W":   single-particle Walsh errors
        "mbstrongkick_F":   many-body Fourier errors
        "mbstrongkick_W":   many-body Walsh errors
    }

"""


import os
import sys
import pickle
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------
# Set working directory and import Walsh helper functions
# ---------------------------------------------------------------------
p = Path.cwd()
os.chdir(p)
p = p.parent / "src"   # go up one level, then into src
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *  # Import the library for Walsh functions

# ---------------------------------------------------------------------
# Common n-values
# ---------------------------------------------------------------------
# Single-particle truncations
nvals = np.arange(2, 11, 1)      # n = 2,...,10
# Many-body truncations for L=3
nvalsmb = np.arange(3, 9, 1)    # n = 3,...,10


# =====================================================================
# 1. Single-particle, up–down kick drive (π/2 kick)
# =====================================================================

H0 = 1 * PAULI_Z
H1 = (np.pi / 2) * PAULI_X
omega = 10.0
T = 2 * np.pi / omega

# Untruncated Floquet spectrum
numsol = 1.0j * np.log(UF_kick(omega, H0, H1)) / T

EF_errors = []
EW_errors = []

for n in nvals:
    # Setup the matrices
    N = 2 ** n
    maxmode = int((N - 2) / 2)

    drive_fourier = (omega / np.pi) * compute_kick_fourier(2 * maxmode + 1)
    drive_walsh = (omega / np.pi) * compute_kick_walsh(n)

    # Extended Hamiltonians and physical blocks
    EF = np.linalg.eigh(
        Q_fourier(H0, H1, drive_fourier, omega, maxmode)
    )[0][2 * maxmode : 2 * maxmode + 2]

    EW = np.linalg.eigh(
        Q_walsh(H0, H1, drive_walsh, omega, n)
    )[0][N - 2 : N]


    # (Debug line – not used further)
    ErrorW = EW - Eexact
    _ = np.average(np.abs(ErrorW))

    # Take the average of the absolute values of EF, EW, and Eexact
    EF_avg = np.mean(np.abs(EF))
    EW_avg = np.mean(np.abs(EW))
    numsol_avg = np.mean(np.abs(numsol))

    # Errors |EF_avg - Eexact_avg| and |EW_avg - Eexact_avg|
    EF_error = np.abs(EF_avg - numsol_avg)
    EW_error = np.abs(EW_avg - numsol_avg)

    # Append the errors for plotting
    EF_errors.append(EF_error)
    EW_errors.append(EW_error)

# Convert lists to numpy arrays for easier plotting
spstrongkick_F = np.array(EF_errors)
spstrongkick_W = np.array(EW_errors)


# =====================================================================
# 2. Many-body (L=3), up–down kick drive (π/4 * 1.1 transverse field)
# =====================================================================

# Basic parameters
nvalsmb = np.arange(3, 9, 1)
L = 3
omega = 10.0
T = 2 * np.pi / omega

# Further tuning
J = 1.0
hz = omega / 2 + 0.5
hx = (np.pi / 4) * 1.1

# Hamiltonians
H0 = H_ising_zz_zfield(J, hz, L)
H1 = H_ising_xfield(hx, L)

# Exact many-body Floquet spectrum
numsol = 1.0j * np.log(UF_kick(omega, H0, H1)) / T #Numerical solution (untruncated)

EF_errorsmb = []
EW_errorsmb = []

for n in nvalsmb:
    # Setup the matrices
    N = 2 ** n
    d = 2 ** L
    halfW = int((N - 1) / 2)
    maxmode = int((N - 2) / 2)

    drive_fourier = (omega / np.pi) * compute_kick_fourier(2 * maxmode + 1)
    drive_walsh = (omega / np.pi) * compute_kick_walsh(n)

    # Calculate each solution for the energies
    EW = np.linalg.eigh(
        Q_walsh(H0, H1, drive_walsh, omega, n)
    )[0][halfW * d : (halfW + 1) * d]

    EF = np.linalg.eigh(
        Q_fourier(H0, H1, drive_fourier, omega, maxmode)
    )[0][d * maxmode : d * maxmode + d]

    # Mean error per eigenstate
    EF_errormb = np.sum(np.abs(EF - numsol)) / (2 ** L)
    EW_errormb = np.sum(np.abs(EW - numsol)) / (2 ** L)

    # Append the errors for plotting
    EF_errorsmb.append(EF_errormb)
    EW_errorsmb.append(EW_errormb)

# Convert lists to numpy arrays for easier plotting
mbstrongkick_F = np.array(EF_errorsmb)
mbstrongkick_W= np.array(EW_errorsmb)


# =====================================================================
# SAVE RESULTS
# =====================================================================
filename = "strong_kick_error_scaling.pkl"
processed_path = Path("..") / ".." / "data" / "processed"
os.makedirs(processed_path, exist_ok=True)
save_path = processed_path / filename

to_save = {
    "nvals": nvals,
    "nvalsmb": nvalsmb,
    "spstrongkick_F": spstrongkick_F,
    "spstrongkick_W": spstrongkick_W,
    "mbstrongkick_F": mbstrongkick_F,
    "mbstrongkick_W": mbstrongkick_W,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)





# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #
