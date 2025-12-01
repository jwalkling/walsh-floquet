"""
Generate Floquet spectral and error data for multi-spin up–down kick drives
==========================================================================

This script produces two sets of processed data for the up–down kick drive for L spins:

1. L = 3 single-parameter sweep (spectrum)
   -----------------------------------------------------------------------
   - System: 3-spin Ising chain (nearest-neighbour J, transverse x-field hx,
     longitudinal z-field hz).
   - Drive: up–down kick, implemented in both Walsh and Fourier extended spaces.
   - Task: For fixed J and hz, sweep over hx and compute:
       * Floquet quasienergies in the Walsh basis extended Hilbert space.
       * Floquet quasienergies in the Fourier basis extended Hilbert space.
       * Numerically exact 2^L × 2^L Floquet spectrum from UF_kick.
   - Output file:
       ../../data/processed/3_spin_kick_spectrum.pkl
     containing:
       * hxs_W, hxs_F, hxs_num
       * evals_W_flat, evals_F_flat, numsols_flat

2. L = 6 two-parameter error scan for J = 1 and J = 3
   ---------------------------------------------------
   - System: 6-spin Ising chain (nearest-neighbour J, transverse hx, longitudinal hz).
   - Drive: up–down kick in Walsh and Fourier bases via compute_MFIM_kick_errors.
   - Task: Over a grid of (hz, hx) values, compute:
       * Fourier-basis Floquet errors Errors_F_J1, Errors_F_J3
       * Walsh-basis Floquet errors Errors_W_J1, Errors_W_J3
       * Log10 error ratios log10(Errors_F / Errors_W) for J = 1 and J = 3.
   - Output file:
       ../../data/processed/6_spin_kick_errors.pkl
     containing:
       * hxs, hzs
       * logErrorJ1, logErrorJ3
       * Errors_W_J1, Errors_W_J3
"""

# --------------------------------------------------------------------------- #
# Imports
# --------------------------------------------------------------------------- #

import os
import sys
import pickle
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors

# Move into src/ directory to import walsh_floquet
p = Path.cwd()
p = p.parent / "src"   # go up one level, then into src
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *   # Import the library for Walsh functions


# =========================================================================== #
# PART 1: L = 3, SINGLE-SLICE SPECTRUM VS hx
# =========================================================================== #

# --------------------------------------------------------------------------- #
# System parameters
# --------------------------------------------------------------------------- #

L = 3                      # Number of spins
n = 5                      # Walsh order
N = 2 ** n                 # Number of Walsh modes
maxmode = int((N - 2) / 2) # Fourier cutoff

omega = 10
T = 2 * np.pi / omega

J = 1
hz = omega / 2 + 0.5       # Fixed static field for spectrum slice

# Sweep parameter
hxs = np.linspace(0, np.pi / 2, 600)


# --------------------------------------------------------------------------- #
# Drive matrices
# --------------------------------------------------------------------------- #

drive_walsh = omega / np.pi * compute_kick_walsh(n)
drive_fourier   = omega / np.pi * compute_kick_fourier(2 * maxmode + 1)


# --------------------------------------------------------------------------- #
# Compute spectra
# --------------------------------------------------------------------------- #

evals_W = []     # Walsh-basis extended-space eigenvalues
evals_F = []     # Fourier-basis extended-space eigenvalues
numsols = []     # Exact Floquet eigenvalues from UF_kick

for hx in hxs:
    H0 = H_ising_zz_zfield(J, hz, L)
    H1 = H_ising_xfield(hx, L)

    eval_W = np.linalg.eigh(Q_walsh(H0, H1, drive_walsh, omega, n))[0]
    eval_F = np.linalg.eigh(Q_fourier(H0, H1, drive_fourier, omega, maxmode))[0]
    numsol = 1.j * np.log(UF_kick(omega, H0, H1)) * omega / (2 * np.pi)

    evals_W.append(eval_W)
    evals_F.append(eval_F)
    numsols.append(numsol)


# --------------------------------------------------------------------------- #
# Flatten and arrange for saving
# --------------------------------------------------------------------------- #

evals_W_flat = np.array(evals_W).flatten()
evals_F_flat = np.array(evals_F).flatten()
numsols_flat = np.array(numsols).flatten()

hxs_W   = np.repeat(hxs, 2 ** L * N)
hxs_F   = np.repeat(hxs, 2 ** L * (2 * maxmode + 1))
hxs_num = np.repeat(hxs, 2 ** L)


# --------------------------------------------------------------------------- #
# Save L = 3 spectrum data
# --------------------------------------------------------------------------- #

processed_path = os.path.join("..", "..", "data", "processed")
os.makedirs(processed_path, exist_ok=True)

filename_3spin = "3_spin_kick_spectrum.pkl"
save_path = os.path.join(processed_path, filename_3spin)

to_save = {
    "hxs_W":        hxs_W,
    "hxs_F":        hxs_F,
    "hxs_num":      hxs_num,
    "evals_W_flat": evals_W_flat,
    "evals_F_flat": evals_F_flat,
    "numsols_flat": numsols_flat,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)


# =========================================================================== #
# PART 2: L = 6, TWO-PARAMETER ERROR SCAN (J = 1, 3)
# =========================================================================== #

# --------------------------------------------------------------------------- #
# Scan parameters and system setup
# --------------------------------------------------------------------------- #

# Note: This data is expensive to generate:
#   At J = 1, 3, numvals = 50:
#   ~1 minute per data point ⇒ ~2500 minutes in total.

omega = 10
T = 2 * np.pi / omega

L = 6        # System size 6 spins
n = 6        # Number of Walsh modes
numvals = 50 # Resolution of scan grid 50

hzs = omega * np.linspace(0, 1, numvals)
hxs = hzs * T / 4


# --------------------------------------------------------------------------- #
# Compute Walsh vs Fourier errors for J = 1 and J = 3
# --------------------------------------------------------------------------- #

# J = 1
Errors_F_J1, Errors_W_J1 = compute_MFIM_kick_errors(omega, 1, hzs, hxs, L, n)

# J = 3
Errors_F_J3, Errors_W_J3 = compute_MFIM_kick_errors(omega, 3, hzs, hxs, L, n)

# Log10 error ratios (Fourier / Walsh)
log_valsJ1 = np.log10(Errors_F_J1 / Errors_W_J1)
log_valsJ3 = np.log10(Errors_F_J3 / Errors_W_J3)


# --------------------------------------------------------------------------- #
# Save L = 6 error scan data
# --------------------------------------------------------------------------- #

filename_6spin = "6_spin_kick_errors.pkl"
save_path = os.path.join(processed_path, filename_6spin)

to_save = {
    "hxs":         hxs,
    "hzs":         hzs,
    "logErrorJ1":  log_valsJ1,
    "logErrorJ3":  log_valsJ3,
    "Errors_W_J1": Errors_W_J1,
    "Errors_W_J3": Errors_W_J3,
}

os.makedirs(processed_path, exist_ok=True)

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)

# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #
