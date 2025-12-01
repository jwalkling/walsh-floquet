"""
Generate Floquet spectral data and Walsh/Fourier error comparison
=================================================================

For the case of an up–down kick drive on a single spin-1/2 (L=1), this script computes:
1. The Floquet quasienergies of a single-spin up–down kick drive,
   evaluated in both the Walsh and Fourier extended Hilbert spaces.

2. A dense two-parameter scan (hx, hz) comparing the Walsh-basis
   and Fourier-basis Floquet errors relative to the numerically exact solution.

All results are written into the repository-standard directory:
    ../../data/processed/

The output files are:
    - spin_kick_spectrum.pkl
    - spin_kick_errors.pkl
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
p = p.parent / "src"
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *   # Import Walsh basis utilities


# --------------------------------------------------------------------------- #
# Constants and Setup
# --------------------------------------------------------------------------- #

omega = 10                           # Drive frequency
T = 2 * np.pi / omega                # Drive period

n = 5                                # Walsh order (N = 2^n)
N = 2 ** n
maxmode = int((N - 2) / 2)           # Fourier truncation
L = 1                                # One physical spin

hz_fixed = omega / 2 - 0.5           # Fixed hz for spectral slice

# Parameter range for spectrum slice
hxs = np.linspace(0, np.pi / 2, 600)


# --------------------------------------------------------------------------- #
# Drive matrices in both bases
# --------------------------------------------------------------------------- #

drive_walsh = omega / np.pi * compute_kick_walsh(n)
drive_fourier   = omega / np.pi * compute_kick_fourier(2 * maxmode + 1)


# --------------------------------------------------------------------------- #
# Floquet Spectrum at Fixed hz
# --------------------------------------------------------------------------- #

evals_W = []        # Walsh eigenvalues
evals_F = []        # Fourier eigenvalues
numsols = []        # Exact numerical solution

for hx in hxs:

    H0 = hz_fixed * PAULI_Z
    H1 = hx * PAULI_X

    # Extended-space eigenvalues
    eval_W = np.linalg.eigh(Q_walsh(H0, H1, drive_walsh, omega, n))[0]
    eval_F = np.linalg.eigh(Q_fourier(H0, H1, drive_fourier, omega, maxmode))[0]

    # Exact 2×2 Floquet spectrum
    numsol = 1.j * np.log(UF_kick(omega, H0, H1)) * omega / (2 * np.pi)

    evals_W.append(eval_W)
    evals_F.append(eval_F)
    numsols.append(numsol)

# Flatten arrays for saving/plotting
evals_W_flat = np.array(evals_W).flatten()
evals_F_flat = np.array(evals_F).flatten()
numsols_flat = np.array(numsols).flatten()

hxs_W   = np.repeat(hxs, 2 ** L * N)
hxs_F   = np.repeat(hxs, 2 ** L * (2 * maxmode + 1))
hxs_num = np.repeat(hxs, 2 ** L)


# --------------------------------------------------------------------------- #
# Save 1D Spectral Slice
# --------------------------------------------------------------------------- #

processed_path = os.path.join("..", "..", "data", "processed")
os.makedirs(processed_path, exist_ok=True)

save_path = os.path.join(processed_path, "spin_kick_spectrum.pkl")

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


# --------------------------------------------------------------------------- #
# Two-Parameter Error Scan (hx, hz)
# --------------------------------------------------------------------------- #

numvals = 80

hxs = np.linspace(0, np.pi / 2, numvals)
hzs = omega * np.linspace(0.01, 0.5, numvals)

M = 2 ** L * N     # Extended space dimension

Errors_F = np.empty((numvals, numvals))
Errors_W = np.empty((numvals, numvals))

for i, hz in enumerate(hzs):
    for j, hx in enumerate(hxs):

        H0 = hz * PAULI_Z
        H1 = hx * PAULI_X

        # Extended Hamiltonians
        H_Fourier = Q_fourier(H0, H1, drive_fourier, omega, maxmode)
        H_Walsh   = Q_walsh(H0, H1, drive_walsh, omega, n)

        sol_F = np.linalg.eigh(H_Fourier)[0]
        sol_W = np.linalg.eigh(H_Walsh)[0]

        # Exact two-level Floquet result
        exactsol = 1.j * np.log(UF_kick(omega, H0, H1))

        # Extract relevant Floquet bands
        E_W = sol_W[int(M / 2) - 2 ** L : int(M / 2)] * T
        E_F = sol_F[maxmode * 2 ** L : (maxmode + 1) * 2 ** L] * T

        # Errors
        Errors_W[i, j] = np.average(np.abs(E_W - exactsol))
        Errors_F[i, j] = np.average(np.abs(E_F - exactsol))


# --------------------------------------------------------------------------- #
# Compute log10(Errors_F / Errors_W)
# --------------------------------------------------------------------------- #

vals     = np.divide(Errors_F, Errors_W, out=np.zeros_like(Errors_F), where=(Errors_W != 0))
log_vals = np.zeros_like(vals)

mask = vals > 0
log_vals[mask] = np.log10(vals[mask])


# --------------------------------------------------------------------------- #
# Save 2D Error Data
# --------------------------------------------------------------------------- #

save_path = os.path.join(processed_path, "spin_kick_errors.pkl")

to_save = {
    "hxs":      hxs,
    "hzs":      hzs,
    "logError": log_vals,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)

# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #





