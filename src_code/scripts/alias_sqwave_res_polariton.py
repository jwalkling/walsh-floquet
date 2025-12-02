"""
Data for square wave errors and localisation, resonance effects for kicked drive and Walsh polaritons.
=========================================================================

This script produces three processed datasets for L = 1 spin systems:

1. Resonance effects (kick drive other component of spin localisation, Walsh vs Fourier)
   ----------------------------------------------------------------
   System:
       - Single spin (L = 1), driven by an up–down kick.
   Task:
       - For a grid of (h_z, h_x) values, compute the participation entropies
         S_W, S_F of the *spin-down* Floquet mode in Walsh vs Fourier bases.
       - Store the difference:
             Loc_Diff_down = S_F - S_W
         which is positive when the Walsh representation is more localised.

   Output:
       ../../data/processed/spin_kick_loc_down.pkl
   Keys:
       - hxs, hzs
       - Loc_Diff_down

2. Square-wave drive (localisation & quasienergy errors, Walsh vs Fourier)
   -----------------------------------------------------------------------
   System:
       - Single spin with a square wave drive.
   Task:
       - For a grid of (h_z, h_x) values, compute:
           * Photon participation entropies in Walsh vs Fourier extended spaces:
                 PhotonPE_W, PhotonPE_F
           * Their ratio:
                 PhotonPE_ratio = PhotonPE_F / PhotonPE_W
             ( > 1 where Fourier is more delocalised than Walsh )
           * Mean quasienergy error vs exact 2×2 Floquet spectrum:
                 Errors_W, Errors_F
             and their ratio:
                 Errors_ratio = Errors_F / Errors_W

   Output:
       ../../data/processed/spin_square_loc_errors.pkl
   Keys:
       - hxs, hzs
       - Loc_Ratio (PhotonPE_ratio)
       - Errors_Ratio

3. Walsh polariton time evolution
   -------------------------------
   System:
       - Single spin with strong kick drive (h_x ~ π/2).
   Task:
       - Time evolution of an up-component (real) and down-component (imag)
         Floquet mode over one period.
       - The other parts of these complex numbers are approximately zero.
       - This is used to visualise the “Walsh polariton” time evolution.

   Output:
       ../../data/processed/Walsh_polariton_t.pkl
   Keys:
       - ts
       - ut_up_real
       - ut_down_imag
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

# Move into repo/src to import walsh_floquet
p = Path.cwd()
os.chdir(p)
p = p.parent / "src"   # go up one level, then into src
print(p)
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *   # Import the library for Walsh functions


# =========================================================================== #
# PART 1: KICKED DRIVE — SPIN-DOWN LOCALISATION (L = 1)
# =========================================================================== #

# --------------------------------------------------------------------------- #
# System parameters
# --------------------------------------------------------------------------- #

L = 1             # Single-spin system
n = 7
N = 2**n          # Number of Walsh modes

omega = 10
T = 2 * np.pi / omega

steps = 30
hzs = (omega / 2) * np.linspace(0, 1, steps)
hxs = hzs * T / 2


# --------------------------------------------------------------------------- #
# Localisation for spin-down component in both bases
# --------------------------------------------------------------------------- #

LocFd = np.empty((steps, steps))
LocWd = np.empty((steps, steps))

for i, hz in enumerate(hzs):
    for j, hx in enumerate(hxs):
        ut, Wm = compute_Floquet_mode_t_down(n, omega, hz, hx)
        Fm = np.fft.fft(ut) / 2**n

        # Fix log(0) by adding a very small value
        PW = np.abs(Wm)**2 + 10**(-15)
        PF = np.abs(Fm)**2 + 10**(-15)

        # Normalise to interpret as PDFs
        PW /= np.sum(PW)
        PF /= np.sum(PF)

        # Participation entropies
        SW = -np.sum(PW * np.log(PW))
        SF = -np.sum(PF * np.log(PF))

        LocFd[i, j] = SF
        LocWd[i, j] = SW

Loc_Diff_down = LocFd - LocWd
# S higher when delocalised; positive if Walsh is more localised


# --------------------------------------------------------------------------- #
# Save kicked spin-down localisation data
# --------------------------------------------------------------------------- #

filename       = "spin_kick_loc_down.pkl"
processed_path = os.path.join("..", "..", "data", "processed")
save_path      = os.path.join(processed_path, filename)

to_save = {
    "hxs":           hxs,
    "hzs":           hzs,
    "Loc_Diff_down": Loc_Diff_down,
}

os.makedirs(processed_path, exist_ok=True)

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)


# =========================================================================== #
# PART 2: SQUARE-WAVE DRIVE — LOCALISATION & ERROR SCANS (L = 1)
# =========================================================================== #

# --------------------------------------------------------------------------- #
# Truncation parameters & scan grid
# --------------------------------------------------------------------------- #

L = 1                            # Single spin
n = 5
N = 2**n                         # Number of Walsh drive modes
maxmode = int((N - 2) / 2)       # Fourier cutoff index
M = N * 2**L                     # Extended space dimension

omega = 10.0
T = 2 * np.pi / omega

steps = 80
hxs = omega * np.linspace(0, 6, steps)
hzs = hxs.copy()                 # Use same grid for h_z


# --------------------------------------------------------------------------- #
# Square-wave drives in Fourier and Walsh bases
# --------------------------------------------------------------------------- #

drive_F = compute_sqwave_fourier(2 * maxmode + 1)
drive_W = compute_sqwave_walsh(n)


# --------------------------------------------------------------------------- #
# Localisation: photon participation entropy in both bases
# --------------------------------------------------------------------------- #

PhotonPE_F = np.empty((steps, steps))   # Fourier-based drive
PhotonPE_W = np.empty((steps, steps))   # Walsh-based drive

for i, hz in enumerate(hzs):
    for j, hx in enumerate(hxs):

        H0 = hz * PAULI_Z
        H1 = hx * PAULI_X

        H_F = Q_fourier(H0, H1, drive_F, omega, maxmode)
        H_W = Q_walsh(H0, H1, drive_W, omega, n)

        PhotonPE_F[i, j] = compute_PPE(H_F, 2)
        PhotonPE_W[i, j] = compute_PPE(H_W, 2)

PhotonPE_ratio = PhotonPE_F / PhotonPE_W
# > 1 where Fourier is more delocalised than Walsh


# --------------------------------------------------------------------------- #
# Quasienergy errors vs exact 2×2 Floquet spectrum
# --------------------------------------------------------------------------- #

Errors_F = np.empty((steps, steps))
Errors_W = np.empty((steps, steps))

for i, hz in enumerate(hzs):
    for j, hx in enumerate(hxs):

        H0 = hz * PAULI_Z
        H1 = hx * PAULI_X

        H_F = Q_fourier(H0, H1, drive_F, omega, maxmode)
        H_W = Q_walsh(H0, H1, drive_W, omega, n)

        evals_F, _ = np.linalg.eigh(H_F)
        evals_W, _ = np.linalg.eigh(H_W)

        # Exact single-period Floquet phases (2×2 system)
        theta_num = 1.0j * np.log(UF_sqwave(omega, H0, H1))

        # Select the physical block in extended space
        theta_W = evals_W[int(M / 2) - 2**L : int(M / 2)] * T
        theta_F = evals_F[maxmode * 2**L : (maxmode + 1) * 2**L] * T

        Errors_W[i, j] = np.mean(np.abs(theta_W - theta_num))
        Errors_F[i, j] = np.mean(np.abs(theta_F - theta_num))

Errors_ratio = Errors_F / Errors_W   # > 1 where Fourier is worse than Walsh


# --------------------------------------------------------------------------- #
# Save square-wave localisation + error data
# --------------------------------------------------------------------------- #

filename       = "spin_square_loc_errors.pkl"
processed_path = os.path.join("..", "..", "data", "processed")
save_path      = os.path.join(processed_path, filename)

to_save = {
    "hxs":          hxs,
    "hzs":          hzs,
    "Loc_Ratio":    PhotonPE_ratio,
    "Errors_Ratio": Errors_ratio,
}

os.makedirs(processed_path, exist_ok=True)

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)


# =========================================================================== #
# PART 3: WALSH POLARITON TIME EVOLUTION (L = 1, n = 9)
# =========================================================================== #

# --------------------------------------------------------------------------- #
# System parameters
# --------------------------------------------------------------------------- #

L = 1
n = 9
N = 2**n

omega = 100
T = 2 * np.pi / omega

hz = 1.0
hx = np.pi / 2 - 0.01   # Strong kick strength ≫ hz·T


# --------------------------------------------------------------------------- #
# Time evolution of Walsh polariton components
# --------------------------------------------------------------------------- #

ut_up,   _ = compute_Floquet_mode_t_up(n, omega, hz, hx)
ut_down, _ = compute_Floquet_mode_t_down(n, omega, hz, hx)

ts = np.linspace(0.0, T, N, endpoint=False)

ut_up_real = np.real(ut_up)     # spin-up response (imag part ~ 0)
ut_down_imag  = np.imag(ut_down)   # spin-down response (real part ~ 0)


# --------------------------------------------------------------------------- #
# Save Walsh polariton time-evolution data
# --------------------------------------------------------------------------- #

filename       = "walsh_polariton_t.pkl"
raw_path = os.path.join("..", "..", "data", "raw")
save_path      = os.path.join(raw_path, filename)

os.makedirs(raw_path, exist_ok=True)

to_save = {
    "ts":            ts,
    "ut_up_real":  ut_up_real,
    "ut_down_imag":   ut_down_imag,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)

# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #

