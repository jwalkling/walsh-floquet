"""
Time-domain Floquet response and Walsh/Fourier localisation & error scans
========================================================================

This script generates two sets of processed data for a driven single spin (L = 1):

1. Time evolution of the driven Floquet mode (N = 2^9)
   ------------------------------------------------------------------------
   Computes:
     * u(t) over one period for ω and for 0.04·ω
     * Walsh-basis coefficients Wm_high, Wm_low
     * Fourier-basis coefficients Fm_high, Fm_low
   These correspond to the top-right time-series panel and its Walsh/Fourier
   mode decompositions.

   Output:
       ../../data/processed/time_evolve_modes.pkl
       ../../data/raw/ut_low_high_freq.pkl
   containing:
       - ts
       - ut_high_real, ut_low_real
       - Fm_high, Fm_low
       - Wm_high, Wm_low
       - Raw data has ts, ut_low, ut_high

2. Two-parameter localisation and error scans (N = 2^6)
   ---------------------------------------------------
   Computes:
     * Participation entropies S_W, S_F   (localisation)
     * Errors_F − Errors_W in the FFBZ    (Floquet error comparison)

   Output:
       ../../data/processed/spin_kick_loc_errors.pkl
   containing:
       - hxs, hzs
       - Loc_Diff
       - Errors_Diff
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
p = p.parent / "src"
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *    # Import Walsh basis utilities


# =========================================================================== #
# PART 1: TIME EVOLUTION + WALSH/FOURIER MODES (L = 1, n = 9)
# =========================================================================== #

# --------------------------------------------------------------------------- #
# System parameters
# --------------------------------------------------------------------------- #

L = 1
n = 9
N = 2**n

omega = 10
T = 2 * np.pi / omega

hz = 1
hx = np.pi / 2


# --------------------------------------------------------------------------- #
# Time-domain Floquet mode + Walsh/Fourier transforms
# --------------------------------------------------------------------------- #

# High-frequency drive
ut_high, Wm_high = compute_Floquet_mode_t_up(n, omega, hz, hx)

# Low-frequency drive
ut_low, Wm_low   = compute_Floquet_mode_t_up(n, 0.04 * omega, hz, hx)

# Fourier coefficients of the response
Fm_high = np.fft.fft(ut_high) / 2**n
Fm_low  = np.fft.fft(ut_low)  / 2**n

# Time grid
ts = np.linspace(0, T, 2**n)

# Keep only real parts of u(t)
ut_high_real = np.real(ut_high)
ut_low_real  = np.real(ut_low)


# --------------------------------------------------------------------------- #
# Save time-evolution dataset
# --------------------------------------------------------------------------- #

processed_path = os.path.join("..", "..", "data", "raw")
os.makedirs(processed_path, exist_ok=True)

save_path = os.path.join(processed_path, "time_evolve_modes.pkl")

to_save = {
    "ts":            ts,
    "ut_high_real":  ut_high_real,
    "ut_low_real":   ut_low_real,
    "Fm_high":       Fm_high,
    "Fm_low":        Fm_low,
    "Wm_high":       Wm_high,
    "Wm_low":        Wm_low,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

#Save just the full raw data for the modes separately too.
raw_path = os.path.join("..", "..", "data", "raw")
os.makedirs(raw_path, exist_ok=True)

save_path = os.path.join(raw_path, "ut_low_high_freq.pkl")

to_save = {
    "ts":            ts,
    "ut_high":  ut_high,
    "ut_low":   ut_low,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)


# =========================================================================== #
# PART 2: TWO-PARAMETER LOCALISATION & FLOQUET ERROR SCANS (L = 1)
# =========================================================================== #

# --------------------------------------------------------------------------- #
# System parameters
# --------------------------------------------------------------------------- #

L = 1
n = 7
N = 2**n

omega = 10
T = 2 * np.pi / omega

steps = 30
hzs = (omega / 2) * np.linspace(0, 1, steps)
hxs = hzs * T / 2


# --------------------------------------------------------------------------- #
# Localisation scan: Walsh vs Fourier participation entropy
# --------------------------------------------------------------------------- #

Loc_F_up = np.empty((steps, steps))
Loc_W_up = np.empty((steps, steps))

for i, hz in enumerate(hzs):
    for j, hx in enumerate(hxs):

        ut, Wm = compute_Floquet_mode_t_up(n, omega, hz, hx)
        Fm = np.fft.fft(ut) / 2**n

        # Prevent log(0)
        PW = np.abs(Wm)**2 + 1e-15
        PF = np.abs(Fm)**2 + 1e-15

        PW /= np.sum(PW)
        PF /= np.sum(PF)

        SW = -np.sum(PW * np.log(PW))
        SF = -np.sum(PF * np.log(PF))

        Loc_W_up[i, j] = SW
        Loc_F_up[i, j] = SF

Loc_Diff = Loc_F_up - Loc_W_up   # Positive means Fourier is more delocalised


# --------------------------------------------------------------------------- #
# Error scan: Walsh vs Fourier quasienergy errors
# --------------------------------------------------------------------------- #

L = 1
n = 6
N = 2**n
maxmode = int((N - 2) / 2)
M = 2**L * N

omega = 10
T = 2 * np.pi / omega

steps = 30
hzs = (omega / 2) * np.linspace(0, 1, steps)
hxs = hzs * T / 2

drive_fourier   = omega / np.pi * compute_kick_fourier(2 * maxmode + 1)
drive_walsh = omega / np.pi * compute_kick_walsh(n)

Errors_F = np.empty((steps, steps))
Errors_W = np.empty((steps, steps))

for i, hz in enumerate(hzs):
    for j, hx in enumerate(hxs):

        H0 = hz * PAULI_Z
        H1 = hx * PAULI_X

        H_Fourier = Q_fourier(H0, H1, drive_fourier, omega, maxmode)
        H_Walsh   = Q_walsh(H0, H1, drive_walsh, omega, n)

        evals_F = np.linalg.eigh(H_Fourier)[0]
        evals_W = np.linalg.eigh(H_Walsh)[0]

        numsol = 1.j * np.log(UF_kick(omega, H0, H1))

        # Relevant Floquet Brillouin zone pieces
        evals_W_FFBZ = evals_W[int(M/2) - 2**L : int(M/2)] * T
        evals_F_FFBZ = evals_F[maxmode*2**L : (maxmode+1)*2**L] * T

        Errors_W[i, j] = np.average(np.abs(evals_W_FFBZ - numsol))
        Errors_F[i, j] = np.average(np.abs(evals_F_FFBZ - numsol))

Errors_Diff = Errors_F - Errors_W


# --------------------------------------------------------------------------- #
# Save localisation + error scan dataset
# --------------------------------------------------------------------------- #

save_path = os.path.join(processed_path, "spin_kick_loc_errors.pkl")

to_save = {
    "hxs":        hxs,
    "hzs":        hzs,
    "Loc_Diff":   Loc_Diff,
    "Errors_Diff": Errors_Diff,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)

# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #
