"""
Error scaling with truncation order n for square-wave and kick drives
in both Fourier and Walsh bases, single-particle and many-body.
=================================================================

This script calculates the error scaling relative to the exact numerical solution.
For both the Walsh basis and the Fourier basis it exports:

- "nvals"    : np.ndarray of single-particle truncation n values (np.arange(2,11))
- "mbkick_F" : EFmbkick — many-body (L=3) kick-drive Fourier error vs nvalsmb
- "spkick_F" : EFspkick — single-particle kick-drive Fourier error vs nvals
- "mbsq_F"   : EFmbsq  — many-body (L=3) square-wave Fourier error vs nvalsmb
- "spsq_F"   : EFspsq  — single-particle square-wave Fourier error vs nvals
- "mbkick_W" : EWmbkick — many-body (L=3) kick-drive Walsh error vs nvalsmb
- "spkick_W" : EWspkick — single-particle kick-drive Walsh error vs nvals
- "mbsq_W"   : EWmbsq  — many-body (L=3) square-wave Walsh error vs nvalsmb
- "spsq_W"   : EWspsq  — single-particle square-wave Walsh error vs nvals

Definitions
-----------
Relevant parameter choices used in this script are listed next to their use
in the main body (nvals, nvalsmb, drive frequencies and Hamiltonians).

Output
------
The script saves a pickle file:

    ../../data/processed/square_kick_error_scaling.pkl

"""


import os
import sys
import pickle
from pathlib import Path

import numpy as np

# =====================================================================
# Working directory and imports
# =====================================================================

p = Path.cwd()
os.chdir(p)
p = p.parent / "src"          # move into repo/src
print("Changing to:", p)
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *    # Walsh/Fourier helpers


# =====================================================================
# Parameters: n-values (single-particle and many-body)
# =====================================================================

# Single-particle truncations
nvals = np.arange(2, 11, 1)          # n = 2,...,10

# Many-body truncations for L = 3
nvalsmb = np.arange(3, 10, 1)        # n = 3,...,10


# =====================================================================
# 1. Single-particle, square-wave drive
#    UF_sqwave  →  EFspsq, EWspsq
# =====================================================================

omega = 50.0
T = 2 * np.pi / omega

H0 = PAULI_Z
H1 = 6 * PAULI_X

EF_errorssq = []
EW_errorssq = []

for n in nvals:
    N = 2 ** n
    maxmode = int((N - 2) / 2)

    drive_fourier = compute_sqwave_fourier(2 * maxmode + 1)
    drive_walsh   = compute_sqwave_walsh(n)

    # Extended Hamiltonians
    H_F = Q_fourier(H0, H1, drive_fourier, omega, maxmode)
    H_W = Q_walsh  (H0, H1, drive_walsh,   omega, n)

    # Physical 2×2 photon block
    EF = np.linalg.eigh(H_F)[0][2 * maxmode : 2 * maxmode + 2]
    EW = np.linalg.eigh(H_W)[0][N - 2 : N]

    # Exact spectrum
    Eexact = 1.0j * np.log(UF_sqwave(omega, H0, H1)) / T

    # Mean absolute energy difference
    EF_avg = np.mean(np.abs(EF))
    EW_avg = np.mean(np.abs(EW))
    Eex_avg = np.mean(np.abs(Eexact))

    EF_errorssq.append(np.abs(EF_avg - Eex_avg))
    EW_errorssq.append(np.abs(EW_avg - Eex_avg))

EFspsq = np.array(EF_errorssq)
EWspsq = np.array(EW_errorssq)


# =====================================================================
# 2. Single-particle, up–down kick drive
#    UF_kick  →  EFspkick, EWspkick
# =====================================================================

omega = 10.0
T = 2 * np.pi / omega

H0 = (omega / 2 + 0.5) * PAULI_Z
H1 = 0.2 * PAULI_X

EF_errors = []
EW_errors = []

for n in nvals:
    N = 2 ** n
    maxmode = int((N - 2) / 2)

    drive_fourier = (omega / np.pi) * compute_kick_fourier(2 * maxmode + 1)
    drive_walsh   = (omega / np.pi) * compute_kick_walsh(n)

    H_F = Q_fourier(H0, H1, drive_fourier, omega, maxmode)
    H_W = Q_walsh  (H0, H1, drive_walsh,   omega, n)

    EF = np.linalg.eigh(H_F)[0][2 * maxmode : 2 * maxmode + 2]
    EW = np.linalg.eigh(H_W)[0][N - 2 : N]

    Eexact = 1.0j * np.log(UF_kick(omega, H0, H1)) / T

    EF_avg = np.mean(np.abs(EF))
    EW_avg = np.mean(np.abs(EW))
    Eex_avg = np.mean(np.abs(Eexact))

    EF_errors.append(np.abs(EF_avg - Eex_avg))
    EW_errors.append(np.abs(EW_avg - Eex_avg))

EFspkick = np.array(EF_errors)
EWspkick = np.array(EW_errors)


# =====================================================================
# 3. Many-body (L=3), square-wave drive
#    UF_sqwave  →  EFmbsq, EWmbsq
# =====================================================================

L = 3
omega = 50.0
T = 2 * np.pi / omega

J = 1.0
hz = 1.0
hx = 6.0

H0 = H_ising_zz_zfield(J, hz, L)
H1 = H_ising_xfield(hx, L)

Eexact = 1.0j * np.log(UF_sqwave(omega, H0, H1)) / T

EF_errorsssqmb = []
EW_errorsssqmb = []

for n in nvalsmb:
    N = 2 ** n
    d = 2 ** L
    halfW = int((N - 1) / 2)
    maxmode = int((N - 2) / 2)

    drive_fourier = compute_sqwave_fourier(2 * maxmode + 1)
    drive_walsh   = compute_sqwave_walsh(n)

    H_W = Q_walsh  (H0, H1, drive_walsh,   omega, n)
    H_F = Q_fourier(H0, H1, drive_fourier, omega, maxmode)

    EW = np.linalg.eigh(H_W)[0][halfW * d : (halfW + 1) * d]
    EF = np.linalg.eigh(H_F)[0][d * maxmode : d * maxmode + d]

    EF_errorsssqmb.append(np.sum(np.abs(EF - Eexact)) / d)
    EW_errorsssqmb.append(np.sum(np.abs(EW - Eexact)) / d)

EFmbsq = np.array(EF_errorsssqmb)
EWmbsq = np.array(EW_errorsssqmb)


# =====================================================================
# 4. Many-body (L=3), up–down kick drive
#    UF_kick  →  EFmbkick, EWmbkick
# =====================================================================

L = 3
omega = 10.0
T = 2 * np.pi / omega

J = 1.0
hz = omega / 2 + 0.5
hx = 0.2

H0 = H_ising_zz_zfield(J, hz, L)
H1 = H_ising_xfield(hx, L)

Eexact = 1.0j * np.log(UF_kick(omega, H0, H1)) / T

EF_errorsmb = []
EW_errorsmb = []

for n in nvalsmb:
    N = 2 ** n
    d = 2 ** L
    halfW = int((N - 1) / 2)
    maxmode = int((N - 2) / 2)

    drive_fourier = (omega / np.pi) * compute_kick_fourier(2 * maxmode + 1)
    drive_walsh   = (omega / np.pi) * compute_kick_walsh(n)

    H_W = Q_walsh  (H0, H1, drive_walsh,   omega, n)
    H_F = Q_fourier(H0, H1, drive_fourier, omega, maxmode)

    EW = np.linalg.eigh(H_W)[0][halfW * d : (halfW + 1) * d]
    EF = np.linalg.eigh(H_F)[0][d * maxmode : d * maxmode + d]

    EF_errorsmb.append(np.sum(np.abs(EF - Eexact)) / d)
    EW_errorsmb.append(np.sum(np.abs(EW - Eexact)) / d)

EFmbkick = np.array(EF_errorsmb)
EWmbkick = np.array(EW_errorsmb)


# =====================================================================
# SAVE RESULTS
# =====================================================================

filename = "square_kick_error_scaling.pkl"
processed_path = Path("..") / ".." / "data" / "processed"
os.makedirs(processed_path, exist_ok=True)

save_path = processed_path / filename

to_save = {
    "nvals":    nvals,
    "mbkick_F": EFmbkick,
    "spkick_F": EFspkick,
    "mbsq_F":   EFmbsq,
    "spsq_F":   EFspsq,
    "mbkick_W": EWmbkick,
    "spkick_W": EWspkick,
    "mbsq_W":   EWmbsq,
    "spsq_W":   EWspsq,
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)


# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #
