"""
Extended-space drive and derivative matrices for the up–down kick in Fourier vs. Walsh
=================================================================

This script constructs the extended-Hilbert-space operators to be plotted schematically.
For both the Walsh basis and the Fourier basis it exports:

    • drive_W :   up–down kick drive operator in the Walsh basis
    • drive_F :   up–down kick drive operator in the Fourier basis
    • deriv_W :   time-translation generator in the Walsh basis
    • deriv_F :   time-translation generator in the Fourier basis

Definitions
-----------
The extended Hilbert space for a single spin is:

    H_ext = H_spin ⊗ L_drive

where H_drive is either the Walsh mode space (dimension 2^n)
or the Fourier truncated space (dimension 2*maxmode + 1).

For the up–down kick drive:

    drive_W  = (ω/π) compute_kick_walsh(n)
    drive_F  = (ω/π) compute_kick_fourier(2*maxmode+1)

For the time-translation generator:

    deriv_W = ω * D_W
    deriv_F = ω * diag(-maxmode ... maxmode)

Output
------
The script saves a pickle file:

    ../../data/processed/H_extended_kick.pkl

containing the four matrices above.
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

# Move into repo/src directory to access walsh_floquet
p = Path.cwd()
os.chdir(p)
p = p.parent / "src"
print(p)
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *     # provides compute_kick_walsh, walsh_generator, etc.


# --------------------------------------------------------------------------- #
# Constants and setup
# --------------------------------------------------------------------------- #

omega = 10
T = 2 * np.pi / omega

n = 6
N = 2 ** n                      # number of Walsh modes
maxmode = int((N - 2) / 2)      # Fourier cutoff
L = 1                           # single spin (for completeness)

# Parameters of the physical spin sector (not used directly here)
hz = omega / 2 - 0.5
hx = 5


# --------------------------------------------------------------------------- #
# Extended-space drive matrices (Walsh and Fourier bases)
# --------------------------------------------------------------------------- #

drive_W = (omega / np.pi) * compute_kick_walsh(n)
drive_F = (omega / np.pi) * compute_kick_fourier(2 * maxmode + 1)


# --------------------------------------------------------------------------- #
# Extended-space derivative (time-translation) generators
# --------------------------------------------------------------------------- #

# Walsh basis derivative operator
deriv_W = omega * walsh_generator(N)

# Fourier basis time-translation generator: diagonal with integer weights
deriv_F = omega * np.diag(np.arange(-maxmode, maxmode + 1, 1))


# --------------------------------------------------------------------------- #
# Save results to processed/ directory
# --------------------------------------------------------------------------- #

filename       = "H_extended_kick.pkl"
processed_path = os.path.join("..", "..", "data", "processed")
save_path      = os.path.join(processed_path, filename)

os.makedirs(processed_path, exist_ok=True)

to_save = {
    "drive_W": drive_W,     # Walsh up–down kick drive
    "drive_F": drive_F,     # Fourier up–down kick drive
    "deriv_W": deriv_W,     # Walsh time-translation generator
    "deriv_F": deriv_F,     # Fourier time-translation generator
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)

# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #



