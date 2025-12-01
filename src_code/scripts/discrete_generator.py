"""
Generator of discrete time translations (real-space derivative operator)
========================================================================

This script constructs the real-space generator of cyclic translations
for a discretised interval of size N = 2^n.

Definitions
-----------
Let T be the N×N cyclic shift operator acting on a real-space grid:

    T_{i,j} = δ_{i, j+1 mod N}

Then the generator G of discrete translations is defined as:

    T = exp(G)

so that

    G = log(T)

This object is used as the infinitesimal generator of time translations
on a discrete real-space grid and appears naturally in the Walsh formulation
of Floquet extended spaces with discrete times.

Output
------
The script saves a pickle file:

    ../../data/processed/discrete_generator.pkl

containing:
    - "N": number of discrete time steps
    - "G": complex-valued generator matrix (NxN)
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
import scipy as sp

# Move into repo/src directory to import walsh_floquet
p = Path.cwd()
os.chdir(p)
p = p.parent / "src"
os.chdir(p)

import sys
sys.path.append(str(p))
from walsh_floquet import *   # Import Walsh utilities


# --------------------------------------------------------------------------- #
# Constants and setup
# --------------------------------------------------------------------------- #

n = 6                   # number of bits → N = 2^n discrete time samples
N = 2 ** n              # grid size


# --------------------------------------------------------------------------- #
# Construct cyclic translation operator and its generator
# --------------------------------------------------------------------------- #

# Cyclic shift matrix T: shifts vector elements one step left with wrap-around
I = np.eye(N)
T = np.roll(I, shift=-1, axis=0)

# Generator G = log(T), computed using scipy.linalg.logm (via walsh_floquet import)
G = sp.logm(T)


# --------------------------------------------------------------------------- #
# Save generator
# --------------------------------------------------------------------------- #

filename       = "discrete_generator.pkl"
processed_path = os.path.join("..", "..", "data", "processed")
save_path      = os.path.join(processed_path, filename)

os.makedirs(processed_path, exist_ok=True)

to_save = {
    "N": N,    # Number of discrete time steps
    "G": G,    # Generator of cyclic translations
}

with open(save_path, "wb") as f:
    pickle.dump(to_save, f)

print("Saved data to:", save_path)

# --------------------------------------------------------------------------- #
# End of Script
# --------------------------------------------------------------------------- #





