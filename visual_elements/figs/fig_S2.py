# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
import pickle
import numpy as np
import matplotlib

matplotlib.use("Agg")  # render without display backend
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# ============================================================================ #
# LATEX / STYLE SETTINGS
# ============================================================================ #

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 9.54
plt.rcParams["axes.titlesize"] = 9.54
plt.rcParams["text.latex.preamble"] = (
    r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{mathptmx}"
)


# ============================================================================ #
# PATHS & FILENAMES
# ============================================================================ #

processed_path = os.path.join("..", "..", "data", "processed")
filename = "discrete_generator"   # readable name for saving figures


# ============================================================================ #
# LOAD PRECOMPUTED DATA
# ============================================================================ #

load_path = os.path.join(processed_path, "discrete_generator.pkl")
with open(load_path, "rb") as f:
    D = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

N = D["N"]  # Number of Walsh modes
G = D["G"]  # Time translation generator matrix in real space


# ============================================================================ #
# COLOURS & CUSTOM COLORMAP
# ============================================================================ #

ORANGE = "#E97C4A"
WHITE  = "#FFFFFF"
GREEN  = "#4B8B3B"

OWGcmap = mcolors.LinearSegmentedColormap.from_list(
    "CustomMap",
    [ORANGE, WHITE, GREEN],
)


# ============================================================================ #
# FIGURE GEOMETRY (APS 2-column, journal-friendly)
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6                 # total width in cm
fig_width    = fig_width_cm * cm   # in inches

aspect_ratio = 2.0 / 3.0
fig_height   = fig_width * aspect_ratio


# ============================================================================ #
# LAYOUT
# ============================================================================ #

fig, axs = plt.subplots(
    1,
    2,
    figsize=(fig_width, fig_height),
    constrained_layout=True,
)


# ============================================================================ #
# PLOTTING — REAL AND IMAGINARY PARTS OF G
# ============================================================================ #

# Common colour limits
vmin = -1
vmax = 1

# Real part
im0 = axs[0].imshow(
    np.real(G),
    vmin=vmin,
    vmax=vmax,
    cmap="viridis",
)
axs[0].set_title(r"\textbf{Re}($\hat{G}$)")

# Imaginary part
im1 = axs[1].imshow(
    np.imag(G),
    vmin=vmin,
    vmax=vmax,
    cmap="viridis",
)
axs[1].set_title(r"\textbf{Im}($\hat{G}$)")

# Axis ticks and limits
for ax in axs:
    k = 0.8
    ax.set_xlim(-k, N - 1 + k)
    ax.set_ylim(N - 1 + k, -k)
    ax.set_xticks([0, N - 1])
    ax.set_xticklabels([r"$0$", r"$N$"])
    ax.set_yticks([0, N - 1])
    ax.set_yticklabels([r"$0$", r"$N$"])


# ============================================================================ #
# COLORBAR
# ============================================================================ #

cbar = fig.colorbar(
    im1,
    ax=axs,
    orientation="vertical",
    fraction=0.022,
)
cbar.ax.tick_params(labelsize=9.54)
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels([r"$-N/T$", r"$0$", r"$N/T$"])


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches='tight')
