# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
import pickle
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap


# ============================================================================ #
# LATEX / STYLE SETTINGS
# ============================================================================ #

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10 * 0.94
plt.rcParams["axes.titlesize"] = 10 * 0.94
plt.rcParams["text.latex.preamble"] = (
    r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{mathptmx}"
)


# ============================================================================ #
# PATHS & FILENAMES
# ============================================================================ #

processed_path = os.path.join("..", "..", "data", "processed")

filename = "H_extended_kick"   # readable name for saving figures


# ============================================================================ #
# LOAD PRECOMPUTED DATA
# ============================================================================ #

load_path = os.path.join(processed_path, "H_extended_kick.pkl")
with open(load_path, "rb") as f:
    D = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

drive_W = D["drive_W"]
drive_F = D["drive_F"]
deriv_W = D["deriv_W"]
deriv_F = D["deriv_F"]

N    = drive_W.shape[0]  # number of modes
half = N // 2


# ============================================================================ #
# COLOURS & CUSTOM COLORMAPS
# ============================================================================ #

# Red-only colormap: white → light pink → strong red
red_cmap = LinearSegmentedColormap.from_list(
    "red_only",
    [
        (0.0, "#ffffff"),  # white
        (0.3, "#ffe5e5"),  # very light pink
        (1.0, "#cc0000"),  # strong red
    ],
)

# Blue-only colormap: white → light blue → strong blue
blue_cmap = LinearSegmentedColormap.from_list(
    "blue_only",
    [
        (0.0, "#ffffff"),  # white
        (0.3, "#e5f0ff"),  # very light blue
        (1.0, "#0033aa"),  # strong blue
    ],
)


# ============================================================================ #
# FIGURE GEOMETRY (APS 2-column, journal-friendly)
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6
fig_width    = fig_width_cm * cm

width_in  = fig_width
height_in = width_in * 2.0 / 3.0


# ============================================================================ #
# LAYOUT
# ============================================================================ #

fig, axs = plt.subplots(
    1,
    2,
    figsize=(width_in, height_in),
    constrained_layout=True,
)


# ============================================================================ #
# SHARED SCALING
# ============================================================================ #

# Use magnitudes; set common scales so Fourier and Walsh are comparable
drive_max = max(np.abs(drive_F).max(), np.abs(drive_W).max())
deriv_max = max(np.abs(deriv_F).max(), np.abs(deriv_W).max())

vmin_drive, vmax_drive = 0.0, drive_max
vmin_deriv, vmax_deriv = 0.0, deriv_max


# ============================================================================ #
# PANEL (a): Q̄ IN FOURIER BASIS
# ============================================================================ #

ax = axs[0]

# Drive (H) in red
ax.imshow(
    np.abs(drive_F),
    cmap=red_cmap,
    vmin=vmin_drive,
    vmax=vmax_drive,
    interpolation="nearest",
)

# Derivative (-i ∂_t) in blue (scaled by 10)
ax.imshow(
    10 * np.abs(deriv_F),
    cmap=blue_cmap,
    vmin=vmin_deriv,
    vmax=vmax_deriv,
    interpolation="nearest",
    alpha=0.9,
)

ax.set_title(r"$\bar{Q}$ in Fourier basis", pad=5)

ax.set_xticks([0, half, N - 1])
ax.set_xticklabels([r"$0$", r"$N/2$", r"$N$"])
ax.set_yticks([0, half, N - 1])
ax.set_yticklabels([r"$0$", r"$N/2$", r"$N$"])

ax.text(0.15 * N, 0.85 * N, r"$H$",            color="red")
ax.text(0.80 * N, 0.15 * N, r"$H$",            color="red")
ax.text(0.60 * N, 0.55 * N, r"$-i\partial_t$", color="blue")
ax.text(-0.25 * N, -0.10 * N, r"$(a)$")


# ============================================================================ #
# PANEL (b): Q̄ IN WALSH BASIS
# ============================================================================ #

ax = axs[1]

# Drive (H) in red
ax.imshow(
    np.abs(drive_W),
    cmap=red_cmap,
    vmin=vmin_drive,
    vmax=vmax_drive,
    interpolation="nearest",
)

# Derivative (-i Ĝ) in blue (scaled by 10)
ax.imshow(
    10 * np.abs(deriv_W),
    cmap=blue_cmap,
    vmin=vmin_deriv,
    vmax=vmax_deriv,
    interpolation="nearest",
    alpha=0.9,
)

ax.set_title(r"$\bar{Q}$ in Walsh basis", pad=5)

ax.set_xticks([0, half, N - 1])
ax.set_xticklabels([r"$0$", r"$N/2$", r"$N$"])
ax.set_yticks([0, half, N - 1])
ax.set_yticklabels([r"$0$", r"$N/2$", r"$N$"])

ax.text(0.15 * N, 0.85 * N, r"$H$",       color="red")
ax.text(0.80 * N, 0.15 * N, r"$H$",       color="red")
ax.text(0.18 * N, 0.12 * N, r"$-i\hat{G}$", color="blue")
ax.text(-0.25 * N, -0.10 * N, r"$(b)$")

for ax in axs:
    ax.tick_params(direction="out", length=2)


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches="tight")
