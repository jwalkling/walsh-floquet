# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
import pickle
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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
data_filename  = "strong_kick_error_scaling.pkl"
filename       = "error_scaling_strong_kick"   # base name for saving figure


# ============================================================================ #
# LOAD PRECOMPUTED ERROR-SCALING DATA
# ============================================================================ #

load_path = os.path.join(processed_path, data_filename)
with open(load_path, "rb") as f:
    D = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

nvals   = D["nvals"]
nvalsmb = D["nvalsmb"]

spstrongkick_F = D["spstrongkick_F"]
spstrongkick_W = D["spstrongkick_W"]

mbstrongkick_F = D["mbstrongkick_F"]
mbstrongkick_W = D["mbstrongkick_W"]


# ============================================================================ #
# DERIVED PARAMETERS
# ============================================================================ #

# Period corresponding to omega = 10.0
T = 2 * np.pi / 10.0


# ============================================================================ #
# COLOURS, MARKERS & LABELS
# ============================================================================ #

FC = "coral"        # Fourier colour
WC = "forestgreen"  # Walsh colour

marker_handles = [
    Line2D(
        [0], [0],
        marker="o",
        color=FC,
        linestyle="None",
        markersize=5,
        label="Fourier",
    ),
    Line2D(
        [0], [0],
        marker="o",
        color=WC,
        linestyle="None",
        markersize=5,
        label="Walsh",
    ),
]

plot_style = dict(marker="o", linestyle="-", markersize=3)

title_sp = r"Single-particle, $h_x=\pi/2$"
title_mb = r"Many-body, $h_x=1.1\pi/4$"

legend_kwargs = dict(
    handlelength=0.1,
    borderpad=0.3,
    handletextpad=0.3,
)


# ============================================================================ #
# FIGURE GEOMETRY (TWO SIDE-BY-SIDE SQUARE PANELS)
# ============================================================================ #

cm = 1 / 2.54
pane_width_cm = 4.3
gap_cm        = 0.2
height_scale  = 1.00

fig_width  = (2 * pane_width_cm + gap_cm) * cm
fig_height = (pane_width_cm * height_scale) * cm


# ============================================================================ #
# LAYOUT
# ============================================================================ #

fig, (ax_sp, ax_mb) = plt.subplots(
    1,
    2,
    figsize=(fig_width, fig_height),
    constrained_layout=True,
)

# Keep each axes approximately square
ax_sp.set_box_aspect(1)
ax_mb.set_box_aspect(1)


# ============================================================================ #
# SINGLE-PARTICLE PANEL
# ============================================================================ #

ax_sp.loglog(
    2**nvals,
    T * spstrongkick_F,
    label="Fourier",
    color=FC,
    **plot_style,
)
ax_sp.loglog(
    2**nvals,
    T * spstrongkick_W,
    label="Walsh",
    color=WC,
    **plot_style,
)

ax_sp.set_xlabel(r"Number of Modes, $N$")
ax_sp.set_ylabel(r"Error, $\Delta \theta$")
ax_sp.set_yticks([10**-1, 10**-2, 10**-3])
ax_sp.grid(True)
ax_sp.set_title(title_sp)
ax_sp.legend(**legend_kwargs)


# ============================================================================ #
# MANY-BODY PANEL
# ============================================================================ #

ax_mb.loglog(
    2**nvalsmb,
    T * mbstrongkick_F,
    label="Fourier",
    color=FC,
    **plot_style,
)
ax_mb.loglog(
    2**nvalsmb,
    T * mbstrongkick_W,
    label="Walsh",
    color=WC,
    **plot_style,
)

ax_mb.set_xlabel(r"Number of Modes, $N$")
ax_mb.set_ylabel(r"Error, $\langle \Delta \theta\rangle$")
ax_mb.set_yticks([10**-1, 10**-2, 10**-3])
ax_mb.grid(True)
ax_mb.set_title(title_mb)
ax_mb.legend(**legend_kwargs)


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches="tight", pad_inches=0)
