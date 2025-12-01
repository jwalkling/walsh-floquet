# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
from pathlib import Path
import pickle

import numpy as np
import matplotlib

matplotlib.use("Agg")  # render without display backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import AnnotationBbox, OffsetImage


# ============================================================================ #
# LATEX / STYLE SETTINGS
# ============================================================================ #

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["text.latex.preamble"] = (
    r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{mathptmx}"
)


# ============================================================================ #
# PATHS & FILENAMES
# ============================================================================ #

base_path      = Path("..") / ".."
processed_path = base_path / "data" / "processed"

fig_filename = "alias_sqwave_res_polariton.pdf"


# ============================================================================ #
# LOAD PRECOMPUTED DATA
# ============================================================================ #

# Resonance effect on localisation of spin down for kick drive
load_path = processed_path / "spin_kick_loc_down.pkl"
with open(load_path, "rb") as f:
    D = pickle.load(f)

# Square-drive localisation / error plots
load_path2 = processed_path / "spin_square_loc_errors.pkl"
with open(load_path2, "rb") as f:
    D2 = pickle.load(f)

# Walsh polariton time evolution
load_path3 = processed_path / "walsh_polariton_t.pkl"
with open(load_path3, "rb") as f:
    D3 = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

# Resonance for spin-down localisation
hxs_r         = D["hxs"]
hzs_r         = D["hzs"]
Loc_Diff_down = D["Loc_Diff_down"]

# Error and localisation data (square wave)
hxs               = D2["hxs"]
hzs               = D2["hzs"]
Sq_PhotonPE_Ratio = D2["Loc_Ratio"]
Sq_Errors_Ratio   = D2["Errors_Ratio"]

# Mode and time data for polariton evolution
ts           = D3["ts"]
ut_high_real = D3["ut_high_real"]
ut_low_imag  = D3["ut_low_imag"]


# ============================================================================ #
# DERIVED PARAMETERS
# ============================================================================ #

omega = 10.0
T     = 2 * np.pi / omega
steps = 30

# Note: mode_labels depends on Wm_high in other scripts; kept untouched here
# to preserve original logic, even though Wm_high is not defined in this file.
# mode_labels = np.arange(1, len(Wm_high) + 1)


# ============================================================================ #
# COLOURS & CUSTOM COLORMAP
# ============================================================================ #

color_F = "#E97C4A"  # Fourier (not used directly, but kept for consistency)
color_W = "#4B8B3B"  # Walsh

ORANGE = "#E97C4A"
WHITE  = "#FFFFFF"
GREEN  = "#4B8B3B"

OWGcmap = mcolors.LinearSegmentedColormap.from_list(
    "CustomMap",
    [ORANGE, WHITE, GREEN],
)


# ============================================================================ #
# FIGURE GEOMETRY
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6 * 2
fig_width    = fig_width_cm * cm

aspect_ratio = 0.55
fig_height   = fig_width * aspect_ratio

fig = plt.figure(figsize=(fig_width, fig_height))

gs = GridSpec(
    2, 3,
    figure=fig,
    width_ratios=[0.3, 1.7, 1.7],
    wspace=0.35,
    hspace=0.5,  # small but nonzero to avoid visual overlap
)

ax_top_left  = fig.add_subplot(gs[0, 1])  # S_F / S_W
ax_top_right = fig.add_subplot(gs[0, 2])  # |Δθ_F| / |Δθ_W|
ax_bot_left  = fig.add_subplot(gs[1, 1])  # ΔS (resonance)
ax_bot_right = fig.add_subplot(gs[1, 2])  # time series


# ============================================================================ #
# TOP-LEFT PANEL — S_F / S_W RATIO
# ============================================================================ #

numvals = Sq_PhotonPE_Ratio.shape[0]
hxs_top = omega * np.linspace(0, 6, numvals)
hzs_top = hxs_top.copy()

extent_top = [
    hxs_top[0] * T,
    hxs_top[-1] * T,
    hzs_top[0] * T,
    hzs_top[-1] * T,
]

xtick_values = np.arange(0, 13 * np.pi, 6 * np.pi)
xtick_labels = [rf"${i}\pi$" if i > 0 else r"$0$" for i in range(0, 13, 6)]

corr = np.nan_to_num(Sq_PhotonPE_Ratio, nan=0.0)

im1 = ax_top_left.imshow(
    corr,
    origin="lower",
    cmap=OWGcmap,
    vmin=0,
    vmax=2,
    extent=extent_top,
    aspect="equal",
)

ax_top_left.set_xlabel(r"$h_x T$")
ax_top_left.set_ylabel(r"$h_z T$")

ax_top_left.set_xticks(xtick_values)
ax_top_left.set_xticklabels(xtick_labels)
ax_top_left.set_yticks(xtick_values)
ax_top_left.set_yticklabels(xtick_labels)

cbar1 = fig.colorbar(im1, ax=ax_top_left, fraction=0.046, pad=0.04)
cbar1.set_ticks([0, 1, 2])
t1 = cbar1.ax.set_title(r"$S_F/S_W$", pad=5, loc="left")
t1.set_x(-4)


# ============================================================================ #
# TOP-RIGHT PANEL — |Δθ_F| / |Δθ_W| RATIO
# ============================================================================ #

im2 = ax_top_right.imshow(
    Sq_Errors_Ratio,
    origin="lower",
    cmap=OWGcmap,
    vmin=0,
    vmax=2,
    extent=extent_top,
    aspect="equal",
)

ax_top_right.set_xlabel(r"$h_x T$")
ax_top_right.set_ylabel(r"$h_z T$")

ax_top_right.set_xticks(xtick_values)
ax_top_right.set_xticklabels(xtick_labels)
ax_top_right.set_yticks(xtick_values)
ax_top_right.set_yticklabels(xtick_labels)

cbar2 = fig.colorbar(im2, ax=ax_top_right, fraction=0.046, pad=0.04)
cbar2.set_ticks([0, 1, 2])
cbar2.set_ticklabels([r"$0$", r"$1$", r"$2$"])
t2 = cbar2.ax.set_title(
    r"$|\Delta\theta_F|/|\Delta\theta_W|$",
    pad=5,
    loc="left",
)
t2.set_x(-8)


# ============================================================================ #
# BOTTOM-LEFT PANEL — ΔS RESONANCE MAP
# ============================================================================ #

steps_loc = Loc_Diff_down.shape[0]
hzs_loc   = (omega / 2) * np.linspace(0, 1, steps_loc)
hxs_loc   = hzs_loc * T / 2

extent_loc = [
    hxs_loc[0],
    hxs_loc[-1],
    hzs_loc[0] * T,
    hzs_loc[-1] * T,
]

im3 = ax_bot_left.imshow(
    Loc_Diff_down,
    origin="lower",
    extent=extent_loc,
    aspect=0.5,
    cmap=OWGcmap,
)

cbar3 = fig.colorbar(im3, ax=ax_bot_left, fraction=0.046, pad=0.04)
cbar3.set_ticks([-1, 0, 1])
cbar3.set_ticklabels([r"$-1$", r"$0$", r"$1$"])
cbar3.ax.set_title(r"$\Delta S$", pad=5)

ax_bot_left.set_xlabel(r"$h_x$")
ax_bot_left.set_ylabel(r"$h_z T$")

ax_bot_left.set_xticks([0, np.pi / 4, np.pi / 2])
ax_bot_left.set_xticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$"])

ax_bot_left.set_yticks([0, np.pi / 2, np.pi])
ax_bot_left.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"])

ax_bot_left.text(
    0.15,
    0.95,
    "Resonance",
    transform=ax_bot_left.transAxes,
    ha="left",
    va="top",
    fontsize=9,
    fontweight="bold",
)


# ============================================================================ #
# BOTTOM-RIGHT PANEL — POLARITON TIME EVOLUTION
# ============================================================================ #

T_plot = ts[-1]

# This dummy call was present in the original code; kept for fidelity
ax_bot_right.plot(aspect=0.2)

ax_bot_right.plot(
    ts,
    ut_high_real,
    label=r"Re$\{ \langle \uparrow| u(t) \rangle \}$",
    color="orchid",
)
ax_bot_right.plot(
    ts,
    ut_low_imag,
    label=r"Im$\{ \langle \downarrow| u(t) \rangle \}$",
    linestyle="dashed",
    color="purple",
)

ax_bot_right.set_xticks([0, T_plot / 2, T_plot])
ax_bot_right.set_xticklabels([r"$0$", r"$T/2$", r"$T$"])

ax_bot_right.set_yticks([-1, 0, 1])
ax_bot_right.set_yticklabels([r"$-1$", r"$0$", r"$1$"])

ax_bot_right.legend(
    loc="lower right",
    bbox_to_anchor=(1.03, 0.08),
    borderpad=0.0,
    handletextpad=0.2,
    frameon=False,
    handlelength=1.0,
)

# Need a draw to get final positions
fig.canvas.draw()
pos = ax_bot_right.get_position()

# Shrink time-series panel vertically, keeping its top fixed
shrink_factor = 0.5
new_h = pos.height * shrink_factor
new_y = pos.y0 + (pos.height - new_h)

ax_bot_right.set_position([pos.x0, new_y, pos.width, new_h])


# ============================================================================ #
# EQUATION PANEL (TEXT REGION)
# ============================================================================ #

# Equation placed in figure coordinates (kept as in original)
fig.text(
    0.75,
    0.15,
    r"$|u_n(t)\rangle = \frac{1}{\sqrt{2}}"
    r"\begin{pmatrix} 1 \\ -i\, W_{N/2}(t) \end{pmatrix}$",
    ha="center",
    va="center",
    fontsize=10,
)


# ============================================================================ #
# PANEL LABELS & EMBEDDED IMAGE
# ============================================================================ #

fig.text(0.0,  0.93, r"(a)", fontsize=10, va="top", ha="left")
fig.text(0.28, 0.93, r"(b)", fontsize=10, va="top", ha="left")
fig.text(0.28, 0.47, r"(c)", fontsize=10, va="top", ha="left")
fig.text(0.58, 0.47, r"(d)", fontsize=10, va="top", ha="left")

arr = mpimg.imread("dft_vs_cft.png") #import the image for aliasing
imagebox = OffsetImage(arr, zoom=0.15)
ab = AnnotationBbox(
    imagebox,
    (0.16, 0.5),
    frameon=False,
    xycoords="figure fraction",
)
fig.add_artist(ab)


# ============================================================================ #
# SAVE
# ============================================================================ #

fig.savefig(fig_filename, dpi=600, bbox_inches="tight")
