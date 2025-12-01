# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
import pickle
import numpy as np
import matplotlib

matplotlib.use("Agg")      # render without display backend
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm


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

processed_path = os.path.join("..", "..", "data", "processed")
filename = "spin_kick"

omega = 10.0


# ============================================================================ #
# LOAD PRECOMPUTED SPECTRUM DATA
# ============================================================================ #

with open(os.path.join(processed_path, "spin_kick_spectrum.pkl"), "rb") as f:
    D = pickle.load(f)

hxs_W        = D["hxs_W"]
hxs_F        = D["hxs_F"]
hxs_num      = D["hxs_num"]
evals_W_flat = D["evals_W_flat"]
evals_F_flat = D["evals_F_flat"]
numsols_flat = D["numsols_flat"]


# ============================================================================ #
# LOAD PRECOMPUTED ERROR DATA
# ============================================================================ #

with open(os.path.join(processed_path, "spin_kick_errors.pkl"), "rb") as f:
    D2 = pickle.load(f)

hxs           = D2["hxs"]
hzs           = D2["hzs"]
logErrorRatio = D2["logError"]


# ============================================================================ #
# DERIVED PARAMETERS
# ============================================================================ #

T = 2 * np.pi / omega
steps = 80


# ============================================================================ #
# COLOURS & CUSTOM COLORMAP
# ============================================================================ #

ORANGE = "#E97C4A"
WHITE  = "#FFFFFF"
GREEN  = "#4B8B3B"

OWGcmap = mcolors.LinearSegmentedColormap.from_list(
    "CustomMap",
    [ORANGE, WHITE, GREEN]
)


# ============================================================================ #
# FIGURE GEOMETRY (APS 2-column, journal-friendly)
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6                    # typical APS 2-column figure width
fig_width    = fig_width_cm * cm

# left panel : right panel width ratio
main_ratio = 0.7
side_ratio = 1.0

side_width = fig_width * side_ratio / (main_ratio + side_ratio)
fig_height = 0.78 * side_width        # tuned for aesthetics

fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)

# modern API replacing deprecated set_constrained_layout_pads
fig.get_layout_engine().set(
    h_pad=0.0,
    w_pad=0.05,
    hspace=0.0,
    wspace=0.0
)


# ============================================================================ #
# LAYOUT
# ============================================================================ #

gs = fig.add_gridspec(
    nrows=1,
    ncols=2,
    width_ratios=[main_ratio, side_ratio],
    wspace=0.001
)

ax_spec = fig.add_subplot(gs[0, 0])   # left: quasienergy spectrum
ax_err  = fig.add_subplot(gs[0, 1])   # right: Walsh vs Fourier error


# ============================================================================ #
# RIGHT PANEL — log10(ΔθF / ΔθW)
# ============================================================================ #

im = ax_err.imshow(
    logErrorRatio,
    origin="lower",
    extent=[hxs[0], hxs[-1], hzs[0] * T, hzs[-1] * T],
    aspect=0.5,
    cmap=OWGcmap,
    vmin=-1,
    vmax=1,
)

cbar = fig.colorbar(im, ax=ax_err, fraction=0.05, pad=0.0)

title = ax_err.set_title(r"$\log_{10}(\Delta \theta_F / \Delta \theta_W)$")
offset = mtransforms.ScaledTranslation(0.15, 0, ax_err.figure.dpi_scale_trans)
title.set_transform(title.get_transform() + offset)

ax_err.set_xlabel(r"$\mathrm{Kick\ Field}\ h_x$")
ax_err.set_ylabel(r"$\mathrm{Static\ Field}\ h_z T$")

ax_err.set_xticks([0, np.pi/4, np.pi/2],
                  [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"])
ax_err.set_yticks([0, np.pi/2, np.pi],
                  [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])


# ============================================================================ #
# LEFT PANEL — QUASIENERGY SPECTRUM
# ============================================================================ #

# Subsample exact solutions for clarity
idx  = np.arange(1200) % 41
mask = (idx < 2)

ax_spec.scatter(hxs_W, evals_W_flat * T, s=0.05, label="Walsh",   c=GREEN,  rasterized=True)
ax_spec.scatter(hxs_F, evals_F_flat * T, s=0.05, label="Fourier", c=ORANGE, rasterized=True)
ax_spec.scatter(hxs_num[mask], np.real(numsols_flat[mask]) * T,
                s=0.05, label="Exact", c="black", marker="x")

ax_spec.set_xlabel(r"Kick Field, $h_x$")
ax_spec.set_ylabel(r"Quasienergy Phase, $\theta_n$", labelpad=-1)

ax_spec.set_xticks([0, np.pi/4, np.pi/2],
                   [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"])
ax_spec.set_yticks([-np.pi, 0, np.pi],
                   [r"$-\pi$", r"$0$", r"$\pi$"])

ax_spec.set_xlim([0, np.pi/2])
ax_spec.set_ylim([-np.pi, np.pi])

ax_spec.legend(
    loc="lower left",
    markerscale=5,
    handlelength=0.5,
    handletextpad=0.1,
    bbox_to_anchor=(-0.05, 0.24),
    labelspacing=0.0,
    borderpad=0.2
)


# ============================================================================ #
# PANEL LABELS
# ============================================================================ #

fig.text(0.07, 0.95, r"(a)", fontsize=10)
fig.text(0.54, 0.95, r"(b)", fontsize=10)


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches="tight", pad_inches=0)





