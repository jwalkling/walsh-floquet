# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
import pickle
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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

processed_path = Path("..") / ".." / "data" / "processed"
data_filename  = "walsh_series.pkl"
filename       = "walsh_series"   # base name for saving figure

load_path = processed_path / data_filename


# ============================================================================ #
# LOAD PRECOMPUTED DATA
# ============================================================================ #

with open(load_path, "rb") as f:
    D = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

tvals     = D["tvals"]
W13       = D["W13"]
omegas    = D["omegas"]
errors_0  = D["errors_0"]
errors_01 = D["errors_01"]


# ============================================================================ #
# FIGURE GEOMETRY (APS 2-COLUMN WIDTH)
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6
fig_width    = fig_width_cm * cm  # inches

top_ratio    = 0.2
bottom_ratio = 0.7
total_ratio  = top_ratio + bottom_ratio

fig_height = fig_width * total_ratio


# ============================================================================ #
# LAYOUT
# ============================================================================ #

fig, (ax_top, ax_bottom) = plt.subplots(
    2,
    1,
    figsize=(fig_width, fig_height),
    gridspec_kw={"height_ratios": [top_ratio, bottom_ratio]},
    constrained_layout=True,
)


# ============================================================================ #
# PANEL (a): DRIVE WAVEFORM
# ============================================================================ #

ax_top.step(tvals, W13, where="post")

# Use actual time range from data
ax_top.set_xlim(tvals[0] - 1, tvals[-1] + 1)
ax_top.set_ylim(-1.2, 1.2)

# ticks at t = 0 and t = T
ax_top.set_xticks([tvals[0], tvals[-1]])
ax_top.set_xticklabels([r"0", r"$T$"])

ax_top.set_yticks([-1, 0, 1])
ax_top.set_yticklabels([r"$-1$", r"$0$", r"$1$"])

ax_top.tick_params(axis="x", direction="out")
ax_top.tick_params(axis="y", direction="out")

ax_top.set_xlabel("")  # no centered xlabel
ax_top.set_ylabel("")  # no y-label

# bottom-right x-label inside axes
ax_top.text(
    1,
    -0.05,
    r"$t$",
    transform=ax_top.transAxes,
    ha="right",
    va="top",
    fontsize=10,
)

# panel label (a)
ax_top.text(
    -0.12,
    1.05,
    r"(a)",
    transform=ax_top.transAxes,
    fontsize=10,
    va="top",
    ha="left",
)


# ============================================================================ #
# PANEL (b): ERROR VS FREQUENCY SCALING
# ============================================================================ #

x  = np.log10(omegas)
y1 = np.log10(errors_0)
y2 = np.log10(errors_01)

# reference lines ∼ ω^{-1} and ∼ ω^{-2}
ax_bottom.plot(
    x,
    np.log10(omegas**(-1)) + 1.15,
    linestyle="--",
    color="lightblue",
    linewidth=1,
)
ax_bottom.plot(
    x,
    np.log10(omegas**(-2)) + 1.1,
    linestyle="--",
    color="orange",
    linewidth=1,
)

# numerical data
ax_bottom.plot(
    x,
    y1,
    label=r"$H_F - H_\mathrm{eff}^{(0)}$",
    color="#1f77b4",
)
ax_bottom.plot(
    x,
    y2,
    label=r"$H_F - H_\mathrm{eff}^{(0+1)}$",
    color="#ff7f0e",
)

ax_bottom.set_xlabel(r"$\log_{10}(\omega)$")
ax_bottom.set_ylabel(r"$\log_{10}(|\Delta \varepsilon|)$")
ax_bottom.grid(True)

# power-law annotations
ax_bottom.text(
    0.6,
    0.72,
    r"$\sim \omega^{-1}$",
    transform=ax_bottom.transAxes,
    fontsize=10,
)
ax_bottom.text(
    0.6,
    0.12,
    r"$\sim \omega^{-2}$",
    transform=ax_bottom.transAxes,
    fontsize=10,
)

# legend
ax_bottom.legend(frameon=True, loc="lower left")

# panel label (b)
ax_bottom.text(
    -0.12,
    1.05,
    r"(b)",
    transform=ax_bottom.transAxes,
    fontsize=10,
    va="top",
    ha="left",
)


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches="tight")
