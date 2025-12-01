# ============================================================================ #
# IMPORTS & BACKEND
# ============================================================================ #

import os
import pickle
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms


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
filename = "error_scaling"   # readable name for saving figures


# ============================================================================ #
# LOAD PRECOMPUTED ERROR-SCALING DATA
# ============================================================================ #

load_path = os.path.join(processed_path, "square_kick_error_scaling.pkl")
with open(load_path, "rb") as f:
    D = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

# Range of modes (as stored, then overridden as in original script)
nvals   = D["nvals"]
nvalsmb = D["nvalsmb"]

# Explicit ranges used in the plotting (kept exactly as in original code)
nvals   = np.arange(2, 11, 1)   # n = 2,...,10
nvalsmb = np.arange(3, 10, 1)   # many-body truncations for L = 3

# Errors
mb_kick_F = D["mbkick_F"]
sp_kick_F = D["spkick_F"]
mb_kick_W = D["mbkick_W"]
sp_kick_W = D["spkick_W"]
mb_sq_F   = D["mbsq_F"]
sp_sq_F   = D["spsq_F"]
mb_sq_W   = D["mbsq_W"]
sp_sq_W   = D["spsq_W"]


# ============================================================================ #
# COLOURS
# ============================================================================ #

color_W = "#4B8B3B"   # Walsh (green)
color_F = "#E97C4A"   # Fourier (orange)

plot_style = dict(marker="o", linestyle="-", markersize=3)


# ============================================================================ #
# FIGURE GEOMETRY (APS 2-column, JOURNAL-FRIENDLY)
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6
fig_width    = fig_width_cm * cm

# left column (single-particle) : right column (many-body)
main_ratio = 1.0
side_ratio = 1.0

# side panels square → figure height equals side panel width
side_width = fig_width * side_ratio / (main_ratio + side_ratio)
fig_height = side_width

fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
fig.set_constrained_layout_pads(
    h_pad=0.0,
    w_pad=0.05,
    hspace=0.0,
    wspace=0.0,
)


# ============================================================================ #
# LAYOUT
# ============================================================================ #

gs = fig.add_gridspec(
    nrows=2,
    ncols=2,
    width_ratios=[main_ratio, side_ratio],
    height_ratios=[1.0, 1.0],
    wspace=0.001,
)

ax_1 = fig.add_subplot(gs[0, 0])  # single-particle, square drive
ax_2 = fig.add_subplot(gs[1, 0])  # single-particle, kick drive
ax_4 = fig.add_subplot(gs[0, 1])  # many-body, square drive
ax_3 = fig.add_subplot(gs[1, 1])  # many-body, kick drive


# ============================================================================ #
# DERIVED PARAMETERS
# ============================================================================ #

x_sp = 2 ** nvals
x_mb = 2 ** nvalsmb


# ============================================================================ #
# SINGLE-PARTICLE PANELS (LEFT COLUMN)
# ============================================================================ #

# (a) Single-particle square drive
ax_1.loglog(x_sp, sp_sq_F, label="Fourier", color=color_F, **plot_style)
ax_1.loglog(x_sp, sp_sq_W, label="Walsh",   color=color_W, **plot_style)
ax_1.set_ylabel(r"Error, $\Delta \theta$")
ax_1.grid(True)

ax_1.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax_1.set_yticks([10**-2, 10**-4, 10**-6, 10**-8])

# (c) Single-particle kick drive
ax_2.loglog(x_sp, sp_kick_F, label="Fourier", color=color_F, **plot_style)
ax_2.loglog(x_sp, sp_kick_W, label="Walsh",   color=color_W, **plot_style)
ax_2.set_xlabel(r"Number of Modes, $N$")
ax_2.grid(True)
ax_2.set_yticks([10**-1, 10**-2, 10**-3])


# ============================================================================ #
# MANY-BODY PANELS (RIGHT COLUMN)
# ============================================================================ #

# (d) Many-body kick drive
ax_3.loglog(x_mb, mb_kick_F, label="Fourier", color=color_F, **plot_style)
ax_3.loglog(x_mb, mb_kick_W, label="Walsh",   color=color_W, **plot_style)
ax_3.set_xlabel(r"Number of Modes, $N$")
ax_3.grid(True)
ax_3.set_yticks([10**-1, 10**-2, 10**-3])

# (b) Many-body square drive
ax_4.loglog(x_mb, mb_sq_F, label="Fourier", color=color_F, **plot_style)
ax_4.loglog(x_mb, mb_sq_W, label="Walsh",   color=color_W, **plot_style)
ax_4.set_ylabel(r"Error, $\langle \Delta \theta \rangle$")
ax_4.grid(True)

ax_4.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax_4.set_yticks([10**-2, 10**-4, 10**-6])


# ============================================================================ #
# AXIS LABEL NUDGES
# ============================================================================ #

dx_in = 0.03
dy_in = -0.3

for ax in (ax_1, ax_4):
    label = ax.yaxis.label
    label.set_transform(
        label.get_transform()
        + mtransforms.ScaledTranslation(dx_in, dy_in, ax.figure.dpi_scale_trans)
    )


# ============================================================================ #
# IN-PANEL TEXT LABELS & TITLES
# ============================================================================ #

for x, y, text in [
    (0.00, 0.95, r"(a)"),
    (0.51, 0.95, r"(b)"),
    (0.20, 0.62, r"Square"),
    (0.72, 0.62, r"Square"),
    (0.20, 0.25, r"Kick"),
    (0.72, 0.25, r"Kick"),
]:
    fig.text(x, y, text, fontsize=10)

ax_1.set_title("Single-particle")
ax_4.set_title("Many-body")


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches="tight")
