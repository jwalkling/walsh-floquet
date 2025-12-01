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
from matplotlib.colors import SymLogNorm  # kept for consistency if needed later


# ============================================================================ #
# LATEX / STYLE SETTINGS
# ============================================================================ #

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = 10
plt.rcParams["figure.titlesize"] = 10
plt.rcParams["text.latex.preamble"] = (
    r"\usepackage{amsmath} \usepackage{amssymb} \usepackage{mathptmx}"
)


# ============================================================================ #
# PATHS & FILENAMES
# ============================================================================ #

processed_path = os.path.join("..", "..", "data", "processed")
filename = "mb_spin_kick"     # readable name for saving figures

omega = 10.0


# ============================================================================ #
# LOAD PRECOMPUTED SPECTRUM DATA
# ============================================================================ #

with open(os.path.join(processed_path, "3_spin_kick_spectrum.pkl"), "rb") as f:
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

with open(os.path.join(processed_path, "6_spin_kick_errors.pkl"), "rb") as f:
    D2 = pickle.load(f)

hxs         = D2["hxs"]
hzs         = D2["hzs"]
logErrorJ1  = D2["logErrorJ1"]
logErrorJ3  = D2["logErrorJ3"]
Errors_W_J1 = D2["Errors_W_J1"]
Errors_W_J3 = D2["Errors_W_J3"]


# ============================================================================ #
# DERIVED PARAMETERS
# ============================================================================ #

T = 2 * np.pi / omega
crit = 0.01   # critical value which defines the 'bad' region (masked)


# ============================================================================ #
# COLOURS & CUSTOM COLORMAP
# ============================================================================ #

ORANGE = "#E97C4A"
WHITE  = "#FFFFFF"
GREEN  = "#4B8B3B"

OWGcmap = mcolors.LinearSegmentedColormap.from_list(
    "OWG",
    [ORANGE, WHITE, GREEN]
)
OWGcmap.set_bad(color="black")


# ============================================================================ #
# FIGURE GEOMETRY (APS 2-column, journal-friendly)
# ============================================================================ #

cm = 1 / 2.54
fig_width_cm = 8.6                    # typical APS 2-column figure width
fig_width    = fig_width_cm * cm      # in inches

# main (left) panel : side (right) column width ratio
main_ratio = 2.0
side_ratio = 1.0

# side column width (inches); height is two stacked squares → 2 * side_width
side_width = fig_width * side_ratio / (main_ratio + side_ratio)
fig_height = 2.0 * side_width

fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)

# modern API replacing set_constrained_layout_pads
fig.get_layout_engine().set(
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

ax_spec   = fig.add_subplot(gs[:, 0])   # left main (spans both rows)
ax_top    = fig.add_subplot(gs[0, 1])   # top right square
ax_bottom = fig.add_subplot(gs[1, 1])   # bottom right square

# tighten the two square panels slightly inside their slots
for ax in (ax_top, ax_bottom):
    pos = ax.get_position()
    new_w = pos.width * 0.90
    new_h = pos.height * 0.90
    new_x0 = pos.x0 + (pos.width - new_w) / 2 + 0.13
    new_y0 = pos.y0 + (pos.height - new_h) / 2
    ax.set_position([new_x0, new_y0, new_w, new_h])

# slight vertical nudge for main axis (keeps labels from kissing page edge)
pos = ax_spec.get_position()
ax_spec.set_position([pos.x0, pos.y0 + 0.015, pos.width, pos.height])


# ============================================================================ #
# LEFT PANEL — QUASIENERGY SPECTRUM
# ============================================================================ #

ax_spec.scatter(
    hxs_W, evals_W_flat * T,
    s=0.05,
    label="Walsh",
    c=GREEN,
    rasterized=True,
)
ax_spec.scatter(
    hxs_F, evals_F_flat * T,
    s=0.05,
    label="Fourier",
    c=ORANGE,
    rasterized=True,
)
ax_spec.scatter(
    hxs_num, np.real(numsols_flat) * T,
    s=0.005,
    label="Exact",
    c="black",
    alpha=1.0,
    rasterized=True,
)

ax_spec.set_xlabel(r"Kick Field, $h_x$")
ax_spec.set_ylabel(r"Quasienergy Phase, $\theta_n$", labelpad=-3)

ax_spec.set_xticks(
    [0, np.pi / 4, np.pi / 2],
    [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"],
)
ax_spec.set_yticks(
    [-np.pi, 0, np.pi],
    [r"$-\pi$", r"$0$", r"$\pi$"],
)

ax_spec.set_xlim([0, np.pi / 2])
ax_spec.set_ylim([-np.pi, np.pi])

ax_spec.legend(
    loc="lower right",
    markerscale=10,
    handlelength=0.5,
    handletextpad=0.1,
    bbox_to_anchor=(1.03, -0.03),
)


# ============================================================================ #
# RIGHT PANELS — LOG10 ERROR RATIO FOR J=1, J=3
# ============================================================================ #

# Mask the large-error region: Δθ_W / 2π > crit
mlogErrorJ1 = np.ma.array(
    logErrorJ1,
    mask=(Errors_W_J1 / (2 * np.pi) > crit),
)
mlogErrorJ3 = np.ma.array(
    logErrorJ3,
    mask=(Errors_W_J3 / (2 * np.pi) > crit),
)

# Top panel: J = 1
im_top = ax_top.imshow(
    mlogErrorJ1,
    origin="lower",
    extent=[hxs[0], hxs[-1], hzs[0] * T, hzs[-1] * T],
    cmap=OWGcmap,
    aspect="auto",
    vmin=-1,
    vmax=1,
)

# Bottom panel: J = 3
im_bot = ax_bottom.imshow(
    mlogErrorJ3,
    origin="lower",
    extent=[hxs[0], hxs[-1], hzs[0] * T, hzs[-1] * T],
    cmap=OWGcmap,
    aspect="auto",
    vmin=-1,
    vmax=1,
)


# ============================================================================ #
# COLORBAR (SHARED, HORIZONTAL, TOP OF RIGHT PANELS)
# ============================================================================ #

cbar = fig.colorbar(
    im_top,
    ax=[ax_top, ax_bottom],
    orientation="horizontal",
    location="top",
    fraction=0.10,
    pad=0.03,
    label=r"$\log_{10}(\Delta \theta_F / \Delta \theta_W)$",
)
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels([r"$-1$", r"$0$", r"$1$"])
cbar.ax.xaxis.set_ticks_position("top")
cbar.ax.xaxis.set_label_position("top")

# Pull colorbar label slightly left
dx_in = -0.088
label = cbar.ax.xaxis.get_label()
label.set_transform(
    label.get_transform()
    + mtransforms.ScaledTranslation(dx_in, 0, cbar.ax.figure.dpi_scale_trans)
)

# Tighten tick labels to bar
cbar.ax.tick_params(axis="x", pad=0)


# ============================================================================ #
# AXIS LABELS & TICKS — RIGHT PANELS
# ============================================================================ #

ax_top.set_ylabel(r"Static Field, $h_z T$")
# Nudge y-label slightly
ax_top.yaxis.label.set_transform(
    ax_top.yaxis.label.get_transform()
    + mtransforms.ScaledTranslation(0.05, -0.4, ax_top.figure.dpi_scale_trans)
)

ax_top.set_xticks([0, np.pi / 4, np.pi / 2])
ax_top.set_xticklabels([])
ax_top.set_yticks(
    [0, np.pi, 2 * np.pi],
    [r"$0$", r"$\pi$", r"$2\pi$"],
)

ax_bottom.set_xlabel(r"Kick Field, $h_x$")
ax_bottom.set_xticks(
    [0, np.pi / 4, np.pi / 2],
    [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"],
)
ax_bottom.set_yticks(
    [0, np.pi, 2 * np.pi],
    [r"$0$", r"$\pi$", r"$2\pi$"],
)


# ============================================================================ #
# IN-PANEL TEXT LABELS (J=1 / J=3) & ERROR REGION LABELS
# ============================================================================ #

ax_top.text(
    0.05,
    0.90,
    r"$J=1$",
    transform=ax_top.transAxes,
    fontsize=10,
    va="top",
    ha="left",
)

ax_bottom.text(
    0.05,
    0.90,
    r"$J=3$",
    transform=ax_bottom.transAxes,
    fontsize=10,
    va="top",
    ha="left",
)

for ax in (ax_top, ax_bottom):
    ax.text(
        0.60,
        0.50,
        r"$\frac{\Delta \theta_W}{2\pi} > 0.01$",
        transform=ax.transAxes,
        rotation=70,
        va="center",
        ha="left",
        color="white",
        fontsize=5,
    )


# ============================================================================ #
# PANEL LABELS
# ============================================================================ #

fig.text(0.01, 0.92, r"(a)", fontsize=10)
fig.text(0.66, 0.92, r"(b)", fontsize=10)


# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches="tight", pad_inches=0)
