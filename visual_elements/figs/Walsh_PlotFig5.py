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
from matplotlib.colors import SymLogNorm  # kept for potential future use
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


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
filename = "localisation_vs_error"

omega = 10.0


# ============================================================================ #
# LOAD PRECOMPUTED DATA
# ============================================================================ #

# Time-evolution modes
with open(os.path.join(processed_path, "time_evolve_modes.pkl"), "rb") as f:
    D = pickle.load(f)

# Error and localisation data
with open(os.path.join(processed_path, "spin_kick_loc_errors.pkl"), "rb") as f:
    D2 = pickle.load(f)


# ============================================================================ #
# UNPACK DATA
# ============================================================================ #

# Mode and time data
ts          = D["ts"]
ut_high_real = D["ut_high_real"]
ut_low_real  = D["ut_low_real"]
Fm_high      = D["Fm_high"]
Fm_low       = D["Fm_low"]
Wm_high      = D["Wm_high"]
Wm_low       = D["Wm_low"]

# Error and localisation data
hxs        = D2["hxs"]
hzs        = D2["hzs"]
Loc_Diff   = D2["Loc_Diff"]
Errors_Diff = D2["Errors_Diff"]


# ============================================================================ #
# DERIVED PARAMETERS
# ============================================================================ #

T = 2 * np.pi / omega
steps = 30  # kept for possible future use
mode_labels = np.arange(1, len(Wm_high) + 1)  # mode indices starting from 1


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

# --- Figure sizing in cm (journal-friendly) ----------------------------------
cm = 1 / 2.54
fig_width_cm = 8.6                      # total width in cm
fig_width = fig_width_cm * cm           # in inches

# Relative widths (left main : right side column)
main_ratio = 2.5
side_ratio = 1.0

# Side column width (right) in inches
side_width = fig_width * side_ratio / (main_ratio + side_ratio)

# Thin top row for an invisible header axis
top_ratio = 0.3

# Ensure right-hand panels are square:
#   fig_height * (1 / (top_ratio + 2)) = side_width
#   → fig_height = side_width * (top_ratio + 2)
fig_height = side_width * (top_ratio + 2)

fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
fig.set_constrained_layout_pads(
    h_pad=0.0,
    w_pad=0.02,
    hspace=0.0,
    wspace=0.0,
)

# --- Layout ------------------------------------------------------------------
# 3 rows × 2 columns:
#   row 0: invisible header spanning both columns
#   row 1: main left + top-right
#   row 2: main left (continued) + bottom-right
gs = fig.add_gridspec(
    3,
    2,
    width_ratios=[main_ratio, side_ratio],
    height_ratios=[top_ratio, 1, 1],
    wspace=0.001,
)

ax_header = fig.add_subplot(gs[0, :])
ax_header.axis("off")

# Main plot (spans rows 1 and 2, left)
ax_main = fig.add_subplot(gs[1:, 0])

# Two side panels (right column)
ax_top = fig.add_subplot(gs[1, 1])
ax_bottom = fig.add_subplot(gs[2, 1])


# =============================================================================
# 1) Main panel: log-log spectra + legends + time-evolution inset
# =============================================================================

def safe_log10(arr):
    """Return log10(arr[mask]) and the mask, ignoring nonpositive entries."""
    arr = np.asarray(arr, dtype=float)
    mask = arr > 0
    return np.log10(arr[mask]), mask


# --- Spectral power series ---------------------------------------------------
P_W1 = np.sort(np.abs(Wm_high)**2)[::-1]
P_F1 = np.sort(np.abs(Fm_high)**2)[::-1]
P_W2 = np.sort(np.abs(Wm_low)**2)[::-1]
P_F2 = np.sort(np.abs(Fm_low)**2)[::-1]

x_all_log = np.log10(mode_labels)

y_W1_log, m_W1 = safe_log10(P_W1)
y_F1_log, m_F1 = safe_log10(P_F1)
y_W2_log, m_W2 = safe_log10(P_W2)
y_F2_log, m_F2 = safe_log10(P_F2)

# --- Plot spectra ------------------------------------------------------------
ax_main.plot(
    x_all_log[m_W1],
    y_W1_log,
    label=r'Walsh, $\omega = 10 h_z$',
    color='forestgreen',
)
ax_main.plot(
    x_all_log[m_F1],
    y_F1_log,
    label=r'Fourier, $\omega = 10 h_z$',
    color='coral',
)
ax_main.plot(
    x_all_log[m_W2],
    y_W2_log,
    label=r'Walsh, $\omega = 0.4 h_z$',
    color='darkolivegreen',
    linestyle='dashed',
)
ax_main.plot(
    x_all_log[m_F2],
    y_F2_log,
    label=r'Fourier, $\omega = 0.4 h_z$',
    color='peru',
    linestyle='dashed',
)

ax_main.set_xlabel(r'Log mode label, $\log_{10}(m)$')
ax_main.set_ylabel(r'Log mode power, $\log_{10}(|\tilde{u}_m|^2)$')

ax_main.set_xlim([0, 4])
ax_main.set_xticks([0, 1, 2, 3, 4])
ax_main.set_xticklabels([r"$0$", r"$1$", r"$2$", r"$3$", r"$4$"])

ax_main.set_ylim([-5, 0])
ax_main.set_yticks([-5, -4, -3, -2, -1, 0])
ax_main.set_yticklabels([r"$-5$", r"$-4$", r"$-3$", r"$-2$", r"$-1$", r"$0$"])

# --- Legends -----------------------------------------------------------------
style_legend = [
    Line2D([0], [0], color='black', linestyle='--',
           label=r'$\frac{\omega}{h_z} = 0.4$'),
    Line2D([0], [0], color='black', linestyle='-',
           label=r'$\frac{\omega}{h_z} = 10$'),
]

color_legend = [
    Line2D([0], [0], color='forestgreen', linestyle='-', label='Walsh'),
    Line2D([0], [0], color='coral', linestyle='-', label='Fourier'),
]

legend1 = ax_main.legend(
    handles=style_legend,
    loc='upper right',
    bbox_to_anchor=(1.05, 1.05),
    handletextpad=0,
    handlelength=2,
    handler_map={Line2D: HandlerLine2D(numpoints=2)},
    framealpha=1.0,
)

legend2 = ax_main.legend(
    handles=color_legend,
    loc='upper right',
    handlelength=0.5,
    handletextpad=0.1,
    bbox_to_anchor=(0.53, 1.05),
    framealpha=1.0,
)

ax_main.add_artist(legend1)
ax_main.add_artist(legend2)

# --- Inset: time evolution ---------------------------------------------------
ax_inset = inset_axes(
    ax_main,
    width="54%",
    height="36%",
    loc='lower right',
    bbox_to_anchor=(0.02, 0.13, 1, 1),
    bbox_transform=ax_main.transAxes,
)

ax_inset.plot(ts, ut_high_real, label=r'$\omega = 10 h_z$', color='orchid')
ax_inset.plot(
    ts,
    ut_low_real,
    label=r'$\omega = 0.4 h_z$',
    linestyle='dashed',
    color='purple',
)

ax_inset.set_title(r'$\mathrm{Re}\{\langle \uparrow | u(t) \rangle\}$')

ax_inset.spines['bottom'].set_position(('data', 0))

for spine in ax_inset.spines.values():
    spine.set_linewidth(1)
    spine.set_color("black")

ax_inset.set_xticks([0, T])
ax_inset.set_xticklabels([r"$0$", r"$T$"])
ax_inset.set_yticks([])

ax_inset.annotate(
    "",
    xy=(0.65 * T, -0.40),
    xytext=(0.2 * T, -0.40),
    arrowprops=dict(arrowstyle="->", lw=1, color="black"),
)
ax_inset.annotate(
    "",
    xy=(-0.045, -1.15),
    xytext=(1.08 * T, -1.15),
    arrowprops=dict(arrowstyle="-", lw=1, color="black"),
    annotation_clip=False,
)
ax_inset.text(0.4 * T, -0.8, r"Time, $t$", ha="center", va="center")

# =============================================================================
# 2) Top-right panel: 10^2 <Δθ>/π error heat map + colorbar
# =============================================================================

omega = 10
T = 2 * np.pi / omega

numvals = 40
hxs = np.linspace(0, np.pi / 2, numvals)
hzs = omega * np.linspace(0, 0.5, numvals)

error = ax_top.imshow(
    100 * Errors_Diff / np.pi,
    origin='lower',
    extent=[hxs[0], hxs[-1], hzs[0] * T, hzs[-1] * T],
    aspect=0.5,
    cmap=OWGcmap,
    vmin=-1,
    vmax=1,
)

ax_top.set_xticks([])
ax_top.set_xticklabels([])

ax_top.set_yticks([0, np.pi / 2, np.pi])
ax_top.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])

ax_top.set_ylabel(r'Static Field, $h_z T$')
dy_in = -0.4
dx_in = 0.05
label = ax_top.yaxis.label
label.set_transform(
    label.get_transform()
    + mtransforms.ScaledTranslation(dx_in, dy_in, ax_top.figure.dpi_scale_trans)
)

# =============================================================================
# 3) Bottom-right panel: 10^2 ΔS map
# =============================================================================

deloc = ax_bottom.imshow(
    100 * (Loc_Diff),
    origin='lower',
    extent=[hxs[0], hxs[-1], hzs[0] * T, hzs[-1] * T],
    aspect=0.5,
    vmin=-1,
    vmax=1,
    cmap=OWGcmap,
)

ax_bottom.set_xlabel(r'Kick Field, $h_x$', labelpad=2.2)

ax_bottom.set_yticks([0, np.pi / 2, np.pi])
ax_bottom.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])

ax_bottom.set_xticks([0, np.pi / 4, np.pi / 2])
ax_bottom.set_xticklabels([r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"])

# =============================================================================
# 4) Panel labels and annotations
# =============================================================================

fig.text(0.065, 0.94, r"(a)", fontsize=10)
fig.text(0.67,  0.94, r"(b)", fontsize=10)

vert_shift = 0.04
fig.text(
    0.83, 0.75 + vert_shift,
    r"$10^2 \frac{\langle \Delta \theta \rangle}{\pi}$",
    fontsize=10,
    bbox=dict(
        facecolor="white",
        alpha=0.5,
        edgecolor="none",
        boxstyle="round,pad=0.2",
    ),
)
fig.text(
    0.83, 0.395 + vert_shift,
    r"$10^2 \Delta S$",
    fontsize=10,
    bbox=dict(
        facecolor="white",
        alpha=0.5,
        edgecolor="none",
        boxstyle="round,pad=0.2",
    ),
)

# Colorbar inset matching width of ax_top
cax = inset_axes(
    ax_top,
    width="100%",
    height="5%",
    bbox_to_anchor=(0, 0.123123, 1, 1),
    bbox_transform=ax_top.transAxes,
    borderpad=0,
)

cbar = fig.colorbar(error, cax=cax, orientation="horizontal")
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels([r"$-1$", r"$0$", r"$1$"])
cbar.ax.xaxis.set_ticks_position("top")
cbar.ax.xaxis.set_label_position("top")
cbar.ax.tick_params(axis="x", pad=-0.5)

cax.set_in_layout(False)

# ============================================================================ #
# SAVE
# ============================================================================ #

plt.savefig(filename + ".pdf", dpi=600, bbox_inches=None)





