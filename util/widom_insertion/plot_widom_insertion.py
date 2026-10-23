"""
Plot Widom insertion results from widom_insertion.type_<X>.dat.

Usage:
    python3 plot_widom.py [datafile] [temperature_K]

Defaults:
    datafile      = widom_insertion.type_C.dat
    temperature_K = config.TEMP (read from config.py in the current directory)

Expected columns:
    frame iter x_new y_new z_new boltzmann_factor weighted_boltzmann_factor

Plots use the weighted Boltzmann factor (last column), which is the
quantity averaged to get the chemical potential.

Output:
    widom_insertion_convergence.png  three side-by-side panels:
        1. per-frame running averages, all frames overlaid
        2. histogram of per-frame chemical potentials
        3. running average over all insertions (insertion-major order)
"""

import os
import sys
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.normpath(os.getcwd()))

datafile = sys.argv[1] if len(sys.argv) > 1 else "widom_insertion.type_C.dat"

if len(sys.argv) > 2:
    temperature = float(sys.argv[2])
else:
    import config
    temperature = config.TEMP

KB = 8.617333262e-5 * 23.06    # kcal/mol/K, must match widom_insertion.py
BETA = 1.0 / (KB * temperature)

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------
data = np.loadtxt(datafile, comments="#")

frame_col = data[:, 0].astype(int)
iter_col = data[:, 1].astype(int)
all_w = data[:, 6]

frame_ids = np.unique(frame_col)
nframes = len(frame_ids)

per_frame_weights = [all_w[frame_col == f] for f in frame_ids]
per_frame_averages = np.array([np.mean(w) for w in per_frame_weights])

overall_avg = np.mean(all_w)
overall_mu = -(1.0 / BETA) * np.log(overall_avg)

# Per-frame chemical potential (frames with non-positive average are skipped)
per_frame_mu = np.full(nframes, np.nan)
pos = per_frame_averages > 0.0
per_frame_mu[pos] = -(1.0 / BETA) * np.log(per_frame_averages[pos])
per_frame_mu = per_frame_mu[np.isfinite(per_frame_mu)]

print(f"Read {len(all_w)} insertions from {nframes} frames.")

# ---------------------------------------------------------
# Running average over all insertions, ordered insertion-major:
# insertion 1 of frame 1, insertion 1 of frame 2, ..., insertion 1 of the
# last frame, then insertion 2 of frame 1, and so on.
# ---------------------------------------------------------
order = np.lexsort((frame_col, iter_col))      # sort by iter, then frame
ordered_weights = all_w[order]
running_insertion_average = (
    np.cumsum(ordered_weights) / np.arange(1, len(ordered_weights) + 1)
)

# ---------------------------------------------------------
# Single figure with three side-by-side panels
# ---------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

# 1. All per-frame running averages
for weights in per_frame_weights:
    n_ins = np.arange(1, len(weights) + 1)
    ax1.plot(n_ins, np.cumsum(weights) / n_ins, alpha=0.4, linewidth=0.8)
ax1.axhline(overall_avg, color="k", linestyle="--", label="Overall average")
ax1.set_xlabel("Number of insertions")
ax1.set_ylabel(r"Running $\langle e^{-\beta \Delta E}\rangle$")
ax1.set_title("Per-frame running averages")
ax1.legend()

# 2. Histogram of per-frame chemical potentials
ax2.hist(per_frame_mu, bins=30)
ax2.axvline(overall_mu, color="k", linestyle="--", label="Overall estimate")
ax2.set_xlabel("Per-frame chemical potential")
ax2.set_ylabel("Count (frames)")
ax2.set_title("Per-frame chemical potentials")
ax2.legend()

# 3. Running average over all insertions
ax3.plot(np.arange(1, len(ordered_weights) + 1), running_insertion_average)
ax3.axhline(overall_avg, linestyle="--", color="k", label="Overall average")
ax3.set_xlabel(f"Number of insertions (cycling through {nframes} frames)")
ax3.set_ylabel(r"$\langle e^{-\beta \Delta E}\rangle$")
ax3.set_title("Running average, all insertions")
ax3.legend()

fig.tight_layout()
fig.savefig("widom_insertion_convergence.png", dpi=200)
plt.close(fig)

print("\nPlot saved: widom_insertion_convergence.png")