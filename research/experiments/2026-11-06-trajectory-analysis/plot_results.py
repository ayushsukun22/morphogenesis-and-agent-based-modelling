"""Rebuild descriptive geometry figure from this folder's archived CSV."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
with (root / "cell_geometry.csv").open(newline="", encoding="utf-8") as stream:
    rows = [r for r in csv.DictReader(stream) if float(r["time"]) == 12.0]
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for ax, (key, ylabel) in zip(axes, [
        ("radial_rank", "Radial rank (0 = central; 1 = peripheral)"),
        ("neighbors", "Proximity neighbors (factor 1.25)")]):
    for high, offset, label, color in [(False, -0.17, "Other cells", "#4682b4"),
                                       (True, 0.17, "High SNAIL (>0.8)", "#d87524")]:
        data = [[float(r[key]) for r in rows if r["seed"] == seed
                 and r["high_snail"] == str(high)] for seed in ["42", "43", "44"]]
        plots = ax.boxplot(data, positions=np.arange(3) + offset, widths=0.27,
                          patch_artist=True, manage_ticks=False, showfliers=False,
                          medianprops={"color": "black"})
        for patch in plots["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.65)
        ax.plot([], [], color=color, linewidth=7, label=label)
    ax.set(xticks=np.arange(3), xticklabels=["42", "43", "44"], xlabel="Seed", ylabel=ylabel)
    ax.grid(axis="y", alpha=0.2)
axes[0].set_ylim(0, 1)
axes[0].legend(fontsize=8)
fig.suptitle("Hour 12: high-SNAIL cells tend to be more central and have more neighbors")
fig.savefig(root / "cell_geometry.png", dpi=180)
