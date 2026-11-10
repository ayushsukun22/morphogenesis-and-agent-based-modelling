"""Rebuild exploratory timestep comparison from the archived seed-44 runs."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
summary = []
for folder, dt in [("seed44", .05), ("seed44_dt0025", .025)]:
    base = root / folder
    with (base / "seed_44.csv").open(newline="", encoding="utf-8") as stream:
        rows = [r for r in csv.DictReader(stream) if float(r["time"]) > 10 - 1e-8]
    time = [float(r["time"]) for r in rows]
    rms = []
    for path in sorted(base.glob("*.npz"))[20:]:
        with np.load(path) as data:
            p = data["positions"] - data["positions"].mean(axis=0)
            rms.append(float(np.sqrt(np.mean(np.sum(p*p, axis=1)))))
    for ax, key, ylabel in [(axes[0, 0], "high_snail_cells", "High-SNAIL cells"),
                            (axes[1, 0], "mean_force_proxy", "Mean net-force magnitude"),
                            (axes[1, 1], "elongation", "Elongation")]:
        ax.plot(time, [float(r[key]) for r in rows], label=f"dt={dt}")
        ax.set(xlabel="Model time (h)", ylabel=ylabel)
        ax.grid(alpha=.2)
    change = 100 * (np.asarray(rms)/rms[0]-1)
    axes[0, 1].plot(time, change, label=f"dt={dt}")
    summary.append(dict(seed=44, dt=dt, initial_rms_radius=rms[0], final_rms_radius=rms[-1],
                        rms_change_percent=float(change[-1]),
                        high_snail_at_20=int(rows[-1]["high_snail_cells"]),
                        peak_high_snail=max(int(r["high_snail_cells"]) for r in rows),
                        maximum_saved_mean_force=max(float(r["mean_force_proxy"]) for r in rows),
                        outside_field_max=max(int(r["outside_field_cells"]) for r in rows)))
axes[0, 0].legend()
axes[0, 1].set(xlabel="Model time (h)", ylabel="RMS radius change since hour 10 (%)")
axes[0, 1].grid(alpha=.2)
fig.suptitle("Exploratory timestep sensitivity: different random draws and hour-10 states")
fig.savefig(root / "timestep_sensitivity.png", dpi=170)
with (root / "timestep_summary.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)
for row in summary:
    print(row)
