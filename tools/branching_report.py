"""Plot baseline CSVs and final cell centers without inferring branch counts."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for folder in args.folders:
        for path in sorted(folder.glob("seed_*.csv")):
            with path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
            records.append((folder, path.stem, rows))
    if not records:
        parser.error("no seed CSV files found")
    args.output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    measures = [("cells", "Cell count"), ("max_snail", "Maximum SNAIL"),
                ("largest_component_fraction", "Largest proximity component / cells"),
                ("elongation", "Elongation (not branch count)")]
    for _, name, rows in records:
        time = [float(row["time"]) for row in rows]
        for ax, (key, label) in zip(axes.flat, measures):
            ax.plot(time, [float(row[key]) for row in rows], label=name)
            ax.set(xlabel="Model time (h)", ylabel=label)
            ax.grid(alpha=0.2)
    axes[0, 1].axhline(0.8, color="black", linestyle="--", linewidth=1, label="SNAIL threshold")
    axes[0, 0].legend()
    axes[0, 1].legend()
    fig.suptitle("Baseline observations (experiment settings accompany source CSVs)")
    fig.savefig(args.output / "timecourses.png", dpi=170)
    plt.close(fig)

    snapshots = []
    for folder, name, _ in records:
        path = sorted(folder.glob(name + "_*.npz"))[-1]
        with np.load(path) as data:
            snapshots.append((name, float(data["time"]), data["positions"].copy(), data["snail"].copy()))
    extent = max(1.0, max(float(np.abs(p).max()) for _, _, p, _ in snapshots)) * 1.05
    fig = plt.figure(figsize=(5 * len(snapshots), 4.5), constrained_layout=True)
    for index, (name, time, positions, snail) in enumerate(snapshots):
        ax = fig.add_subplot(1, len(snapshots), index + 1, projection="3d")
        points = ax.scatter(*positions.T, c=snail, cmap="viridis", vmin=0, vmax=2, s=8, alpha=0.7)
        ax.set(xlim=(-extent, extent), ylim=(-extent, extent), zlim=(-extent, extent),
               xlabel="x (model units)", ylabel="y", zlabel="z",
               title=f"{name}, {time:.2f} h, n={len(positions)}")
        ax.set_box_aspect((1, 1, 1))
        fig.colorbar(points, ax=ax, shrink=0.6, label="SNAIL")
    fig.suptitle("Final cell centers; equal axis scales; markers do not depict cell radii")
    fig.savefig(args.output / "final_geometry.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
