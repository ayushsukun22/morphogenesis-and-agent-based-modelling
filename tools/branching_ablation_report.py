"""Paired descriptive summary of delayed force ablation and archived controls."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", type=Path, nargs="+", required=True)
    parser.add_argument("--treatment", type=Path, required=True)
    parser.add_argument("--control-analysis", type=Path, required=True)
    parser.add_argument("--treatment-analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    control_paths = {p.name: p for folder in args.controls for p in folder.glob("seed_*.csv")}
    conditions = {"control": read(args.control_analysis / "lineage_intervals.csv"),
                  "force_off": read(args.treatment_analysis / "lineage_intervals.csv")}
    paired = []
    intervals = []
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    colors = {"42": "#3575af", "43": "#d57824", "44": "#348d50"}
    for treatment_path in sorted(args.treatment.glob("seed_*.csv")):
        seed = treatment_path.stem.split("_")[1]
        control_path = control_paths[treatment_path.name]
        final = {"control": read(control_path)[-1], "force_off": read(treatment_path)[-1]}
        result = dict(seed=seed)
        for condition, csv_path in [("control", control_path), ("force_off", treatment_path)]:
            row = final[condition]
            data = np.load(csv_path.with_name(csv_path.stem + "_0024.npz"))
            centered = data["positions"] - data["positions"].mean(axis=0)
            rms = float(np.sqrt(np.mean(np.sum(centered ** 2, axis=1))))
            result[condition + "_rms_radius"] = rms
            for key in ["cells", "high_snail_cells", "components", "largest_component_fraction", "elongation", "outside_field_cells"]:
                result[condition + "_" + key] = float(row[key])
            subset = [r for r in conditions[condition] if r["seed"] == seed]
            timepoints = sorted({float(r["start_time"]) for r in subset})
            velocity = []
            for time in timepoints:
                same_interval = [r for r in subset if float(r["start_time"]) == time]
                value = float(np.mean([float(r["outward_velocity"]) for r in same_interval]))
                velocity.append(value)
                intervals.append(dict(seed=seed, condition=condition, start_time=time,
                                      start_lineages=len(same_interval), mean_outward_velocity=value))
            result[condition + "_mean_interval_velocity"] = float(np.mean(velocity))
            ax = axes[int(seed) - 42]
            ax.plot(timepoints, velocity, marker="o", color=colors[seed],
                    linestyle="-" if condition == "control" else "--", label=condition)
            ax.set(title=f"Seed {seed}", xlabel="Interval start (model h)",
                   ylabel="Mean outward velocity (model units/h)")
            ax.axhline(0, color="gray", linewidth=0.7)
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        result["velocity_difference_off_minus_control"] = result["force_off_mean_interval_velocity"] - result["control_mean_interval_velocity"]
        paired.append(result)
        print(result)
    for name, rows in [("paired_summary.csv", paired), ("interval_summary.csv", intervals)]:
        with (args.output / name).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    fig.suptitle("Protrusion force switched off at hour 10; all starting lineages per interval")
    fig.savefig(args.output / "paired_motion.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
