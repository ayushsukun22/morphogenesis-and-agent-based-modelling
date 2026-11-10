"""Summarize four matched intervention conditions from an explicit config."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def save_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    base = args.config.resolve().parent
    args.output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey="row", constrained_layout=True)
    rows, time_rows, contrasts = [], [], []
    inputs = {str(args.config): hashlib.sha256(args.config.read_bytes()).hexdigest()}
    colors = {"normal": "#3477ad", "protrusion_off": "#d97f25",
              "growth_off": "#389454", "both_off": "#9464ad"}
    for condition, paths in config["conditions"].items():
        analysis = base / paths["analysis"] / "lineage_intervals.csv"
        intervals = read_csv(analysis)
        inputs[str(analysis)] = hashlib.sha256(analysis.read_bytes()).hexdigest()
        for folder in paths["runs"]:
            folder = base / folder
            for csv_path in sorted(folder.glob("seed_*.csv")):
                seed = csv_path.stem.split("_")[1]
                inputs[str(csv_path)] = hashlib.sha256(csv_path.read_bytes()).hexdigest()
                observations = read_csv(csv_path)
                final = observations[-1]
                series = []
                for path in sorted(folder.glob(csv_path.stem + "_*.npz")):
                    with np.load(path) as data:
                        time = float(data["time"])
                        if time < config["onset"] - 1e-8:
                            continue
                        centered = data["positions"] - data["positions"].mean(axis=0)
                        rms = float(np.sqrt(np.mean(np.sum(centered ** 2, axis=1))))
                    series.append((time, rms))
                    inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
                subset = [r for r in intervals if r["seed"] == seed]
                timepoints = sorted({float(r["start_time"]) for r in subset})
                mean_velocity = []
                for time in timepoints:
                    same = [r for r in subset if float(r["start_time"]) == time]
                    value = float(np.mean([float(r["outward_velocity"]) for r in same]))
                    mean_velocity.append(value)
                    time_rows.append(dict(condition=condition, seed=seed, start_time=time,
                                          lineages=len(same), mean_outward_velocity=value))
                result = dict(condition=condition, seed=seed, cells=int(final["cells"]),
                              mean_outward_velocity=float(np.mean(mean_velocity)),
                              initial_rms_radius=series[0][1], final_rms_radius=series[-1][1],
                              rms_change_percent=100 * (series[-1][1] / series[0][1] - 1),
                              high_snail=int(final["high_snail_cells"]),
                              components=int(final["components"]),
                              largest_component_fraction=float(final["largest_component_fraction"]),
                              elongation=float(final["elongation"]),
                              outside_field=int(final["outside_field_cells"]))
                rows.append(result)
                column = config["seeds"].index(int(seed))
                axes[0, column].plot([t for t, _ in series],
                                     [100 * (r / series[0][1] - 1) for _, r in series],
                                     label=condition, color=colors[condition])
                axes[1, column].plot(timepoints, mean_velocity, marker="o",
                                     label=condition, color=colors[condition])
    for seed in config["seeds"]:
        group = {r["condition"]: r["mean_outward_velocity"] for r in rows if r["seed"] == str(seed)}
        on = group["normal"] - group["protrusion_off"]
        off = group["growth_off"] - group["both_off"]
        contrasts.append(dict(seed=seed, protrusion_contrast_growth_on=on,
                              protrusion_contrast_growth_off=off, interaction=on-off,
                              growth_contrast_protrusion_on=group["normal"]-group["growth_off"],
                              growth_contrast_protrusion_off=group["protrusion_off"]-group["both_off"]))
    for i, seed in enumerate(config["seeds"]):
        axes[0, i].set(title=f"Seed {seed}", xlabel="Model time (h)", ylabel="RMS radius change since hour 10 (%)")
        axes[1, i].set(xlabel="Interval start (model h)", ylabel="Mean outward velocity (model units/h)")
        for ax in axes[:, i]:
            ax.axhline(0, color="gray", linewidth=.7)
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Matched interventions begin at hour 10; same starting tissue within each seed")
    fig.savefig(args.output / "factorial_comparison.png", dpi=170)
    plt.close(fig)
    save_csv(args.output / "condition_summary.csv", rows)
    save_csv(args.output / "interval_summary.csv", time_rows)
    save_csv(args.output / "factorial_contrasts.csv", contrasts)
    (args.output / "report_provenance.json").write_text(json.dumps(dict(
        inputs=inputs, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2), encoding="utf-8")
    for row in rows:
        print(row)
    for row in contrasts:
        print(row)


if __name__ == "__main__":
    main()
