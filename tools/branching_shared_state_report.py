"""Compare matched-state timestep arms without calling the finest arm truth."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output folder")
    arms = []
    inputs = {}
    for folder in args.folders:
        settings = json.loads((folder/"settings.json").read_text())
        summary = json.loads((folder/"summary.json").read_text())
        frames = []
        for path in sorted(folder.glob("sample_*.npz")):
            with np.load(path) as data:
                frames.append({k: data[k].copy() for k in data.files})
            inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        for name in ["settings.json", "summary.json", "timecourses.csv"]:
            path = folder/name
            inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        arms.append(dict(dt=settings["dt"], settings=settings, summary=summary,
                         frames=frames, rows=read(folder/"timecourses.csv")))
    arms.sort(key=lambda arm: arm["dt"], reverse=True)
    assert len({a["settings"]["initial_state_sha256"] for a in arms}) == 1
    assert all(a["settings"]["source_sha256"] == arms[0]["settings"]["source_sha256"] for a in arms)
    initial = arms[0]["frames"][0]
    assert all(np.array_equal(initial[k], arm["frames"][0][k]) for arm in arms for k in initial)
    pairs, finals = [], []
    for i, coarse in enumerate(arms[:-1]):
        for fine in arms[i+1:]:
            assert len(coarse["frames"]) == len(fine["frames"])
            for a, b in zip(coarse["frames"], fine["frames"]):
                assert abs(float(a["time"]-b["time"])) < 1e-8
                np.testing.assert_array_equal(a["ids"], b["ids"])
                error = float(np.sqrt(np.mean(np.sum((a["positions"]-b["positions"])**2, axis=1))))
                pairs.append(dict(coarse_dt=coarse["dt"], fine_dt=fine["dt"], time=float(b["time"]),
                                  position_rms_distance=error,
                                  position_error_in_mean_radii=error/float(b["radii"].mean()),
                                  grn_rmse=float(np.sqrt(np.mean((a["genes"]-b["genes"])**2))),
                                  pressure_rmse=float(np.sqrt(np.mean((a["pressure"]-b["pressure"])**2))),
                                  polarity_vector_rmse=float(np.sqrt(np.mean((a["polarity"]-b["polarity"])**2))),
                                  snail_rmse=float(np.sqrt(np.mean((a["genes"][:,1]-b["genes"][:,1])**2))),
                                  oxygen_rmse=float(np.sqrt(np.mean((a["oxygen_grid"]-b["oxygen_grid"])**2))),
                                  morphogen_rmse=float(np.sqrt(np.mean((a["morphogen_grid"]-b["morphogen_grid"])**2)))))
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for arm in arms:
        rows = arm["rows"]
        time = [float(r["time"]) for r in rows]
        for ax, key, label in [(axes[0,0], "rms_radius", "RMS radius (model units)"),
                               (axes[0,1], "high_snail_cells", "High-SNAIL cells")]:
            ax.plot(time, [float(r[key]) for r in rows], label=f"dt={arm['dt']}")
            ax.set(xlabel="Model time (h)", ylabel=label)
            ax.grid(alpha=.2)
        last = rows[-1]
        finals.append(dict(dt=arm["dt"], final_time=float(last["time"]),
                           rms_radius=float(last["rms_radius"]), high_snail=int(last["high_snail_cells"]),
                           mean_snail=float(last["mean_snail"]), mean_force=float(last["mean_force"]),
                           elongation=float(last["elongation"]), components=int(last["components"]),
                           largest_component_fraction=float(last["largest_component_fraction"]),
                           outside_field_max=max(int(r["outside_field"]) for r in rows),
                           steps=arm["summary"]["steps"]))
    for i, coarse in enumerate(arms[:-1]):
        for fine in arms[i+1:]:
            subset = [r for r in pairs if r["coarse_dt"] == coarse["dt"] and r["fine_dt"] == fine["dt"]]
            for ax, key, label in [(axes[1,0], "position_rms_distance", "Per-cell position RMS difference"),
                                   (axes[1,1], "grn_rmse", "Six-node GRN RMSE")]:
                ax.plot([r["time"] for r in subset], [r[key] for r in subset],
                        label=f"{coarse['dt']} vs {fine['dt']}")
                ax.set(xlabel="Model time (h)", ylabel=label)
                ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    axes[1,0].legend(fontsize=8)
    fig.suptitle("Identical complete hour-10 state; post-onset GRN noise zero; growth arrested")
    args.output.mkdir(parents=True)
    fig.savefig(args.output/"timestep_comparison.png", dpi=170)
    plt.close(fig)
    write(args.output/"final_summary.csv", finals)
    write(args.output/"pairwise_errors.csv", pairs)
    provenance = dict(inputs=inputs, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    for row in finals:
        print(row)
    for row in pairs:
        if abs(row["time"]-20)<1e-8:
            print(row)


if __name__ == "__main__":
    main()
