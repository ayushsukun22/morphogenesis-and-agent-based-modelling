"""Describe long fixed-population runs; peripheral screening is not branch detection."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from branching_baseline import measure
from branching_trajectories import geometry
from types import SimpleNamespace


def finite_mean(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(values.mean()) if len(values) else float("nan")


def window_metrics(frames, factor=1.25):
    """Frames must have equal IDs in equal order; return every cell's window metrics."""
    assert all(np.array_equal(frames[0]["ids"], f["ids"]) for f in frames)
    features = [geometry(f, factor) for f in frames]
    centered = np.array([x[0] for x in features])
    steps = np.linalg.norm(np.diff(centered, axis=0), axis=2)
    path = steps.sum(axis=0)
    displacement = np.linalg.norm(centered[-1] - centered[0], axis=1)
    straightness = np.divide(displacement, path, out=np.full_like(path, np.nan), where=path > 1e-12)
    high = np.all(np.array([f["snail"] > .8 for f in frames]), axis=0)
    min_rank = np.min([x[2] for x in features], axis=0)
    min_degree = np.min([x[3].sum(axis=1) for x in features], axis=0)
    before, after = features[0][3], features[-1][3]
    degree = before.sum(axis=1)
    retention = np.divide((before & after).sum(axis=1), degree,
                          out=np.full(len(degree), np.nan), where=degree > 0)
    return [dict(cell_id=int(cell_id), all_high_snail=bool(high[i]),
                 min_radial_rank=float(min_rank[i]), min_neighbors=int(min_degree[i]),
                 radial_change=float(features[-1][1][i] - features[0][1][i]),
                 net_displacement=float(displacement[i]), path_length=float(path[i]),
                 straightness=float(straightness[i]), neighbor_retention=float(retention[i]))
            for i, cell_id in enumerate(frames[0]["ids"])]


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("choose a fresh output folder")
    groups, hashes = {}, {}
    for folder in args.folders:
        for path in sorted(folder.glob("seed_*_*.npz")):
            with np.load(path) as data:
                if float(data["time"]) < 10 - 1e-8:
                    continue
                frame = {k: data[k].copy() for k in data.files}
            # All per-cell arrays are aligned; time is scalar.
            order = np.argsort(frame["ids"])
            frame = {k: v[order] if v.ndim else v for k, v in frame.items()}
            groups.setdefault(path.name.split("_")[1], []).append(frame)
            hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    if not groups:
        parser.error("no post-hour-10 snapshots found")
    summaries, times, windows, screens = [], [], [], []
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ordered = sorted(groups)
    extent = max(1., max(float(np.abs(f["positions"]).max()) for frames in groups.values() for f in frames)) * 1.05
    shape = plt.figure(figsize=(12, 12), constrained_layout=True)
    for column, seed in enumerate(ordered):
        frames = sorted(groups[seed], key=lambda f: float(f["time"]))
        assert all(np.array_equal(frames[0]["ids"], f["ids"]) for f in frames)
        assert all(np.array_equal(frames[0]["radii"], f["radii"]) for f in frames)
        series = []
        for f in frames:
            offsets, radial, rank, contact, _ = geometry(f, 1.25)
            high = f["snail"] > .8
            series.append(dict(seed=seed, time=float(f["time"]), cells=len(radial),
                               rms_radius=float(np.sqrt(np.mean(radial ** 2))),
                               high_snail=int(high.sum()), mean_rank_high=finite_mean(rank[high]),
                               high_peripheral=int(np.sum(high & (rank > .8))),
                               connected_high_peripheral=int(np.sum(high & (rank > .8) & (contact.sum(axis=1) > 0)))))
        times.extend(series)
        final = frames[-1]
        result = dict(seed=seed, **{k: v for k, v in series[-1].items() if k != "seed"},
                      rms_change_percent=100 * (series[-1]["rms_radius"] / series[0]["rms_radius"] - 1))
        cells = [SimpleNamespace(position=p, radius=float(r), grn=SimpleNamespace(snail=float(s)))
                 for p, r, s in zip(final["positions"], final["radii"], final["snail"])]
        for factor in [1., 1.25, 1.5]:
            m = measure(cells, factor)
            key = str(factor).replace(".", "_")
            result["components_" + key] = m["components"]
            result["largest_fraction_" + key] = m["largest_component_fraction"]
            result["elongation"] = m["elongation"]
        summaries.append(result)
        for index, start in enumerate(frames):
            if float(start["time"]) < 12 - 1e-8:
                continue
            end = index + 4
            if end >= len(frames):
                continue
            selected = frames[index:end + 1]
            if not np.isclose(float(selected[-1]["time"] - selected[0]["time"]), 2., atol=1e-8):
                raise ValueError("Expected half-hour snapshots spanning two hours")
            for row in window_metrics(selected):
                windows.append(dict(seed=seed, start_time=float(start["time"]),
                                    end_time=float(selected[-1]["time"]), **row))
        for threshold in [.7, .8, .9]:
            candidates = [r for r in windows if r["seed"] == seed and r["all_high_snail"]
                          and r["min_radial_rank"] > threshold and r["min_neighbors"] > 0]
            screens.append(dict(seed=seed, rank_threshold=threshold,
                                unique_cells=len({r["cell_id"] for r in candidates}),
                                cell_windows=len(candidates),
                                mean_radial_change=finite_mean([r["radial_change"] for r in candidates]),
                                mean_straightness=finite_mean([r["straightness"] for r in candidates]),
                                mean_neighbor_retention=finite_mean([r["neighbor_retention"] for r in candidates])))
        x = [r["time"] for r in series]
        for ax, key, ylabel in [(axes[0, 0], "high_snail", "High-SNAIL cells"),
                                (axes[0, 1], "mean_rank_high", "Mean radial rank of high-SNAIL cells"),
                                (axes[1, 1], "connected_high_peripheral", "High + peripheral + connected cells")]:
            ax.plot(x, [r[key] for r in series], label=f"Seed {seed}")
            ax.set(xlabel="Model time (h)", ylabel=ylabel)
            ax.grid(alpha=.2)
        axes[1, 0].plot(x, [100*(r["rms_radius"]/series[0]["rms_radius"]-1) for r in series], label=f"Seed {seed}")
        for row, target in enumerate([10., 15., 20.]):
            f = min(frames, key=lambda f: abs(float(f["time"])-target))
            ax = shape.add_subplot(3, len(ordered), row*len(ordered)+column+1, projection="3d")
            points = ax.scatter(*f["positions"].T, c=f["snail"], cmap="viridis", vmin=0, vmax=2, s=9, alpha=.75)
            ax.set(xlim=(-extent, extent), ylim=(-extent, extent), zlim=(-extent, extent),
                   title=f"Seed {seed}, {float(f['time']):.1f} h")
            ax.set_xticks([-5, 0, 5])
            ax.set_yticks([-5, 0, 5])
            ax.set_zticks([-5, 0, 5])
            ax.tick_params(labelsize=7, pad=0)
            ax.set_box_aspect((1, 1, 1))
    axes[1, 0].set(xlabel="Model time (h)", ylabel="RMS radius change since hour 10 (%)")
    axes[1, 0].grid(alpha=.2)
    axes[0, 1].set_ylim(0, 1)
    axes[0, 0].legend()
    fig.suptitle("Fixed population after hour 10; peripheral = radial rank >0.8")
    shape.suptitle("Cell centers in Cartesian model coordinates, equal scales; markers do not represent radii")
    shape.colorbar(points, ax=shape.axes, shrink=.5, label="SNAIL")
    args.output.mkdir(parents=True)
    fig.savefig(args.output / "timecourses.png", dpi=170)
    shape.savefig(args.output / "geometry.png", dpi=150)
    plt.close("all")
    for name, rows in [("final_summary.csv", summaries), ("timecourses.csv", times),
                       ("cell_windows.csv", windows), ("screen_summary.csv", screens)]:
        save_csv(args.output / name, rows)
    (args.output / "provenance.json").write_text(json.dumps(dict(inputs=hashes,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256={name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                       for name in ["branching_baseline.py", "branching_trajectories.py"]}), indent=2), encoding="utf-8")
    for row in summaries + screens:
        print(row)


if __name__ == "__main__":
    main()
