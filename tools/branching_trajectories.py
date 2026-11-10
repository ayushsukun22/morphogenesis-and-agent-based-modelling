"""Descriptive, lineage-aware reanalysis of stored baseline observations."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def geometry(frame, factor):
    positions = frame["positions"]
    offsets = positions - positions.mean(axis=0)
    radial = np.linalg.norm(offsets, axis=1)
    delta = positions[None, :, :] - positions[:, None, :]
    distance = np.linalg.norm(delta, axis=2)
    contact = distance <= factor * (frame["radii"][:, None] + frame["radii"][None, :])
    np.fill_diagonal(contact, False)
    directions = delta / np.maximum(distance[:, :, None], 1e-12)
    degree = contact.sum(axis=1)
    exposure = np.linalg.norm((directions * contact[:, :, None]).sum(axis=1), axis=1)
    exposure /= np.maximum(degree, 1)
    exposure[degree == 0] = np.nan  # isolated is not proof of being a tissue tip
    # Rank within this frame avoids equating absolute radius with surface.
    rank = np.array([(np.sum(radial < x) + 0.5 * np.sum(radial == x)) / len(radial)
                     for x in radial])
    return offsets, radial, rank, contact, exposure


def ancestry_to_start(end_ids, start_ids, parents):
    """Map observed endpoint cells to a start-frame ancestor when identifiable."""
    result = []
    for cell_id in end_ids:
        current = int(cell_id)
        seen = set()
        while current not in start_ids and current not in seen and current in parents:
            seen.add(current)
            current = parents[current]
        result.append(current if current in start_ids else -1)
    return np.array(result)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contact-factor", type=float, default=1.25)
    args = parser.parse_args()
    if not np.isfinite(args.contact_factor) or args.contact_factor <= 0:
        parser.error("contact-factor must be finite and positive")
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("choose a fresh output folder")
    paths = [p for folder in args.folders for p in sorted(folder.glob("seed_*_*.npz"))]
    if not paths:
        parser.error("no snapshots found")
    groups = {}
    for path in paths:
        with np.load(path) as data:
            frame = {k: data[k].copy() for k in data.files}
        groups.setdefault((str(path.parent), path.name.split("_")[1]), []).append(frame)
    frame_rows, track_rows, summaries = [], [], []
    for (_, seed), frames in groups.items():
        frames.sort(key=lambda f: float(f["time"]))
        parents = {}
        for frame in frames:
            parents.update({int(i): int(p) for i, p in zip(frame["ids"], frame["parent_ids"]) if p >= 0})
        for frame in frames:
            time = float(frame["time"])
            if time < 9.99:
                continue
            _, radial, rank, contact, exposure = geometry(frame, args.contact_factor)
            for i, cell_id in enumerate(frame["ids"]):
                frame_rows.append(dict(seed=seed, time=time, cell_id=int(cell_id),
                                       high_snail=bool(frame["snail"][i] > 0.8),
                                       snail=float(frame["snail"][i]), radial_distance=float(radial[i]),
                                       radial_rank=float(rank[i]), neighbors=int(contact[i].sum()),
                                       local_directional_imbalance=float(exposure[i])))
        for start, end in zip(frames, frames[1:]):
            if float(start["time"]) < 9.99:
                continue
            offsets, _, _, contact, _ = geometry(start, args.contact_factor)
            _, _, _, end_contact, _ = geometry(end, args.contact_factor)
            start_ids = set(map(int, start["ids"]))
            ancestors = ancestry_to_start(end["ids"], start_ids, parents)
            dt = float(end["time"] - start["time"])
            drift = end["positions"].mean(axis=0) - start["positions"].mean(axis=0)
            unmatched = int(np.sum(ancestors == -1))
            for i, cell_id in enumerate(start["ids"]):
                descendants = np.flatnonzero(ancestors == cell_id)
                if not len(descendants):
                    continue
                endpoint = np.average(end["positions"][descendants], axis=0,
                                      weights=end["radii"][descendants] ** 3)
                displacement = endpoint - start["positions"][i] - drift
                radial_unit = offsets[i] / max(np.linalg.norm(offsets[i]), 1e-12)
                old_neighbors = set(map(int, start["ids"][contact[i]]))
                adjacent = np.any(end_contact[descendants], axis=0)
                retained = set(map(int, ancestors[adjacent])) & old_neighbors
                retention = len(retained) / len(old_neighbors) if old_neighbors else np.nan
                track_rows.append(dict(seed=seed, start_time=float(start["time"]),
                                       end_time=float(end["time"]), cell_id=int(cell_id),
                                       high_at_start=bool(start["snail"][i] > 0.8),
                                       high_descendant_fraction=float(np.mean(end["snail"][descendants] > 0.8)),
                                       descendants=len(descendants), unmatched_end_cells=unmatched,
                                       outward_velocity=float(displacement @ radial_unit / dt),
                                       displacement_speed=float(np.linalg.norm(displacement) / dt),
                                       neighbor_retention=float(retention)))
        last_time = float(frames[-1]["time"])
        final = [r for r in frame_rows if r["seed"] == seed and r["time"] == last_time]
        for high in (False, True):
            subset = [r for r in final if r["high_snail"] == high]
            summaries.append(dict(seed=seed, time=last_time, high_snail=high, cells=len(subset),
                                  mean_radial_rank=float(np.mean([r["radial_rank"] for r in subset])),
                                  mean_neighbors=float(np.mean([r["neighbors"] for r in subset])),
                                  mean_local_imbalance=float(np.nanmean([r["local_directional_imbalance"] for r in subset]))))
    args.output.mkdir(parents=True)
    write_csv(args.output / "cell_geometry.csv", frame_rows)
    write_csv(args.output / "lineage_intervals.csv", track_rows)
    write_csv(args.output / "final_groups.csv", summaries)
    provenance = dict(contact_factor=args.contact_factor, snail_threshold=0.8,
                      script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      inputs={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    for row in summaries:
        print(row)
    for seed in sorted({r["seed"] for r in track_rows}):
        for high in (False, True):
            subset = [r for r in track_rows if r["seed"] == seed and r["high_at_start"] == high]
            print(dict(seed=seed, high_at_start=high, intervals=len(subset),
                       mean_outward_velocity=float(np.mean([r["outward_velocity"] for r in subset])) if subset else None,
                       mean_neighbor_retention=float(np.nanmean([r["neighbor_retention"] for r in subset])) if subset else None))


if __name__ == "__main__":
    main()
