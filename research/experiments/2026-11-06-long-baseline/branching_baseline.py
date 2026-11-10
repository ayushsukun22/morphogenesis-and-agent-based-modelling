"""Run the existing cell model without Kit; export observations, not branch claims."""

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import types
from datetime import datetime, timezone
import platform

import numpy as np


def load_simulation():
    # Loading the regular package executes its Kit-only extension import.
    # Give the pure model its own package namespace instead.
    root = Path(__file__).resolve().parents[1]
    folder = root / "source/extensions/embryo_twin.sim/embryo_twin/sim"
    name = "_branching_model"
    package = types.ModuleType(name)
    package.__path__ = [str(folder)]
    sys.modules[name] = package
    spec = importlib.util.spec_from_file_location(name + ".simulation", folder / "simulation.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.EmbryoSimulation


def measure(cells, contact_factor=1.25):
    """Contact graph uses distance <= factor * sum of cell radii.

    This is a geometric proximity proxy, not persistent biological junctions.
    Graph junction degree is deliberately not interpreted as branch count:
    a compact tissue also has many high-degree vertices.
    """
    count = len(cells)
    if not count:
        return dict(cells=0, components=0, largest_component_fraction=0.0,
                    isolated_cells=0, elongation=0.0, high_snail_cells=0)
    positions = np.array([c.position for c in cells])
    radii = np.array([c.radius for c in cells])
    if not np.isfinite(positions).all() or not np.isfinite(radii).all():
        raise ValueError("Non-finite geometry")
    distances = np.linalg.norm(positions[:, None] - positions[None, :], axis=2)
    contact = distances <= contact_factor * (radii[:, None] + radii[None, :])
    np.fill_diagonal(contact, False)
    remaining = set(range(count))
    sizes = []
    while remaining:
        pending = [remaining.pop()]
        size = 0
        while pending:
            index = pending.pop()
            size += 1
            neighbors = remaining.intersection(np.flatnonzero(contact[index]).tolist())
            remaining.difference_update(neighbors)
            pending.extend(neighbors)
        sizes.append(size)
    centered = positions - positions.mean(axis=0)
    eigenvalues = np.linalg.eigvalsh(centered.T @ centered / count)
    # Radius regularization keeps a line or a single cell finite.
    floor = float(np.mean(radii ** 2)) / 5.0
    elongation = np.sqrt((eigenvalues[-1] + floor) / (max(eigenvalues[0], 0.0) + floor))
    return dict(cells=count, components=len(sizes),
                largest_component_fraction=max(sizes) / count,
                isolated_cells=int(np.sum(contact.sum(axis=1) == 0)),
                elongation=float(elongation),
                high_snail_cells=sum(c.grn.snail > 0.8 for c in cells))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hours", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 43, 44])
    parser.add_argument("--sample-every", type=float, default=1.0)
    parser.add_argument("--contact-factor", type=float, default=1.25)
    parser.add_argument("--output", type=Path, default=Path("_branching_baseline"))
    args = parser.parse_args()
    if not all(np.isfinite(x) and x > 0 for x in
               (args.hours, args.dt, args.sample_every, args.contact_factor)):
        parser.error("hours, dt, sample-every and contact-factor must be finite and positive")
    if any(seed < 0 for seed in args.seeds):
        parser.error("seeds must be nonnegative")
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("output folder is not empty; choose a new folder to preserve earlier experiments")
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = vars(args).copy()
    metadata["output"] = str(args.output)
    source = Path(__file__).resolve().parents[1] / "source/extensions/embryo_twin.sim/embryo_twin/sim"
    metadata.update(created_utc=datetime.now(timezone.utc).isoformat(),
                    python=platform.python_version(), numpy=np.__version__,
                    source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(source.glob("*.py"))},
                    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    stop_policy="requested horizon or application cell limit")
    (args.output / "settings.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    simulation_class = load_simulation()
    for seed in args.seeds:
        sim = simulation_class()
        sim.rng = np.random.default_rng(seed)
        sim._spawn_seed_cell()
        rows = []
        next_sample = 0.0
        snapshot = 0
        previous_ids = set(sim.cells)
        divisions_since_sample = 0
        peak_cells = 0
        csv_stream = (args.output / f"seed_{seed}.csv").open("w", newline="", encoding="utf-8")
        writer = None
        while True:
            at_limit = len(sim.cells) >= sim.MAX_CELLS
            final = sim.time >= args.hours - 1e-10 or at_limit
            if sim.time >= next_sample - 1e-10 or final:
                cells = sim.alive_cells()
                row = dict(time=sim.time, **measure(cells, args.contact_factor))
                values = lambda name: np.array([getattr(c.grn, name) for c in cells])
                pressure = np.array([c.pressure for c in cells])
                positions = np.array([c.position for c in cells]).reshape(-1, 3)
                row.update(max_snail=float(values("snail").max()) if cells else 0.0,
                           mean_snail=float(values("snail").mean()) if cells else 0.0,
                           mean_yap=float(values("yap").mean()) if cells else 0.0,
                           mean_notch=float(values("notch").mean()) if cells else 0.0,
                           mean_sox2=float(values("sox2").mean()) if cells else 0.0,
                           mean_force_proxy=float(pressure.mean()) if cells else 0.0,
                           divisions_since_sample=divisions_since_sample,
                           outside_field_cells=int(np.sum(np.any(
                               (positions < sim.bounds_min) | (positions > sim.bounds_max), axis=1))),
                           at_cell_limit=at_limit)
                rows.append(row)
                if writer is None:
                    writer = csv.DictWriter(csv_stream, fieldnames=list(row))
                    writer.writeheader()
                writer.writerow(row)
                csv_stream.flush()
                print(f"seed={seed} sample: {row}", flush=True)
                np.savez_compressed(args.output / f"seed_{seed}_{snapshot:04d}.npz",
                                    time=sim.time,
                                    ids=np.array([c.cell_id for c in cells]),
                                    positions=np.array([c.position for c in cells]).reshape(-1, 3),
                                    radii=np.array([c.radius for c in cells]),
                                    snail=values("snail"), yap=values("yap"),
                                    notch=values("notch"), sox2=values("sox2"),
                                    pressure=pressure,
                                    polarity=np.array([c.polarity for c in cells]).reshape(-1, 3),
                                    parent_ids=np.array([c.parent_id if c.parent_id is not None else -1 for c in cells]))
                snapshot += 1
                divisions_since_sample = 0
                next_sample += args.sample_every
            if final:
                break
            sim.update(min(args.dt, args.hours - sim.time))
            current_ids = set(sim.cells)
            divisions_since_sample += len(previous_ids - current_ids)
            previous_ids = current_ids
            peak_cells = max(peak_cells, len(current_ids))
        csv_stream.close()
        summary = dict(seed=seed, final_time=sim.time, peak_cells=peak_cells,
                       reason="cell_limit" if at_limit else "requested_horizon")
        (args.output / f"seed_{seed}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"seed={seed}: {rows[-1]}", flush=True)


if __name__ == "__main__":
    main()
