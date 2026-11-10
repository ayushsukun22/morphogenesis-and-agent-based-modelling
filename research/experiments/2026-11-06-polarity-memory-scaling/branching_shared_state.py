"""Capture a complete model state and run deterministic timestep diagnostics."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import platform
import csv

import numpy as np

from branching_baseline import load_simulation, install_growth_arrest, measure


def registry(sim):
    cell = next(iter(sim.cells.values()))
    return {type(x).__name__: type(x) for x in
            [sim, cell, cell.grn, sim.mechanics, sim.oxygen_field, sim.morphogen_field]}


def encode(value, classes):
    if isinstance(value, np.ndarray):
        return dict(kind="array", dtype=value.dtype.str, shape=value.shape, data=value.tolist())
    if isinstance(value, np.random.Generator):
        return dict(kind="rng", engine=type(value.bit_generator).__name__,
                    state=encode(value.bit_generator.state, classes))
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return dict(kind="dict", items=[[encode(k, classes), encode(v, classes)] for k, v in value.items()])
    if isinstance(value, (list, tuple)):
        return dict(kind=type(value).__name__, items=[encode(v, classes) for v in value])
    if type(value).__name__ in classes and type(value) is classes[type(value).__name__]:
        return dict(kind="object", name=type(value).__name__, state=encode(vars(value), classes))
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"Unsupported checkpoint value: {type(value)}")


def decode(value, classes):
    if not isinstance(value, dict):
        return value
    kind = value["kind"]
    if kind == "array":
        return np.array(value["data"], dtype=value["dtype"]).reshape(value["shape"])
    if kind == "dict":
        return {decode(k, classes): decode(v, classes) for k, v in value["items"]}
    if kind in ("list", "tuple"):
        items = [decode(v, classes) for v in value["items"]]
        return tuple(items) if kind == "tuple" else items
    if kind == "object":
        obj = object.__new__(classes[value["name"]])
        obj.__dict__.update(decode(value["state"], classes))
        return obj
    if kind == "rng":
        if value["engine"] != "PCG64":
            raise ValueError("Only the model's PCG64 generator is supported")
        engine = np.random.PCG64()
        engine.state = decode(value["state"], classes)
        return np.random.Generator(engine)
    raise ValueError(f"Unknown checkpoint kind: {kind}")


def canonical_state(sim, classes):
    return json.dumps(encode(sim, classes), allow_nan=False, separators=(",", ":"), sort_keys=True).encode()


def source_hashes():
    source = Path(__file__).resolve().parents[1] / "source/extensions/embryo_twin.sim/embryo_twin/sim"
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.glob("*.py"))}


class MidpointRandom:
    """Zero the existing 0.05*(random()-0.5) term; reject unplanned RNG use."""
    def __init__(self):
        self.calls = 0

    def random(self):
        self.calls += 1
        return 0.5


def scaled_memory(memory, dt, reference_dt):
    """Rescale scalar retention; vector renormalization remains unchanged."""
    if not np.isfinite(memory) or not 0 <= memory <= 1:
        raise ValueError("memory must be finite and between zero and one")
    if not all(np.isfinite(v) and v > 0 for v in [dt, reference_dt]):
        raise ValueError("timesteps must be finite and positive")
    return memory ** (dt / reference_dt)


def install_time_scaled_polarity(sim, reference_dt):
    original_update = sim.mechanics.update
    reference_memory = float(sim.mechanics.polarity_memory)
    # Validate before installing a wrapper.
    scaled_memory(reference_memory, reference_dt, reference_dt)
    sim.polarity_scaled_steps = 0

    def update(cells, dt):
        previous = sim.mechanics.polarity_memory
        try:
            sim.mechanics.polarity_memory = scaled_memory(reference_memory, dt, reference_dt)
            sim.polarity_scaled_steps += 1
            return original_update(cells, dt)
        finally:
            sim.mechanics.polarity_memory = previous

    sim.mechanics.update = update


def capture(args):
    sim = load_simulation()()
    sim.rng = np.random.default_rng(args.seed)
    sim._spawn_seed_cell()
    for _ in range(round(args.onset / args.warmup_dt)):
        sim.update(args.warmup_dt)
    classes = registry(sim)
    payload = canonical_state(sim, classes)
    restored = decode(json.loads(payload), classes)
    assert canonical_state(restored, classes) == payload
    assert restored.cells is not sim.cells
    assert not np.shares_memory(restored.oxygen_field.grid, sim.oxygen_field.grid)
    assert not np.shares_memory(restored.morphogen_field.grid, sim.morphogen_field.grid)
    assert all(not np.shares_memory(restored.cells[k].position, c.position) for k, c in sim.cells.items())
    # Test future stochastic equality without consuming the stored checkpoint RNG.
    assert np.array_equal(restored.rng.random(8), sim.rng.random(8))
    args.output.mkdir(parents=True)
    (args.output / "state.json.gz").write_bytes(gzip.compress(payload, mtime=0))
    metadata = dict(seed=args.seed, onset=args.onset, actual_time=restored.time,
                    warmup_dt=args.warmup_dt, cells=len(restored.cells),
                    state_sha256=hashlib.sha256(payload).hexdigest(), source_sha256=source_hashes(),
                    python=platform.python_version(), numpy=np.__version__,
                    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output / "settings.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("Captured complete model state; exact round-trip, array isolation and future RNG checks passed", metadata, flush=True)


def run(args):
    metadata = json.loads((args.checkpoint / "settings.json").read_text())
    if metadata["source_sha256"] != source_hashes():
        raise ValueError("Model source differs from checkpoint provenance")
    payload = gzip.decompress((args.checkpoint / "state.json.gz").read_bytes())
    if hashlib.sha256(payload).hexdigest() != metadata["state_sha256"]:
        raise ValueError("Checkpoint digest mismatch")
    prototype = load_simulation()()
    classes = registry(prototype)
    sim = decode(json.loads(payload), classes)
    assert canonical_state(sim, classes) == payload
    initial_ids = sorted(sim.cells)
    initial_radii = [sim.cells[k].radius for k in initial_ids]
    steps = round((args.horizon - metadata["onset"]) / args.dt)
    every = round(args.sample_every / args.dt)
    install_growth_arrest(sim, metadata["onset"])
    reference_memory = float(sim.mechanics.polarity_memory)
    if args.polarity_reference_dt is not None:
        install_time_scaled_polarity(sim, args.polarity_reference_dt)
    sim.rng = MidpointRandom()
    args.output.mkdir(parents=True)
    settings = dict(checkpoint=str(args.checkpoint), initial_state_sha256=metadata["state_sha256"],
                    dt=args.dt, horizon=args.horizon, sample_every=args.sample_every,
                    noise="GRN random() fixed at 0.5 after onset", growth="radius growth and division arrested",
                    polarity_reference_dt=args.polarity_reference_dt,
                    reference_polarity_memory=reference_memory,
                    applied_polarity_memory=(scaled_memory(reference_memory,args.dt,args.polarity_reference_dt)
                                             if args.polarity_reference_dt is not None else reference_memory),
                    source_sha256=source_hashes(), python=platform.python_version(), numpy=np.__version__,
                    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    helper_sha256=hashlib.sha256((Path(__file__).parent/"branching_baseline.py").read_bytes()).hexdigest())
    (args.output / "settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")
    with (args.output / "timecourses.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = None
        for step in range(steps + 1):
            if step % every == 0 or step == steps:
                cells = [sim.cells[k] for k in initial_ids]
                assert sorted(sim.cells) == initial_ids
                np.testing.assert_array_equal([c.radius for c in cells], initial_radii)
                positions = np.array([c.position for c in cells])
                genes = np.array([[getattr(c.grn, key) for key in
                                   ["sox2", "snail", "yap", "notch", "cell_cycle", "hif1a"]] for c in cells])
                arrays = dict(ids=np.array(initial_ids), positions=positions, radii=np.array(initial_radii),
                              genes=genes, polarity=np.array([c.polarity for c in cells]),
                              pressure=np.array([c.pressure for c in cells]),
                              oxygen_grid=sim.oxygen_field.grid, morphogen_grid=sim.morphogen_field.grid)
                assert all(np.isfinite(a).all() for a in arrays.values())
                row = dict(time=sim.time, steps=step, **measure(cells),
                           rms_radius=float(np.sqrt(np.mean(np.sum((positions-positions.mean(axis=0))**2,axis=1)))),
                           mean_snail=float(genes[:,1].mean()),
                           mean_force=float(arrays["pressure"].mean()),
                           polarity_scaled_steps=getattr(sim,"polarity_scaled_steps",0),
                           outside_field=int(np.any((positions<sim.bounds_min)|(positions>sim.bounds_max),axis=1).sum()))
                if writer is None:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                writer.writerow(row)
                stream.flush()
                np.savez_compressed(args.output / f"sample_{round(step/every):04d}.npz", time=sim.time, **arrays)
                print(f"dt={args.dt}: {row}", flush=True)
            if step < steps:
                sim.update(args.dt)
    assert sim.growth_arrest_steps == steps
    assert sim.rng.calls == steps * len(initial_ids)
    if args.polarity_reference_dt is not None:
        assert sim.polarity_scaled_steps == steps
    summary = dict(final_time=sim.time, steps=steps, cells=len(sim.cells),
                   growth_arrest_steps=sim.growth_arrest_steps, midpoint_calls=sim.rng.calls,
                   polarity_scaled_steps=getattr(sim,"polarity_scaled_steps",0),
                   initial_state_sha256=metadata["state_sha256"])
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("capture")
    make.add_argument("--seed", type=int, default=44)
    make.add_argument("--onset", type=float, default=10.)
    make.add_argument("--warmup-dt", type=float, default=.05)
    make.add_argument("--output", type=Path, required=True)
    arm = sub.add_parser("run")
    arm.add_argument("--checkpoint", type=Path, required=True)
    arm.add_argument("--dt", type=float, required=True)
    arm.add_argument("--horizon", type=float, default=20.)
    arm.add_argument("--sample-every", type=float, default=.5)
    arm.add_argument("--polarity-reference-dt", type=float, default=None,
                     help="Experimental retention scaling: memory_ref**(dt/reference_dt)")
    arm.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output folder")
    if args.command == "capture":
        if args.seed < 0 or not all(np.isfinite(v) and v > 0 for v in [args.onset,args.warmup_dt]):
            parser.error("seed must be nonnegative; onset and warmup-dt must be positive and finite")
        if not np.isclose(args.onset/args.warmup_dt, round(args.onset/args.warmup_dt)):
            parser.error("onset must align with warmup-dt")
        capture(args)
    else:
        meta = json.loads((args.checkpoint/"settings.json").read_text())
        duration = args.horizon-meta["onset"]
        if not all(np.isfinite(v) and v > 0 for v in [args.dt,args.sample_every,duration]):
            parser.error("dt, sample interval and post-onset duration must be positive and finite")
        if args.polarity_reference_dt is not None and not (np.isfinite(args.polarity_reference_dt) and args.polarity_reference_dt > 0):
            parser.error("polarity-reference-dt must be finite and positive")
        if not all(np.isclose(v/args.dt, round(v/args.dt), rtol=0,atol=1e-8) for v in [duration,args.sample_every]):
            parser.error("duration and sample interval must align with dt")
        run(args)


if __name__ == "__main__":
    main()
