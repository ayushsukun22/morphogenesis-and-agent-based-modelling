"""Refine passive pair mechanics; validate translational equivalence to E010."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from branching_baseline import load_simulation,measure
from branching_shared_state import decode,registry,canonical_state,source_hashes,regulatory_arrays
from branching_passive_report import pair_potential


def pair_step(cells,mechanics,dt):
    """Original pair-force order and arithmetic, with zero protrusive force.

    Omit neighbor-signaling/polarity work: with frozen outputs and no
    protrusion those calculations cannot affect translation.
    """
    if len(cells)<2:
        return
    forces={cell.cell_id:np.zeros(3,dtype=np.float64) for cell in cells}
    for i in range(len(cells)):
        cell_a=cells[i]
        for j in range(i+1,len(cells)):
            cell_b=cells[j]
            delta=cell_b.position-cell_a.position
            distance=np.linalg.norm(delta)
            if distance<1e-8:
                continue
            direction=delta/distance
            contact_distance=cell_a.radius+cell_b.radius
            if distance>mechanics.interaction_radius:
                continue
            overlap=contact_distance-distance
            if overlap>0.:
                repulsion=mechanics.repulsion_strength*overlap
                forces[cell_a.cell_id]-=direction*repulsion
                forces[cell_b.cell_id]+=direction*repulsion
            adhesion=(cell_a.adhesion_strength*cell_b.adhesion_strength
                      *min(1.,mechanics.interaction_radius/max(distance,1e-8)))
            forces[cell_a.cell_id]+=direction*adhesion
            forces[cell_b.cell_id]-=direction*adhesion
    for cell in cells:
        force=forces[cell.cell_id]
        cell.pressure=float(np.linalg.norm(force))
        velocity=force/max(mechanics.damping,1e-8)
        cell.velocity=velocity
        cell.position+=velocity*dt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint",type=Path,required=True)
    parser.add_argument("--dt",type=float,required=True)
    parser.add_argument("--horizon",type=float,default=20.)
    parser.add_argument("--sample-every",type=float,default=.5)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output folder")
    meta=json.loads((args.checkpoint/"settings.json").read_text())
    duration=args.horizon-meta["onset"]
    if not all(np.isfinite(x) and x>0 for x in [args.dt,args.sample_every,duration]):
        parser.error("dt, sample interval and duration must be finite and positive")
    if not all(np.isclose(x/args.dt,round(x/args.dt),rtol=0,atol=1e-8) for x in [duration,args.sample_every]):
        parser.error("duration and sample interval must align with dt")
    if meta["source_sha256"]!=source_hashes():
        raise ValueError("Application source mismatch")
    payload=gzip.decompress((args.checkpoint/"state.json.gz").read_bytes())
    assert hashlib.sha256(payload).hexdigest()==meta["state_sha256"]
    classes=registry(load_simulation()())
    sim=decode(json.loads(payload),classes)
    assert canonical_state(sim,classes)==payload
    cells=sim.alive_cells()  # preserve original force accumulation order
    output_cells=sorted(cells,key=lambda c:c.cell_id)
    initial_ids=np.array([c.cell_id for c in output_cells])
    radii=np.array([c.radius for c in output_cells])
    adhesion=np.array([c.adhesion_strength for c in output_cells])
    nodes,outputs=regulatory_arrays(output_cells)
    initial_center=np.mean([c.position for c in output_cells],axis=0)
    steps=round(duration/args.dt)
    every=round(args.sample_every/args.dt)
    args.output.mkdir(parents=True)
    helpers=["branching_baseline.py","branching_shared_state.py","branching_passive_report.py"]
    settings=dict(dt=args.dt,horizon=args.horizon,sample_every=args.sample_every,
                  initial_state_sha256=meta["state_sha256"],source_sha256=source_hashes(),
                  force_order_ids=[c.cell_id for c in cells],
                  protocol="pair mechanics only; fixed regulation/radii; no protrusion/noise; fields, age, fate and polarity not advanced",
                  python=platform.python_version(),numpy=np.__version__,
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  helper_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in helpers})
    (args.output/"settings.json").write_text(json.dumps(settings,indent=2),encoding="utf-8")
    max_increase=0.
    increases=0
    near_zero_max=0
    drift_max=0.
    previous=None
    with (args.output/"step_diagnostics.csv").open("w",newline="",encoding="utf-8") as trace, \
         (args.output/"timecourses.csv").open("w",newline="",encoding="utf-8") as samples:
        trace_writer=None
        sample_writer=None
        for step in range(steps+1):
            positions=np.array([c.position for c in output_cells])
            assert np.isfinite(positions).all()
            energy,near_zero=pair_potential(positions,radii,adhesion,
                                           sim.mechanics.interaction_radius,sim.mechanics.repulsion_strength)
            change=0. if previous is None else energy-previous
            tolerance=1e-9*max(1.,abs(energy),abs(previous) if previous is not None else 0.)
            if previous is not None and change>tolerance:
                increases+=1
                max_increase=max(max_increase,change)
            previous=energy
            near_zero_max=max(near_zero_max,near_zero)
            drift=float(np.linalg.norm(positions.mean(axis=0)-initial_center))
            drift_max=max(drift_max,drift)
            row=dict(time=sim.time,step=step,potential=energy,potential_change=change,
                     increase_tolerance=tolerance,centroid_drift=drift,near_zero_pairs=near_zero)
            if trace_writer is None:
                trace_writer=csv.DictWriter(trace,fieldnames=list(row))
                trace_writer.writeheader()
            trace_writer.writerow(row)
            if step%every==0 or step==steps:
                trace.flush()
                np.testing.assert_array_equal([c.cell_id for c in output_cells],initial_ids)
                np.testing.assert_array_equal([c.radius for c in output_cells],radii)
                new_nodes,new_outputs=regulatory_arrays(output_cells)
                np.testing.assert_array_equal(new_nodes,nodes)
                np.testing.assert_array_equal(new_outputs,outputs)
                pressure=np.array([c.pressure for c in output_cells])
                velocity=np.array([c.velocity for c in output_cells])
                assert np.isfinite(pressure).all() and np.isfinite(velocity).all()
                observation=dict(time=sim.time,steps=step,**measure(output_cells),potential=energy,
                                 rms_radius=float(np.sqrt(np.mean(np.sum((positions-positions.mean(axis=0))**2,axis=1)))),
                                 mean_force=float(pressure.mean()),centroid_drift=drift,
                                 outside_field=int(np.any((positions<sim.bounds_min)|(positions>sim.bounds_max),axis=1).sum()))
                if sample_writer is None:
                    sample_writer=csv.DictWriter(samples,fieldnames=list(observation))
                    sample_writer.writeheader()
                sample_writer.writerow(observation)
                samples.flush()
                np.savez_compressed(args.output/f"sample_{round(step/every):04d}.npz",time=sim.time,
                                    ids=initial_ids,positions=positions,radii=radii,genes=nodes,
                                    regulatory_outputs=outputs,pressure=pressure,velocity=velocity,
                                    mechanical_strengths=np.array([[c.adhesion_strength,c.protrusion_strength] for c in output_cells]),
                                    effective_protrusion_strength=np.zeros(len(cells)))
                print(f"dt={args.dt}: {observation}",flush=True)
            if step<steps:
                pair_step(cells,sim.mechanics,args.dt)
                sim.time+=args.dt
    summary=dict(final_time=sim.time,steps=steps,cells=len(cells),
                 potential_increases_above_tolerance=increases,maximum_step_increase=max_increase,
                 maximum_centroid_drift=drift_max,maximum_near_zero_pairs=near_zero_max,
                 initial_state_sha256=meta["state_sha256"])
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")


if __name__=="__main__":
    main()
