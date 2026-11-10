"""Coupled diagnostic with fixed outer timing and inner mechanical force steps."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from branching_baseline import load_simulation,install_growth_arrest,measure
from branching_shared_state import decode,registry,canonical_state,source_hashes,MidpointRandom,regulatory_arrays


def orient_and_signal(cells,mechanics):
    """Original neighbor/polarity phase, evaluated once per outer update."""
    active=[]
    for cell in cells:
        neighbors=[]
        neighbor_snail=[]
        for other in cells:
            if other.cell_id==cell.cell_id:
                continue
            delta=other.position-cell.position
            distance=np.linalg.norm(delta)
            if distance>1e-8 and distance<=mechanics.interaction_radius:
                neighbors.append(delta/distance)
                neighbor_snail.append(other.grn.snail)
        cell.grn.neighbor_signal=float(np.mean(neighbor_snail)) if neighbor_snail else 0.
        if not neighbors:
            cell.pressure=0.
            active.append(False)
            continue
        adhesion_direction=np.sum(neighbors,axis=0)
        norm=np.linalg.norm(adhesion_direction)
        if norm<1e-8:
            cell.pressure=0.
            active.append(False)
            continue
        adhesion_direction/=norm
        proposed_polarity=-adhesion_direction*cell.grn.polarity_bias
        memory=np.clip(mechanics.polarity_memory,0.,1.)
        cell.polarity=memory*cell.polarity+(1.-memory)*proposed_polarity
        polarity_norm=np.linalg.norm(cell.polarity)
        if polarity_norm>1e-8:
            cell.polarity/=polarity_norm
        active.append(True)
    return active


def force_step(cells,mechanics,dt,active):
    """Original ordered pair forces plus the held outer-step protrusion."""
    forces={c.cell_id:np.zeros(3,dtype=np.float64) for c in cells}
    for i in range(len(cells)):
        a=cells[i]
        for j in range(i+1,len(cells)):
            b=cells[j]
            delta=b.position-a.position
            distance=np.linalg.norm(delta)
            if distance<1e-8:
                continue
            direction=delta/distance
            contact=a.radius+b.radius
            if distance>mechanics.interaction_radius:
                continue
            overlap=contact-distance
            if overlap>0.:
                repulsion=mechanics.repulsion_strength*overlap
                forces[a.cell_id]-=direction*repulsion
                forces[b.cell_id]+=direction*repulsion
            adhesion=(a.adhesion_strength*b.adhesion_strength
                      *min(1.,mechanics.interaction_radius/max(distance,1e-8)))
            forces[a.cell_id]+=direction*adhesion
            forces[b.cell_id]-=direction*adhesion
    for cell,eligible in zip(cells,active):
        if eligible:
            forces[cell.cell_id]+=cell.protrusion_strength*cell.polarity
    for cell in cells:
        force=forces[cell.cell_id]
        cell.pressure=float(np.linalg.norm(force))
        cell.velocity=force/max(mechanics.damping,1e-8)
        cell.position+=cell.velocity*dt


def split_update(cells,mechanics,dt,substeps):
    if len(cells)<2:
        # Preserve the original mechanics' early return.
        pressure=[c.pressure for c in cells]
        return pressure,pressure,[False]*len(cells)
    active=orient_and_signal(cells,mechanics)
    first=None
    for _ in range(substeps):
        force_step(cells,mechanics,dt/substeps,active)
        if first is None:
            first=[c.pressure for c in cells]
    last=[c.pressure for c in cells]
    for cell,pressure in zip(cells,first):
        cell.pressure=pressure  # original outer-start lag for next GRN update
    return first,last,active


def install_substeps(sim,substeps):
    if not isinstance(substeps,int) or substeps<1:
        raise ValueError("substeps must be a positive integer")
    stats=dict(outer_calls=0,inner_calls=0,orientation_cell_evaluations=0,
               last_pressure={c.cell_id:c.pressure for c in sim.cells.values()},active={})
    def update(cells,dt):
        _,last,active=split_update(cells,sim.mechanics,dt,substeps)
        stats["outer_calls"]+=1
        stats["inner_calls"]+=substeps
        stats["orientation_cell_evaluations"]+=len(cells)
        stats["last_pressure"]={c.cell_id:p for c,p in zip(cells,last)}
        stats["active"]={c.cell_id:int(a) for c,a in zip(cells,active)}
    sim.mechanics.update=update
    return stats


def instrument_outer_updates(sim):
    counts=dict(grn_calls=0,oxygen_steps=0,morphogen_steps=0)
    for cell in sim.cells.values():
        original=cell.grn.update
        def grn_update(dt,rng,_original=original):
            counts["grn_calls"]+=1
            return _original(dt,rng)
        cell.grn.update=grn_update
    for key,field in [("oxygen_steps",sim.oxygen_field),("morphogen_steps",sim.morphogen_field)]:
        original=field.step
        def field_update(dt,*args,_original=original,_key=key,**kwargs):
            counts[_key]+=1
            return _original(dt,*args,**kwargs)
        field.step=field_update
    return counts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint",type=Path,required=True)
    parser.add_argument("--outer-dt",type=float,default=.05)
    parser.add_argument("--substeps",type=int,required=True)
    parser.add_argument("--horizon",type=float,default=20.)
    parser.add_argument("--sample-every",type=float,default=.5)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("choose a fresh output folder")
    meta=json.loads((args.checkpoint/"settings.json").read_text())
    duration=args.horizon-meta["onset"]
    if args.substeps<1 or not all(np.isfinite(x) and x>0 for x in [args.outer_dt,args.sample_every,duration]):
        parser.error("positive finite time settings and at least one substep are required")
    if not all(np.isclose(x/args.outer_dt,round(x/args.outer_dt),rtol=0,atol=1e-8) for x in [duration,args.sample_every]):
        parser.error("duration and sample interval must align with outer-dt")
    assert source_hashes()==meta["source_sha256"]
    payload=gzip.decompress((args.checkpoint/"state.json.gz").read_bytes())
    assert hashlib.sha256(payload).hexdigest()==meta["state_sha256"]
    classes=registry(load_simulation()())
    sim=decode(json.loads(payload),classes)
    assert canonical_state(sim,classes)==payload
    ids=sorted(sim.cells)
    radii=np.array([sim.cells[k].radius for k in ids])
    install_growth_arrest(sim,meta["onset"])
    sim.rng=MidpointRandom()
    stats=install_substeps(sim,args.substeps)
    counts=instrument_outer_updates(sim)
    steps=round(duration/args.outer_dt)
    every=round(args.sample_every/args.outer_dt)
    args.output.mkdir(parents=True)
    helpers=["branching_baseline.py","branching_shared_state.py"]
    settings=dict(outer_dt=args.outer_dt,substeps=args.substeps,inner_dt=args.outer_dt/args.substeps,
                  horizon=args.horizon,sample_every=args.sample_every,initial_state_sha256=meta["state_sha256"],
                  source_sha256=source_hashes(),python=platform.python_version(),numpy=np.__version__,
                  polarity_policy="once at outer start; original memory; held direction and eligibility within outer interval",
                  neighbor_policy="sample at outer start for next GRN update",
                  pressure_policy="cell.pressure is outer-start force norm for next GRN; last-inner pressure stored separately",
                  biological_policy="active GRN/fields at outer dt; growth/division arrested; GRN noise zero",
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  helper_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in helpers})
    (args.output/"settings.json").write_text(json.dumps(settings,indent=2),encoding="utf-8")
    with (args.output/"timecourses.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=None
        for step in range(steps+1):
            if step%every==0 or step==steps:
                assert sorted(sim.cells)==ids
                cells=[sim.cells[k] for k in ids]
                np.testing.assert_array_equal([c.radius for c in cells],radii)
                genes,outputs=regulatory_arrays(cells)
                positions=np.array([c.position for c in cells])
                pressure=np.array([c.pressure for c in cells])
                last_pressure=np.array([stats["last_pressure"][k] for k in ids])
                input_pressure=np.array([c.grn.pressure for c in cells])
                arrays=dict(ids=np.array(ids),positions=positions,radii=radii,genes=genes,
                            polarity=np.array([c.polarity for c in cells]),pressure=pressure,
                            regulatory_outputs=outputs,velocity=np.array([c.velocity for c in cells]),
                            last_inner_pressure=last_pressure,grn_pressure_inputs=input_pressure,
                            neighbor_signal=np.array([c.grn.neighbor_signal for c in cells]),
                            propulsion_active=np.array([stats["active"].get(k,-1) for k in ids]),
                            oxygen_grid=sim.oxygen_field.grid,morphogen_grid=sim.morphogen_field.grid)
                assert all(np.isfinite(a).all() for a in arrays.values())
                row=dict(time=sim.time,outer_steps=step,**measure(cells),
                         rms_radius=float(np.sqrt(np.mean(np.sum((positions-positions.mean(axis=0))**2,axis=1)))),
                         mean_snail=float(genes[:,1].mean()),mean_yap=float(genes[:,2].mean()),
                         mean_pressure_buffer=float(pressure.mean()),mean_grn_pressure_input=float(input_pressure.mean()),
                         mean_last_inner_pressure=float(last_pressure.mean()),
                         outside_field=int(np.any((positions<sim.bounds_min)|(positions>sim.bounds_max),axis=1).sum()))
                if writer is None:
                    writer=csv.DictWriter(stream,fieldnames=list(row))
                    writer.writeheader()
                writer.writerow(row)
                stream.flush()
                np.savez_compressed(args.output/f"sample_{round(step/every):04d}.npz",time=sim.time,**arrays)
                print(f"substeps={args.substeps}: {row}",flush=True)
            if step<steps:
                sim.update(args.outer_dt)
    assert sim.growth_arrest_steps==steps
    assert counts==dict(grn_calls=steps*len(ids),oxygen_steps=steps,morphogen_steps=steps)
    assert sim.rng.calls==steps*len(ids)
    assert stats["outer_calls"]==steps and stats["inner_calls"]==steps*args.substeps
    assert stats["orientation_cell_evaluations"]==steps*len(ids)
    summary=dict(final_time=sim.time,cells=len(ids),outer_steps=steps,inner_steps=stats["inner_calls"],
                 orientation_cell_evaluations=stats["orientation_cell_evaluations"],
                 grn_calls=counts["grn_calls"],oxygen_steps=counts["oxygen_steps"],
                 morphogen_steps=counts["morphogen_steps"],midpoint_calls=sim.rng.calls,
                 growth_arrest_steps=sim.growth_arrest_steps,initial_state_sha256=meta["state_sha256"])
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")


if __name__=="__main__":
    main()
