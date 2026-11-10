"""Describe pair-force-only runs, including their model potential and centroid drift."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from branching_baseline import load_simulation
from branching_shared_state import decode,registry,source_hashes


def pair_potential(positions,radii,adhesion,interaction_radius,repulsion_strength):
    """Potential whose radial derivative matches the existing cutoff pair force.

    Zero outside cutoff; continuous at cutoff. The near-zero force skip is
    counted separately rather than interpreted as a well-defined gradient.
    """
    i,j=np.triu_indices(len(positions),1)
    distance=np.linalg.norm(positions[j]-positions[i],axis=1)
    contact=radii[i]+radii[j]
    active=distance<=interaction_radius
    overlap=np.maximum(contact-distance,0.)
    cutoff_overlap=np.maximum(contact-interaction_radius,0.)
    potential=(adhesion[i]*adhesion[j]*(distance-interaction_radius)
               +.5*repulsion_strength*(overlap**2-cutoff_overlap**2))
    return float(potential[active].sum()),int(np.sum(distance<1e-8))


def read(path):
    with path.open(newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write(path,rows):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders",type=Path,nargs="+")
    parser.add_argument("--checkpoint",type=Path,required=True)
    parser.add_argument("--reference-analysis",type=Path,required=True)
    parser.add_argument("--passive-analysis",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("choose a fresh output folder")
    meta=json.loads((args.checkpoint/"settings.json").read_text())
    if meta["source_sha256"]!=source_hashes():
        raise ValueError("Application source differs from checkpoint provenance")
    payload=gzip.decompress((args.checkpoint/"state.json.gz").read_bytes())
    assert hashlib.sha256(payload).hexdigest()==meta["state_sha256"]
    model=decode(json.loads(payload),registry(load_simulation()()))
    radius=model.mechanics.interaction_radius
    repulsion=model.mechanics.repulsion_strength
    times,summaries,contrasts=[],[],[]
    hashes={str(args.checkpoint/"state.json.gz"):hashlib.sha256((args.checkpoint/"state.json.gz").read_bytes()).hexdigest()}
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for folder in args.folders:
        settings=json.loads((folder/"settings.json").read_text())
        assert settings["freeze_regulation"] and settings["protrusion_off"]
        dt=settings["dt"]
        rows=[]
        anchor=None
        for path in sorted(folder.glob("sample_*.npz")):
            with np.load(path) as frame:
                position=frame["positions"]
                center=position.mean(axis=0)
                if anchor is None:
                    anchor=center.copy()
                energy,near_zero=pair_potential(position,frame["radii"],frame["mechanical_strengths"][:,0],radius,repulsion)
                assert np.all(frame["effective_protrusion_strength"]==0)
                rows.append(dict(dt=dt,time=float(frame["time"]),potential=energy,
                                 centroid_drift=float(np.linalg.norm(center-anchor)),near_zero_pairs=near_zero))
            hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        times.extend(rows)
        delta=np.diff([r["potential"] for r in rows])
        tolerance=1e-9*max(1.,max(abs(r["potential"]) for r in rows))
        summaries.append(dict(dt=dt,initial_potential=rows[0]["potential"],final_potential=rows[-1]["potential"],
                              maximum_saved_increase=float(max(0.,delta.max())),
                              saved_increases_above_tolerance=int(np.sum(delta>tolerance)),
                              increase_tolerance=tolerance,
                              maximum_centroid_drift=max(r["centroid_drift"] for r in rows),
                              maximum_near_zero_pairs=max(r["near_zero_pairs"] for r in rows)))
        axes[1,0].plot([r["time"] for r in rows],[r["potential"] for r in rows],label=f"dt={dt}")
        measurements=read(folder/"timecourses.csv")
        axes[1,1].plot([float(r["time"]) for r in measurements],[float(r["rms_radius"]) for r in measurements],label=f"dt={dt}")
        hashes[str(folder/"timecourses.csv")]=hashlib.sha256((folder/"timecourses.csv").read_bytes()).hexdigest()
    errors={name:read(folder/"pairwise_errors.csv") for name,folder in
            [("reference",args.reference_analysis),("passive",args.passive_analysis)]}
    for name,folder in [("reference",args.reference_analysis),("passive",args.passive_analysis)]:
        path=folder/"pairwise_errors.csv"
        hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    for column,(coarse,fine) in enumerate([(.05,.025),(.025,.0125)]):
        summary=dict(coarse_dt=coarse,fine_dt=fine)
        for name in ["reference","passive"]:
            subset=sorted([r for r in errors[name] if float(r["coarse_dt"])==coarse and float(r["fine_dt"])==fine],key=lambda r:float(r["time"]))
            axes[0,column].plot([float(r["time"]) for r in subset],[float(r["position_rms_distance"]) for r in subset],
                                label="Frozen GRN, protrusion on" if name=="reference" else "Frozen GRN, protrusion off")
            summary[name+"_final_position_difference"]=float(subset[-1]["position_rms_distance"])
            for window,lo,hi in [("early",10,15),("late",15,20)]:
                summary[name+"_"+window+"_mean"]=float(np.mean([float(r["position_rms_distance"]) for r in subset if lo+1e-8<float(r["time"])<=hi+1e-8]))
        contrasts.append(summary)
        axes[0,column].set(title=f"dt {coarse} versus {fine}",xlabel="Model time (h)",ylabel="Position RMS difference (model units)")
    axes[1,0].set(xlabel="Model time (h)",ylabel="Total model pair potential (arbitrary units)")
    axes[1,1].set(xlabel="Model time (h)",ylabel="RMS radius (model units)")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    axes[1,0].legend(fontsize=8)
    fig.suptitle("Frozen regulation, adhesion/repulsion only; potential is a mathematical diagnostic")
    args.output.mkdir(parents=True)
    fig.savefig(args.output/"passive_mechanics.png",dpi=170)
    plt.close(fig)
    write(args.output/"potential_timecourses.csv",times)
    write(args.output/"invariants_and_potential.csv",summaries)
    write(args.output/"comparison.csv",contrasts)
    provenance=dict(inputs=hashes,interaction_radius=radius,repulsion_strength=repulsion,
                     helper_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                                    for name in ["branching_baseline.py","branching_shared_state.py"]},
                     script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/"provenance.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    for row in summaries+contrasts:
        print(row)


if __name__=="__main__":
    main()
