"""Compare passive refinement trajectories and every-step consistency checks."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from branching_baseline import load_simulation
from branching_shared_state import decode,registry,source_hashes
from branching_passive_refinement import pair_step


def read(path):
    with path.open(newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write(path,rows):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def current_forces(frame,settings,mechanics):
    """Evaluate original-ordered passive force at current saved positions."""
    index={int(k):i for i,k in enumerate(frame["ids"])}
    cells=[SimpleNamespace(cell_id=k,position=frame["positions"][index[k]].copy(),
                           radius=float(frame["radii"][index[k]]),
                           adhesion_strength=float(frame["mechanical_strengths"][index[k],0]))
           for k in settings["force_order_ids"]]
    pair_step(cells,mechanics,0.)
    by_id={c.cell_id:c for c in cells}
    return np.array([by_id[int(k)].velocity*mechanics.damping for k in frame["ids"]])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders",type=Path,nargs="+")
    parser.add_argument("--checkpoint",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("choose a new output folder")
    metadata=json.loads((args.checkpoint/"settings.json").read_text())
    assert metadata["source_sha256"]==source_hashes()
    payload=gzip.decompress((args.checkpoint/"state.json.gz").read_bytes())
    assert hashlib.sha256(payload).hexdigest()==metadata["state_sha256"]
    model=decode(json.loads(payload),registry(load_simulation()()))
    arms=[]
    inputs={str(args.checkpoint/"state.json.gz"):hashlib.sha256((args.checkpoint/"state.json.gz").read_bytes()).hexdigest()}
    for folder in args.folders:
        settings=json.loads((folder/"settings.json").read_text())
        assert settings["initial_state_sha256"]==metadata["state_sha256"]
        assert settings["source_sha256"]==metadata["source_sha256"]
        frames=[]
        for path in sorted(folder.glob("sample_*.npz")):
            with np.load(path) as data:
                frame={k:data[k].copy() for k in data.files}
            frame["current_force"]=current_forces(frame,settings,model.mechanics)
            frames.append(frame)
            inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        for name in ["settings.json","summary.json","step_diagnostics.csv","timecourses.csv"]:
            path=folder/name
            inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        arms.append(dict(dt=settings["dt"],frames=frames,rows=read(folder/"timecourses.csv"),
                         trace=read(folder/"step_diagnostics.csv"),summary=json.loads((folder/"summary.json").read_text())))
    arms.sort(key=lambda a:a["dt"],reverse=True)
    pairs,finals,windows=[],[],[]
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for arm in arms:
        dt=arm["dt"]
        rows=arm["rows"]
        axes[0,0].plot([float(r["time"]) for r in rows],[float(r["rms_radius"]) for r in rows],label=f"dt={dt}")
        axes[0,1].plot([float(r["time"]) for r in arm["trace"]],[float(r["potential"]) for r in arm["trace"]],label=f"dt={dt}")
        last=rows[-1]
        finals.append(dict(dt=dt,rms_radius=float(last["rms_radius"]),potential=float(last["potential"]),
                           latest_step_mean_force=float(last["mean_force"]),
                           current_state_mean_force=float(np.linalg.norm(arm["frames"][-1]["current_force"],axis=1).mean()),
                           elongation=float(last["elongation"]),components=int(last["components"]),
                           largest_component_fraction=float(last["largest_component_fraction"]),
                           steps=arm["summary"]["steps"],step_potential_increases=arm["summary"]["potential_increases_above_tolerance"],
                           maximum_step_increase=arm["summary"]["maximum_step_increase"],
                           max_centroid_drift=arm["summary"]["maximum_centroid_drift"],
                           max_near_zero_pairs=arm["summary"]["maximum_near_zero_pairs"],
                           max_outside_field=max(int(r["outside_field"]) for r in rows)))
    for i,coarse in enumerate(arms[:-1]):
        for fine in arms[i+1:]:
            assert len(coarse["frames"])==len(fine["frames"])
            subset=[]
            for a,b in zip(coarse["frames"],fine["frames"]):
                assert abs(float(a["time"]-b["time"]))<1e-8
                np.testing.assert_array_equal(a["ids"],b["ids"])
                error=float(np.sqrt(np.mean(np.sum((a["positions"]-b["positions"])**2,axis=1))))
                row=dict(coarse_dt=coarse["dt"],fine_dt=fine["dt"],time=float(b["time"]),
                         position_rms_distance=error,position_error_in_mean_radii=error/float(b["radii"].mean()),
                         current_force_rms_distance=float(np.sqrt(np.mean(np.sum((a["current_force"]-b["current_force"])**2,axis=1)))))
                pairs.append(row)
                subset.append(row)
            for ax,key in [(axes[1,0],"position_rms_distance"),(axes[1,1],"current_force_rms_distance")]:
                ax.plot([r["time"] for r in subset],[r[key] for r in subset],label=f"{coarse['dt']} vs {fine['dt']}")
            summary=dict(coarse_dt=coarse["dt"],fine_dt=fine["dt"],final_position_difference=subset[-1]["position_rms_distance"],
                         final_difference_in_mean_radii=subset[-1]["position_error_in_mean_radii"],
                         final_current_force_difference=subset[-1]["current_force_rms_distance"])
            for window,lo,hi in [("early",10,15),("late",15,20)]:
                summary[window+"_mean_position_difference"]=float(np.mean([r["position_rms_distance"] for r in subset if lo+1e-8<r["time"]<=hi+1e-8]))
            windows.append(summary)
    for ax,label in zip(axes.flat,["RMS radius (model units)","Total pair potential (arbitrary units)",
                                   "Per-cell position RMS difference","Current-state force-vector RMS difference"]):
        ax.set(xlabel="Model time (h)",ylabel=label)
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    axes[1,0].legend(fontsize=7)
    fig.suptitle("Passive mechanical refinement; identical checkpoint; potential checked every update")
    args.output.mkdir(parents=True)
    fig.savefig(args.output/"refinement_comparison.png",dpi=170)
    plt.close(fig)
    for name,rows in [("final_summary.csv",finals),("pairwise_differences.csv",pairs),("window_summary.csv",windows)]:
        write(args.output/name,rows)
    helpers=["branching_baseline.py","branching_shared_state.py","branching_passive_refinement.py"]
    provenance=dict(inputs=inputs,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    helper_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in helpers})
    (args.output/"provenance.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    for row in finals+windows:
        print(row)


if __name__=="__main__":
    main()
