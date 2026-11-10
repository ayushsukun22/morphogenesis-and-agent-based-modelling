"""Compare coupled mechanical resolutions with identical outer clocks."""
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
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("choose a fresh output folder")
    arms=[]
    hashes={}
    for folder in args.folders:
        settings=json.loads((folder/"settings.json").read_text())
        summary=json.loads((folder/"summary.json").read_text())
        frames=[]
        for path in sorted(folder.glob("sample_*.npz")):
            with np.load(path) as data:
                frames.append({k:data[k].copy() for k in data.files})
            hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        for name in ["settings.json","summary.json","timecourses.csv"]:
            path=folder/name
            hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        arms.append(dict(n=settings["substeps"],h=settings["inner_dt"],settings=settings,
                         summary=summary,frames=frames,rows=read(folder/"timecourses.csv")))
    arms.sort(key=lambda a:a["n"])
    assert len({a["settings"]["initial_state_sha256"] for a in arms})==1
    assert len({a["settings"]["outer_dt"] for a in arms})==1
    initial=arms[0]["frames"][0]
    assert all(np.array_equal(initial[k],arm["frames"][0][k]) for arm in arms for k in initial)
    finals,pairs,windows=[],[],[]
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for arm in arms:
        rows=arm["rows"]
        time=[float(r["time"]) for r in rows]
        for ax,key,label in [(axes[0,0],"rms_radius","RMS radius (model units)"),
                             (axes[0,1],"high_snail_cells","High-SNAIL cells"),
                             (axes[0,2],"mean_yap","Mean YAP"),
                             (axes[1,0],"mean_grn_pressure_input","Mean GRN force-proxy input (symlog)")]:
            ax.plot(time,[float(r[key]) for r in rows],label=f"{arm['n']} steps, h={arm['h']}")
            ax.set(xlabel="Model time (h)",ylabel=label)
            ax.grid(alpha=.2)
        last=rows[-1]
        p_input=arm["frames"][-1]["grn_pressure_inputs"]
        finals.append(dict(substeps=arm["n"],inner_dt=arm["h"],
                           rms_radius=float(last["rms_radius"]),high_snail=int(last["high_snail_cells"]),
                           mean_snail=float(last["mean_snail"]),mean_yap=float(last["mean_yap"]),
                           mean_grn_pressure_input=float(last["mean_grn_pressure_input"]),
                           fraction_grn_inputs_above_one=float(np.mean(p_input>1.)),
                           mean_yap_squeeze_input=float(np.minimum(np.maximum(p_input-1.,0.),1.).mean()),
                           mean_pressure_buffer=float(last["mean_pressure_buffer"]),
                           mean_last_inner_pressure=float(last["mean_last_inner_pressure"]),
                           peak_high_snail=max(int(r["high_snail_cells"]) for r in rows),
                           max_mean_grn_pressure_input=max(float(r["mean_grn_pressure_input"]) for r in rows),
                           components=int(last["components"]),largest_component_fraction=float(last["largest_component_fraction"]),
                           elongation=float(last["elongation"]),
                           max_outside_field=max(int(r["outside_field"]) for r in rows),
                           outer_steps=arm["summary"]["outer_steps"],inner_steps=arm["summary"]["inner_steps"]))
    for i,coarse in enumerate(arms[:-1]):
        for fine in arms[i+1:]:
            subset=[]
            assert len(coarse["frames"])==len(fine["frames"])
            for a,b in zip(coarse["frames"],fine["frames"]):
                assert float(a["time"])==float(b["time"])
                np.testing.assert_array_equal(a["ids"],b["ids"])
                error=float(np.sqrt(np.mean(np.sum((a["positions"]-b["positions"])**2,axis=1))))
                r=dict(coarse_substeps=coarse["n"],fine_substeps=fine["n"],time=float(b["time"]),
                       position_rms_distance=error,position_error_in_mean_radii=error/float(b["radii"].mean()),
                       grn_rmse=float(np.sqrt(np.mean((a["genes"]-b["genes"])**2))),
                       grn_pressure_input_rmse=float(np.sqrt(np.mean((a["grn_pressure_inputs"]-b["grn_pressure_inputs"])**2))),
                       pressure_buffer_rmse=float(np.sqrt(np.mean((a["pressure"]-b["pressure"])**2))),
                       neighbor_signal_rmse=float(np.sqrt(np.mean((a["neighbor_signal"]-b["neighbor_signal"])**2))),
                       oxygen_rmse=float(np.sqrt(np.mean((a["oxygen_grid"]-b["oxygen_grid"])**2))),
                       morphogen_rmse=float(np.sqrt(np.mean((a["morphogen_grid"]-b["morphogen_grid"])**2))))
                pairs.append(r)
                subset.append(r)
            for ax,key,label in [(axes[1,1],"position_rms_distance","Per-cell position RMS difference"),
                                 (axes[1,2],"grn_rmse","Six-node GRN RMSE")]:
                ax.plot([r["time"] for r in subset],[r[key] for r in subset],label=f"{coarse['n']} vs {fine['n']} inner steps")
                ax.set(xlabel="Model time (h)",ylabel=label)
                ax.grid(alpha=.2)
            summary=dict(coarse_substeps=coarse["n"],fine_substeps=fine["n"],
                         final_position_difference=subset[-1]["position_rms_distance"],
                         final_difference_in_mean_radii=subset[-1]["position_error_in_mean_radii"],
                         final_grn_rmse=subset[-1]["grn_rmse"],
                         final_grn_pressure_input_rmse=subset[-1]["grn_pressure_input_rmse"])
            for window,lo,hi in [("early",10,15),("late",15,20)]:
                selected=[r for r in subset if lo+1e-8<r["time"]<=hi+1e-8]
                for metric in ["position_rms_distance","grn_rmse"]:
                    summary[window+"_mean_"+metric]=float(np.mean([r[metric] for r in selected]))
            windows.append(summary)
    axes[1,0].set_yscale("symlog",linthresh=1.)
    axes[0,0].legend(fontsize=7)
    axes[1,1].legend(fontsize=7)
    fig.suptitle("Coupled simulation: fixed outer dt=0.05; original outer-lag feedback; held inner polarity")
    args.output.mkdir(parents=True)
    fig.savefig(args.output/"coupled_substeps.png",dpi=170)
    plt.close(fig)
    for name,rows in [("final_summary.csv",finals),("pairwise_differences.csv",pairs),("window_summary.csv",windows)]:
        write(args.output/name,rows)
    (args.output/"provenance.json").write_text(json.dumps(dict(inputs=hashes,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding="utf-8")
    for row in finals+windows:
        print(row)


if __name__=="__main__":
    main()
