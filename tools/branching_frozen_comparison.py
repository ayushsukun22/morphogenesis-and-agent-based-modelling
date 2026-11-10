"""Compare evolving versus frozen GRN timestep disagreement; zero GRN errors are checks."""
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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coupled",type=Path,required=True)
    parser.add_argument("--frozen",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("choose a fresh output folder")
    records={name:read(folder/"pairwise_errors.csv") for name,folder in
             [("coupled",args.coupled),("frozen",args.frozen)]}
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    summary=[]
    for column,(coarse,fine) in enumerate([(.05,.025),(.025,.0125)]):
        selected={name:sorted([r for r in rows if float(r["coarse_dt"])==coarse and float(r["fine_dt"])==fine],
                              key=lambda r:float(r["time"])) for name,rows in records.items()}
        np.testing.assert_allclose([float(r["time"]) for r in selected["coupled"]],
                                   [float(r["time"]) for r in selected["frozen"]],rtol=0,atol=1e-8)
        row=dict(coarse_dt=coarse,fine_dt=fine)
        for name in ["coupled","frozen"]:
            data=selected[name]
            axes[0,column].plot([float(r["time"]) for r in data],
                                [float(r["position_rms_distance"]) for r in data],label=name)
            row[name+"_final_position_error"]=float(data[-1]["position_rms_distance"])
            for window,lo,hi in [("early",10,15),("late",15,20)]:
                values=[float(r["position_rms_distance"]) for r in data if lo+1e-8<float(r["time"])<=hi+1e-8]
                row[name+"_"+window+"_mean_position_error"]=float(np.mean(values))
        row["final_ratio_frozen_to_coupled"]=row["frozen_final_position_error"]/row["coupled_final_position_error"]
        summary.append(row)
        axes[0,column].set(title=f"dt {coarse} versus {fine}",xlabel="Model time (h)",ylabel="Position RMS difference (model units)")
        axes[0,column].grid(alpha=.2)
    for coarse,fine in [(.05,.025),(.025,.0125)]:
        subset=[r for r in records["frozen"] if float(r["coarse_dt"])==coarse and float(r["fine_dt"])==fine]
        for ax,key,label in [(axes[1,0],"pressure_rmse","Frozen-arm force-proxy RMSE"),
                             (axes[1,1],"polarity_vector_rmse","Frozen-arm polarity-vector RMSE")]:
            ax.plot([float(r["time"]) for r in subset],[float(r[key]) for r in subset],label=f"{coarse} vs {fine}")
            ax.set(xlabel="Model time (h)",ylabel=label)
            ax.grid(alpha=.2)
    axes[0,0].legend()
    axes[1,0].legend(fontsize=8)
    fig.suptitle("Same checkpoint and original memory rule; regulatory outputs frozen in one protocol")
    args.output.mkdir(parents=True)
    fig.savefig(args.output/"frozen_regulation_comparison.png",dpi=170)
    with (args.output/"comparison.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    paths=[folder/"pairwise_errors.csv" for folder in [args.coupled,args.frozen]]
    provenance=dict(inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                     script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/"provenance.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    for row in summary:
        print(row)


if __name__=="__main__":
    main()
