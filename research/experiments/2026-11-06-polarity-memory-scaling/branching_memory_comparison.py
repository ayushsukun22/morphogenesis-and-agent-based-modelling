"""Compare E007 and E008 timestep discrepancies for the two adjacent pairs."""
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
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--scaled", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a fresh output folder")
    records = {name: read(folder/"pairwise_errors.csv") for name,folder in
               [("original",args.original),("scaled",args.scaled)]}
    fig, axes = plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    result = []
    for column,(coarse,fine) in enumerate([(.05,.025),(.025,.0125)]):
        pair = {name: sorted([r for r in rows if float(r["coarse_dt"])==coarse and float(r["fine_dt"])==fine],
                             key=lambda r:float(r["time"])) for name,rows in records.items()}
        np.testing.assert_allclose([float(r["time"]) for r in pair["original"]],
                                   [float(r["time"]) for r in pair["scaled"]],rtol=0,atol=1e-8)
        for row,(metric,label) in enumerate([("position_rms_distance","Position RMS difference (model units)"),
                                             ("grn_rmse","Six-node GRN RMSE")]):
            summary=dict(coarse_dt=coarse,fine_dt=fine,metric=metric)
            for name in ["original","scaled"]:
                observations=pair[name]
                axes[row,column].plot([float(r["time"]) for r in observations],
                                       [float(r[metric]) for r in observations],
                                       label="Fixed per-update memory" if name=="original" else "Time-scaled memory")
                summary[name+"_final"]=float(observations[-1][metric])
                for window,lo,hi in [("early",10,15),("late",15,20)]:
                    values=[float(r[metric]) for r in observations if lo+1e-8<float(r["time"])<=hi+1e-8]
                    summary[name+"_"+window+"_mean"]=float(np.mean(values))
            summary["final_ratio_scaled_to_original"]=summary["scaled_final"]/summary["original_final"]
            result.append(summary)
            axes[row,column].set(title=f"dt {coarse} versus {fine}",xlabel="Model time (h)",ylabel=label)
            axes[row,column].grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Polarity-memory interpretation: effect on timestep disagreement")
    args.output.mkdir(parents=True)
    fig.savefig(args.output/"memory_scaling_comparison.png",dpi=170)
    with (args.output/"comparison.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    paths=[folder/name for folder in [args.original,args.scaled] for name in ["pairwise_errors.csv","final_summary.csv"]]
    provenance=dict(inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                     script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/"provenance.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    for r in result:
        print(r)


if __name__=="__main__":
    main()
