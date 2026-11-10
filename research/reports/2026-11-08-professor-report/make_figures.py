"""Create compact report figures from archived results without running models."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).parent/'figures'; OUT.mkdir(exist_ok=True)
EXP=ROOT/'research/experiments'
plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':7})

rows=list(csv.DictReader((EXP/'2026-11-07-pair-specific-adhesion/timecourses.csv').open()))
fig,axs=plt.subplots(1,2,figsize=(6.7,2.5),constrained_layout=True)
for condition,label in [('equal','Equal'),('uniform_reduced','Uniform reduced'),('selective','Selective'),('selective_compensated','Selective compensated')]:
    selected=[r for r in rows if r['condition']==condition and float(r['dt'])==.0015625]
    times=sorted({float(r['time']) for r in selected})
    for ax,key,title in zip(axs,['heterotypic_contact_fraction','rms_radius'],['Unlike contact fraction','RMS radius']):
        a=np.array([[float(r[key]) for r in selected if float(r['time'])==t] for t in times])
        ax.plot(times,a.mean(axis=1),label=label,lw=1.5); ax.fill_between(times,a.min(axis=1),a.max(axis=1),alpha=.1)
        ax.set(xlabel='Model time',ylabel=title); ax.grid(alpha=.2)
axs[0].legend(); fig.savefig(OUT/'adhesion.png',dpi=220); plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(6.7,2.6),constrained_layout=True)
for ax,topology,label in zip(axs,['initial_sheet','contracted_E015'],['Initial sheet','Contracted network']):
    with np.load(EXP/f'2026-11-07-static-tip-selection/seed42_{topology}_E015_dt0.025.npz') as z:
        x=z['positions']; g=z['graph']
        for i,j in zip(*np.where(np.triu(g,1))): ax.plot(x[[i,j],0],x[[i,j],1],color='gray',alpha=.3,lw=.4)
        im=ax.scatter(x[:,0],x[:,1],c=z['dll4'][-1],s=22,vmin=0,vmax=1,cmap='viridis')
        ax.set(title=label,xlabel='x',ylabel='y'); ax.set_aspect('equal')
fig.colorbar(im,ax=axs,label='Abstract activity',shrink=.8); fig.savefig(OUT/'static_signaling.png',dpi=220); plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(6.7,2.6),constrained_layout=True)
for ax,timing in zip(axs,['immediate','prepatterned']):
    with np.load(EXP/f'2026-11-08-junction-sprouting/seed42_{timing}_steps32_both.npz') as z:
        x=z['positions'][-1]; g=z['bonds'][-1]; fixed=z['fixed']
        for i,j in zip(*np.where(np.triu(g,1))): ax.plot(x[[i,j],0],x[[i,j],1],color='steelblue',lw=.6)
        im=ax.scatter(x[:,0],x[:,1],c=z['dll4'][-1],s=22,vmin=0,vmax=1,cmap='viridis')
        ax.scatter(x[fixed,0],x[fixed,1],c='red',marker='x',s=12)
        ax.set(title=timing.capitalize(),xlabel='x',ylabel='y'); ax.set_aspect('equal')
fig.colorbar(im,ax=axs,label='Abstract activity',shrink=.8); fig.savefig(OUT/'combined.png',dpi=220); plt.close(fig)

rows=list(csv.DictReader((EXP/'2026-11-08-traction-budget/final_summary.csv').open()))
fig,axs=plt.subplots(1,2,figsize=(6.7,2.7),constrained_layout=True)
for b,label in [(0.,'Activity weighted'),(.25,'Follower floor'),(1.,'Uniform')]:
    for ax,key in zip(axs,['material_attached_max_x','material_attachment']):
        means=[]; low=[]; high=[]
        for F in [2.5,7.5,15.]:
            a=np.array([float(r[key]) for r in rows if r['condition']==f'F{F}_b{b}' and int(r['substeps'])==32])
            means.append(a.mean()); low.append(a.min()); high.append(a.max())
        means=np.array(means)
        ax.errorbar([2.5,7.5,15],means,yerr=[means-np.array(low),np.array(high)-means],marker='o',capsize=3,lw=1.3,label=label)
        ax.set(xlabel='Total traction budget',xticks=[2.5,7.5,15]); ax.grid(alpha=.2)
axs[0].set_ylabel('Material-connected max x'); axs[1].set_ylabel('Material attachment fraction'); axs[1].set_ylim(-.06,1.06)
axs[0].legend(); fig.savefig(OUT/'force_budget.png',dpi=220); plt.close(fig)
print('Created four report figures from archived observations')
