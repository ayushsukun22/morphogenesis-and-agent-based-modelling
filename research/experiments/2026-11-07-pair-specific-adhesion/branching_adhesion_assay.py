"""Controlled, passive two-identity adhesion assay; no application edits."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_shared_state import source_hashes


def forces_energy(x, labels, same, cross):
    # F_ij points toward j. All pairs in range receive adhesion, including overlaps.
    delta=x[None,:,:]-x[:,None,:]
    d=np.linalg.norm(delta,axis=2)
    valid=(d>1e-8)&(d<=2.5)
    adhesion=.7**2*np.where(labels[:,None]==labels[None,:],same,cross)
    magnitude=np.where(valid,adhesion-4*np.maximum(1-d,0),0.)
    force=np.sum(delta/np.maximum(d[:,:,None],1e-8)*magnitude[:,:,None],axis=1)
    # Continuous cutoff-shifted potential: F = -gradient(U).
    u=np.where(valid,2*np.maximum(1-d,0)**2+adhesion*(d-2.5),0.)
    return force,float(np.triu(u,1).sum())


def metrics(x, labels):
    d=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
    edges=np.triu((d>1e-8)&(d<=1.2),1)  # fixed geometric contact threshold
    mixed=labels[:,None]!=labels[None,:]
    adjacency=edges|edges.T
    remaining=set(range(len(x))); sizes=[]
    while remaining:
        todo=[remaining.pop()]; size=0
        while todo:
            i=todo.pop(); size+=1
            found=set(np.flatnonzero(adjacency[i]))&remaining
            remaining-=found; todo.extend(found)
        sizes.append(size)
    centered=x-x.mean(axis=0)
    eig=np.linalg.eigvalsh(centered.T@centered/len(x))
    return dict(rms_radius=float(np.sqrt(np.mean(np.sum(centered**2,axis=1)))),
                heterotypic_contact_fraction=float((edges&mixed).sum()/edges.sum()) if edges.any() else 0.,
                contact_edges=int(edges.sum()),components=len(sizes),largest_component_fraction=max(sizes)/len(x),
                elongation=float(np.sqrt(eig[-1]/max(eig[0],1e-12))),
                type_centroid_separation=float(np.linalg.norm(x[labels==0].mean(axis=0)-x[labels==1].mean(axis=0))))


def initial(seed):
    rng=np.random.default_rng(seed)
    x=np.array(list(np.ndindex(4,4,4)),dtype=float)*1.05
    x+=rng.normal(0,.025,x.shape); x-=x.mean(axis=0)
    labels=np.repeat([0,1],32); rng.shuffle(labels)
    return x,labels


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    root=Path('research/experiments/2026-11-07-pair-specific-adhesion')
    root.mkdir(parents=True,exist_ok=True)
    if (root/'final_summary.csv').exists():
        raise RuntimeError('Results already exist; preserve prior experiment')
    settings=dict(seeds=[42,43,44],dt=[.003125,.0015625],horizon=5.,cells=64,radius=.5,
                  per_cell_adhesion=.7,cutoff=2.5,repulsion=4.,damping=2.,contact_threshold=1.2,
                  source_sha256=source_hashes(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  protocol='synthetic mixed 4x4x4 aggregate; fixed identities, radii and adhesion; passive mechanics only; no GRN/fields/growth/polarity/propulsion')
    (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    rows=[]; finals=[]; snapshots={}
    for seed in settings['seeds']:
        start,labels=initial(seed)
        distance=np.linalg.norm(start[:,None,:]-start[None,:,:],axis=2)
        pairs=np.triu((distance>1e-8)&(distance<=2.5),1)
        mixed=labels[:,None]!=labels[None,:]
        mean_k=float(np.where(mixed,.2,1.)[pairs].mean())
        conditions=dict(equal=(1.,1.),uniform_reduced=(mean_k,mean_k),
                        selective=(1.,.2),selective_compensated=(1/mean_k,.2/mean_k))
        np.savez_compressed(root/f'initial_seed{seed}.npz',positions=start,labels=labels,
                            mean_selective_multiplier=mean_k)
        for dt in settings['dt']:
            for condition,(same,cross) in conditions.items():
                x=start.copy(); _,previous=forces_energy(x,labels,same,cross)
                increases=0; max_increase=0.; drift=0.
                for step in range(round(5/dt)+1):
                    f,energy=forces_energy(x,labels,same,cross)
                    change=energy-previous
                    if change>1e-9*max(1,abs(previous),abs(energy)):
                        increases+=1; max_increase=max(max_increase,change)
                    previous=energy
                    assert np.isfinite(x).all() and np.isfinite(f).all()
                    drift=max(drift,float(np.linalg.norm(x.mean(axis=0)-start.mean(axis=0))))
                    if step%round(.5/dt)==0:
                        row=dict(seed=seed,dt=dt,condition=condition,time=step*dt,
                                 same_multiplier=same,cross_multiplier=cross,potential=energy,**metrics(x,labels))
                        rows.append(row)
                    if step==round(5/dt):
                        finals.append(dict(row,potential_increases=increases,max_potential_increase=max_increase,max_centroid_drift=drift))
                        np.savez_compressed(root/f'final_seed{seed}_dt{dt}_{condition}.npz',positions=x,labels=labels)
                        snapshots[seed,dt,condition]=x.copy()
                        print(seed,dt,condition,finals[-1],flush=True)
                    else:
                        x+=dt*f/2
    write_csv(root/'timecourses.csv',rows); write_csv(root/'final_summary.csv',finals)
    differences=[]
    for seed in settings['seeds']:
        for condition in conditions:
            a=snapshots[seed,.003125,condition]; b=snapshots[seed,.0015625,condition]
            differences.append(dict(seed=seed,condition=condition,position_rms_difference=float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))))
    write_csv(root/'resolution_differences.csv',differences)
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    for condition in conditions:
        selected=[r for r in rows if r['condition']==condition and r['dt']==.0015625]
        times=sorted({r['time'] for r in selected})
        for ax,key in zip(axes,['heterotypic_contact_fraction','rms_radius','largest_component_fraction']):
            values=np.array([[r[key] for r in selected if r['time']==t] for t in times])
            ax.plot(times,values.mean(axis=1),label=condition)
            ax.fill_between(times,values.min(axis=1),values.max(axis=1),alpha=.12)
            ax.set(xlabel='Model time',ylabel=key.replace('_',' ')); ax.grid(alpha=.2)
    axes[0].legend(fontsize=7); fig.savefig(root/'adhesion_comparison.png',dpi=170); plt.close(fig)
    assert settings['source_sha256']==source_hashes()


if __name__=='__main__':
    main()
