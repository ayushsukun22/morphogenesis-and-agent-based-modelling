"""Factorial VEGF-like guidance / contact inhibition assay, separate from app."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_adhesion_assay import forces_energy
from branching_shared_state import source_hashes


def initial(seed):
    rng=np.random.default_rng(seed)
    x=np.array([[i-5,j-2.5,0.] for i in range(6) for j in range(6)])
    x[:,:2]+=rng.normal(0,.025,(36,2))
    fixed=np.repeat(np.arange(6)<2,6)
    return x,fixed,np.clip(.5+rng.normal(0,.03,36),0,1)


def cue(x):
    # A uniform +x direction, no imposed lateral branch paths.
    v=1/(1+np.exp(-np.clip((x[:,0]+2)/3,-50,50)))
    grad=np.zeros_like(x); grad[:,0]=v*(1-v)/3
    return v,grad


def contact_graph(x,cutoff=1.2):
    d=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
    return (d>1e-8)&(d<=cutoff)


def components(graph,allowed=None):
    remaining=set(range(len(graph))) if allowed is None else set(np.flatnonzero(allowed))
    result=[]
    while remaining:
        todo=[remaining.pop()]; group=[]
        while todo:
            i=todo.pop(); group.append(i)
            found=set(np.flatnonzero(graph[i]))&remaining
            remaining-=found; todo.extend(found)
        result.append(group)
    return result


def regulate(x,dll4,notch,competition,dt):
    graph=contact_graph(x)
    signal=graph@dll4/np.maximum(graph.sum(axis=1),1)
    v,_=cue(x)
    # Abstract lateral inhibition: neighboring ligand suppresses own ligand.
    new_notch=notch+dt*(signal-notch)/.5
    target=v/(1+(notch/.3)**4) if competition else v
    new_dll4=dll4+dt*(target-dll4)/.5
    return new_dll4,new_notch


def traction(x,dll4,guidance):
    if guidance:
        _,direction=cue(x)
    else:
        delta=x[None,:,:]-x[:,None,:]
        dist=np.linalg.norm(delta,axis=2)
        valid=(dist>1e-8)&(dist<=2.5)
        direction=-np.sum(delta/np.maximum(dist[:,:,None],1e-8)*valid[:,:,None],axis=1)
    direction=direction/np.maximum(np.linalg.norm(direction,axis=1,keepdims=True),1e-8)
    return (.15+2.35*dll4)[:,None]*direction


def readout(x,fixed,dll4):
    graph=contact_graph(x)
    groups=components(graph)
    attached=np.zeros(len(x),bool)
    for g in groups:
        if fixed[g].any(): attached[g]=True
    candidates=[]
    for g in components(graph,attached&(x[:,0]>0)):
        if len(g)>=3 and x[g,0].max()>=2 and np.ptp(x[g,1])<=2:
            candidates.append(set(g))
    mobile=~fixed
    covariance=np.cov(x[mobile,:2].T)
    eig=np.linalg.eigvalsh(covariance)
    return dict(max_attached_extension=float(max(0,x[attached,0].max())),
                attached_mobile_fraction=float(attached[mobile].mean()),
                detached_mobile_cells=int((mobile&~attached).sum()),
                components=len(groups),mobile_mean_x=float(x[mobile,0].mean()),
                high_tip_activity_cells=int((dll4>.6).sum()),tip_activity_std=float(dll4.std()),
                narrow_extension_candidates=len(candidates),
                planar_elongation=float(np.sqrt(eig[-1]/max(eig[0],1e-12)))),candidates


def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); root=args.output
    if (root/'settings.json').exists(): raise RuntimeError('Use fresh output')
    root.mkdir(parents=True,exist_ok=True)
    settings=dict(seeds=[42,43,44],outer_dt=.05,substeps=[16,32],horizon=12.,
        conditions=dict(neither=[False,False],guidance=[True,False],competition=[False,True],both=[True,True]),
        cell_count=36,fixed_parent_cells=12,radius=.5,adhesion=.7,cutoff=2.5,repulsion=4,damping=2,
        contact_cutoff=1.2,regulatory_tau=.5,notch_half_inhibition=.3,hill_exponent=4,
        basal_traction=.15,tip_traction_increment=2.35,source_sha256=source_hashes(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256={f:hashlib.sha256((Path(__file__).parent/f).read_bytes()).hexdigest() for f in ['branching_adhesion_assay.py','branching_shared_state.py','branching_baseline.py']})
    (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    rows=[]; finals=[]; states={}; histories={}
    for seed in settings['seeds']:
        start,fixed,start_d=initial(seed)
        np.savez_compressed(root/f'initial_seed{seed}.npz',positions=start,fixed=fixed,dll4=start_d,notch=np.zeros(36))
        for substeps in settings['substeps']:
            for condition,(guidance,competition) in settings['conditions'].items():
                x=start.copy(); dll4=start_d.copy(); notch=np.zeros(36)
                trajectory=[]; ds=[]; ns=[]; times=[]; candidates=[]
                increases=0; max_increase=0.; max_overlap=0.; anchor_drift=0.
                for outer in range(241):
                    if outer%10==0:
                        r,c=readout(x,fixed,dll4)
                        row=dict(seed=seed,substeps=substeps,condition=condition,time=outer*.05,**r)
                        rows.append(row); candidates.append(c)
                        trajectory.append(x.copy()); ds.append(dll4.copy()); ns.append(notch.copy()); times.append(outer*.05)
                    if outer==240: break
                    dll4,notch=regulate(x,dll4,notch,competition,.05)
                    active=traction(x,dll4,guidance); active[fixed]=0
                    _,u=forces_energy(x,np.zeros(36),1,1); previous=u-float(np.sum(active*x))
                    for inner in range(substeps):
                        f,_=forces_energy(x,np.zeros(36),1,1); f+=active; f[fixed]=0
                        x+=.05/substeps*f/2
                        _,u=forces_energy(x,np.zeros(36),1,1); w=u-float(np.sum(active*x))
                        change=w-previous
                        if change>1e-9*max(1,abs(w),abs(previous)):
                            increases+=1; max_increase=max(max_increase,change)
                        previous=w
                        assert np.isfinite(x).all()
                    distance=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
                    max_overlap=max(max_overlap,float(np.maximum(1-distance[np.triu_indices(36,1)],0).max()))
                    anchor_drift=max(anchor_drift,float(np.max(np.abs(x[fixed]-start[fixed]))))
                # Persistence proxy: candidates share >=half the smaller ID set over 3 samples (one time unit).
                persistent_windows=0
                for a,b,c in zip(candidates,candidates[1:],candidates[2:]):
                    match=lambda s,t:len(s&t)>=.5*min(len(s),len(t))
                    persistent_windows+=int(any(match(i,j) and match(j,k) for i in a for j in b for k in c))
                finals.append(dict(row,persistent_candidate_windows=persistent_windows,
                    active_work_potential_increases=increases,max_increase=max_increase,
                    max_overlap=max_overlap,anchor_drift=anchor_drift))
                np.savez_compressed(root/f'seed{seed}_steps{substeps}_{condition}.npz',
                    positions=trajectory,dll4=ds,notch=ns,time=times,fixed=fixed)
                states[seed,substeps,condition]=(x.copy(),dll4.copy(),notch.copy())
                histories[seed,substeps,condition]=(np.array(trajectory),np.array(ds))
                print(finals[-1],flush=True)
    write(root/'timecourses.csv',rows); write(root/'final_summary.csv',finals)
    differences=[]
    for seed in settings['seeds']:
        for condition in settings['conditions']:
            a=states[seed,16,condition]; b=states[seed,32,condition]
            differences.append(dict(seed=seed,condition=condition,
                position_rms_difference=float(np.sqrt(np.mean(np.sum((a[0]-b[0])**2,axis=1)))),
                dll4_rmse=float(np.sqrt(np.mean((a[1]-b[1])**2))),notch_rmse=float(np.sqrt(np.mean((a[2]-b[2])**2)))))
    write(root/'resolution_differences.csv',differences)
    fig,axes=plt.subplots(1,4,figsize=(14,4),constrained_layout=True)
    for ax,condition in zip(axes,settings['conditions']):
        x,d=histories[42,32,condition]; start=x[0]; end=x[-1]
        for i in range(36): ax.plot(x[:,i,0],x[:,i,1],color='gray',alpha=.2,lw=.5)
        ax.scatter(start[:,0],start[:,1],s=10,facecolors='none',edgecolors='gray')
        im=ax.scatter(end[:,0],end[:,1],c=d[-1],vmin=0,vmax=1,cmap='viridis',s=32)
        ax.scatter(end[fixed,0],end[fixed,1],marker='x',c='red',s=15)
        ax.set(title=condition,xlabel='x (model units)',ylabel='y (model units)'); ax.set_aspect('equal')
    fig.colorbar(im,ax=axes,label='Abstract tip/Dll4 activity',shrink=.7)
    fig.savefig(root/'trajectories_seed42.png',dpi=170); plt.close(fig)
    assert settings['source_sha256']==source_hashes()


if __name__=='__main__': main()
