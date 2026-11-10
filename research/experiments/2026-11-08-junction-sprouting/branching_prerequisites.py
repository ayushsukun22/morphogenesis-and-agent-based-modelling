"""Independent static-signaling and loaded-tissue prerequisite experiments."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_sprouting_assay import initial,contact_graph,components,write
from branching_shared_state import source_hashes


def signaling_step(d,n,graph,v,dt,hill=4,half=.3,inhibition=True):
    signal=graph@d/np.maximum(graph.sum(axis=1),1)
    target=v/(1+(n/half)**hill) if inhibition else np.full_like(d,v)
    return d+dt*(target-d)/.5,n+dt*(signal-n)/.5


def signal_metrics(d,graph):
    high=d>.6; low=d<.2
    edges=np.triu(graph,1)
    isolated=high&~np.any(graph&high[None,:],axis=1)
    return dict(activity_std=float(d.std()),high_cells=int(high.sum()),low_cells=int(low.sum()),
                isolated_high_cells=int(isolated.sum()),
                mean_edge_contrast=float(np.abs(d[:,None]-d[None,:])[edges].mean()))


def mechanical_force_energy(x,repulsion,cutoff,bonds=None,rest=None,spring=4.):
    delta=x[None,:,:]-x[:,None,:]; dist=np.linalg.norm(delta,axis=2)
    valid=dist>1e-8
    overlap=np.maximum(1-dist,0)
    if bonds is None:
        active=valid&(dist<=cutoff)
        magnitude=np.where(active,.49-repulsion*overlap,0.)
        u=np.where(active,.5*repulsion*overlap**2+.49*(dist-cutoff),0.)
    else:
        # Repulsion for every overlapping pair; initial material junctions in tension/compression.
        extension=dist-rest
        magnitude=np.where(valid,-repulsion*overlap,0.)+spring*extension*bonds
        u=np.where(valid,.5*repulsion*overlap**2,0.)+.5*spring*extension**2*bonds
    force=np.sum(delta/np.maximum(dist[:,:,None],1e-8)*magnitude[:,:,None],axis=1)
    return force,float(np.triu(u,1).sum())


def attachment(x,fixed,cutoff):
    attached=np.zeros(len(x),bool)
    for group in components(contact_graph(x,cutoff)):
        if fixed[group].any(): attached[group]=True
    return float(attached[~fixed].mean())


def hashes():
    files=['branching_prerequisites.py','branching_sprouting_assay.py','branching_adhesion_assay.py','branching_shared_state.py','branching_baseline.py']
    return {f:hashlib.sha256((Path(__file__).parent/f).read_bytes()).hexdigest() for f in files}


def signaling(root):
    root.mkdir(parents=True,exist_ok=True)
    if (root/'settings.json').exists(): raise RuntimeError('Use fresh output')
    settings=dict(seeds=[42,43,44],dt=[.05,.025],horizon=20,uniform_cue=.7,tau=.5,
        conditions=dict(no_inhibition=[4,.3,False],E015=[4,.3,True],steeper=[8,.2,True]),
        graphs=['initial_sheet','contracted_E015'],source_sha256=source_hashes(),script_sha256=hashes())
    (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    finals=[]; rows=[]; input_hashes={}; histories={}
    for seed in settings['seeds']:
        start,_,d0=initial(seed)
        f=Path('research/experiments/2026-11-07-guidance-tip-competition')/f'seed{seed}_steps32_neither.npz'
        input_hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
        with np.load(f) as z: contracted=z['positions'][-1].copy()
        for topology,x in [('initial_sheet',start),('contracted_E015',contracted)]:
            graph=contact_graph(x)
            for dt in settings['dt']:
                for name,(hill,half,inhibit) in settings['conditions'].items():
                    d=d0.copy(); n=np.zeros(36); ds=[]; ns=[]; times=[]
                    for step in range(round(20/dt)+1):
                        if step%round(.5/dt)==0:
                            rows.append(dict(seed=seed,topology=topology,dt=dt,condition=name,time=step*dt,**signal_metrics(d,graph)))
                            ds.append(d.copy()); ns.append(n.copy()); times.append(step*dt)
                        if step<round(20/dt): d,n=signaling_step(d,n,graph,.7,dt,hill,half,inhibit)
                    assert np.isfinite(ds).all() and np.isfinite(ns).all()
                    history=np.array(ds); stable=(history[np.array(times)>=15]>.6).mean(axis=0)>=.8
                    isolated_stable=stable&~np.any(graph&stable[None,:],axis=1)
                    finals.append(dict(rows[-1],stable_high_cells=int(stable.sum()),stable_isolated_high_cells=int(isolated_stable.sum())))
                    np.savez_compressed(root/f'seed{seed}_{topology}_{name}_dt{dt}.npz',positions=x,graph=graph,
                        time=times,dll4=ds,notch=ns,initial_dll4=d0)
                    histories[seed,topology,name,dt]=history
    write(root/'timecourses.csv',rows); write(root/'final_summary.csv',finals)
    diffs=[]
    for seed in settings['seeds']:
        for topology in settings['graphs']:
            for name in settings['conditions']:
                a=histories[seed,topology,name,.05]; b=histories[seed,topology,name,.025]
                diffs.append(dict(seed=seed,topology=topology,condition=name,final_rmse=float(np.sqrt(np.mean((a[-1]-b[-1])**2))),trajectory_rmse=float(np.sqrt(np.mean((a-b)**2)))))
    write(root/'resolution_differences.csv',diffs)
    (root/'input_hashes.json').write_text(json.dumps(input_hashes,indent=2)+'\n')
    fig,axes=plt.subplots(2,3,figsize=(10,6),constrained_layout=True)
    for row,topology in enumerate(settings['graphs']):
        for col,name in enumerate(settings['conditions']):
            f=root/f'seed42_{topology}_{name}_dt0.025.npz'
            with np.load(f) as z:
                x=z['positions']; graph=z['graph']; d=z['dll4'][-1]
                for i,j in zip(*np.where(np.triu(graph,1))): axes[row,col].plot(x[[i,j],0],x[[i,j],1],color='gray',lw=.3,alpha=.3)
                im=axes[row,col].scatter(x[:,0],x[:,1],c=d,vmin=0,vmax=1,cmap='viridis',s=45)
                axes[row,col].set(title=topology+' / '+name,aspect='equal')
    fig.colorbar(im,ax=axes,label='Abstract tip activity',shrink=.6); fig.savefig(root/'signaling_states.png',dpi=170); plt.close(fig)


def mechanics(root):
    root.mkdir(parents=True,exist_ok=True)
    if (root/'settings.json').exists(): raise RuntimeError('Use fresh output')
    variants=dict(original=[4.,2.5,False],stiffer=[40.,2.5,False],short_range=[40.,1.2,False],material_junctions=[40.,1.2,True])
    settings=dict(seeds=[42,43,44],dt=[.005,.0025],horizon=5.,variants=variants,loads=['none','one_front_cell','all_mobile'],
        traction_per_loaded_cell=2.5,junction_spring=4.,junction_break_distance=1.6,contact_threshold=1.2,
        source_sha256=source_hashes(),script_sha256=hashes())
    (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    finals=[]; rows=[]; states={}
    for seed in settings['seeds']:
        start,fixed,_=initial(seed)
        distance=np.linalg.norm(start[:,None,:]-start[None,:,:],axis=2)
        initial_bonds=contact_graph(start)
        # Reproducible imposed load site: frontmost column, nearest y=0 (tie by ID).
        front=np.arange(30,36); leader=int(front[np.argmin(np.abs(start[front,1]))])
        np.savez_compressed(root/f'initial_seed{seed}.npz',positions=start,fixed=fixed,bonds=initial_bonds,rest=distance,loaded_cell=leader)
        for dt in settings['dt']:
            for name,(repulsion,cutoff,material) in variants.items():
                for load in settings['loads']:
                    x=start.copy(); bonds=initial_bonds.copy() if material else None
                    active=np.zeros_like(x)
                    if load=='one_front_cell': active[leader,0]=2.5
                    if load=='all_mobile': active[~fixed,0]=2.5
                    increases=0; max_overlap=0.; broken=0; min_attach=1.; snapshots=[]; times=[]
                    _,u=mechanical_force_energy(x,repulsion,cutoff,bonds,distance)
                    previous=u-float(np.sum(active*x))
                    for step in range(round(5/dt)+1):
                        dist=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
                        max_overlap=max(max_overlap,float(np.maximum(1-dist[np.triu_indices(36,1)],0).max()))
                        if step%round(.25/dt)==0:
                            attach=attachment(x,fixed,1.2); min_attach=min(min_attach,attach)
                            row=dict(seed=seed,dt=dt,variant=name,load=load,time=step*dt,
                                attached_fraction=attach,attached_fraction_1_6=attachment(x,fixed,1.6),
                                front_displacement=float(x[leader,0]-start[leader,0]),
                                mean_mobile_displacement=float((x[~fixed,0]-start[~fixed,0]).mean()),
                                max_overlap=float(np.maximum(1-dist[np.triu_indices(36,1)],0).max()),
                                broken_initial_junctions=broken)
                            rows.append(row); snapshots.append(x.copy()); times.append(step*dt)
                        if step==round(5/dt): break
                        f,_=mechanical_force_energy(x,repulsion,cutoff,bonds,distance)
                        f+=active; f[fixed]=0; x+=dt*f/2
                        _,u=mechanical_force_energy(x,repulsion,cutoff,bonds,distance)
                        w=u-float(np.sum(active*x))
                        if w-previous>1e-9*max(1,abs(w),abs(previous)): increases+=1
                        # Break after force step; removing strained springs cannot add stored energy.
                        if material:
                            new_dist=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
                            remove=bonds&(new_dist>1.6); broken+=int(np.triu(remove,1).sum()); bonds[remove]=False
                            _,u=mechanical_force_energy(x,repulsion,cutoff,bonds,distance); w=u-float(np.sum(active*x))
                        previous=w
                        assert np.isfinite(x).all()
                        np.testing.assert_array_equal(x[fixed],start[fixed])
                    finals.append(dict(row,max_overlap_over_run=max_overlap,min_sampled_attachment=min_attach,
                        active_work_increases=increases))
                    np.savez_compressed(root/f'seed{seed}_dt{dt}_{name}_{load}.npz',positions=snapshots,time=times,
                        fixed=fixed,active=active,final_bonds=bonds if material else np.zeros((36,36),bool))
                    states[seed,dt,name,load]=x.copy()
            print('mechanics seed',seed,'dt',dt,'completed',flush=True)
    write(root/'timecourses.csv',rows); write(root/'final_summary.csv',finals)
    diffs=[]
    for seed in settings['seeds']:
        for name in variants:
            for load in settings['loads']:
                a=states[seed,.005,name,load]; b=states[seed,.0025,name,load]
                diffs.append(dict(seed=seed,variant=name,load=load,position_rms_difference=float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))))
    write(root/'resolution_differences.csv',diffs)
    fig,axes=plt.subplots(3,4,figsize=(12,8),constrained_layout=True)
    for row,load in enumerate(settings['loads']):
        for col,name in enumerate(variants):
            with np.load(root/f'seed42_dt0.0025_{name}_{load}.npz') as z:
                x=z['positions'][-1]; fixed=z['fixed']
                axes[row,col].scatter(x[~fixed,0],x[~fixed,1],s=20)
                axes[row,col].scatter(x[fixed,0],x[fixed,1],s=15,c='red',marker='x')
                axes[row,col].set(title=name+' / '+load,aspect='equal')
    fig.savefig(root/'mechanics_states.png',dpi=170); plt.close(fig)


def main():
    p=argparse.ArgumentParser(); p.add_argument('mode',choices=['signaling','mechanics']); p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(); before=source_hashes()
    (signaling if args.mode=='signaling' else mechanics)(args.output)
    assert source_hashes()==before


if __name__=='__main__': main()
