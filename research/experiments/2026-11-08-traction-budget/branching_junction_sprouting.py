"""Recombine reference signaling, guidance and material-junction mechanics."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_sprouting_assay import initial,regulate,traction,readout,contact_graph,components,write
from branching_prerequisites import mechanical_force_energy
from branching_shared_state import source_hashes


def parent_attachment(graph,fixed):
    attached=np.zeros(len(graph),bool)
    for group in components(graph):
        if fixed[group].any(): attached[group]=True
    return attached


def mechanical_interval(x,fixed,bonds,rest,active,substeps):
    """Held traction, all-pair repulsion and breakable initial springs."""
    x=x.copy(); bonds=bonds.copy(); active=active.copy(); active[fixed]=0
    _,u=mechanical_force_energy(x,40,1.2,bonds,rest)
    previous=u-float(np.sum(active*x)); broken=0; increases=0; max_overlap=0.
    min_attached=1.
    for inner in range(substeps):
        f,_=mechanical_force_energy(x,40,1.2,bonds,rest)
        f+=active; f[fixed]=0; x+=.05/substeps*f/2
        _,u=mechanical_force_energy(x,40,1.2,bonds,rest)
        w=u-float(np.sum(active*x))
        if w-previous>1e-9*max(1,abs(w),abs(previous)): increases+=1
        dist=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
        remove=bonds&(dist>1.6); broken+=int(np.triu(remove,1).sum()); bonds[remove]=False
        _,u=mechanical_force_energy(x,40,1.2,bonds,rest)
        previous=u-float(np.sum(active*x))
        max_overlap=max(max_overlap,float(np.maximum(1-dist[np.triu_indices(len(x),1)],0).max()))
        # Track material connection each inner step, separate from geometric contact.
        min_attached=min(min_attached,float(parent_attachment(bonds,fixed)[~fixed].mean()))
    return x,bonds,dict(broken=broken,increases=increases,max_overlap=max_overlap,min_material_attachment=min_attached)


def observe(x,fixed,d,n,bonds):
    r,candidates=readout(x,fixed,d)
    graph=contact_graph(x); high=d>.6; mobile=~fixed
    isolated=high&~np.any(graph&high[None,:],axis=1)
    return dict(r,material_attachment=float(parent_attachment(bonds,fixed)[mobile].mean()),
        mobile_tip_mean=float(d[mobile].mean()),mobile_tip_std=float(d[mobile].std()),
        mobile_high=int((high&mobile).sum()),mobile_isolated_high=int((isolated&mobile).sum()),
        mobile_notch_mean=float(n[mobile].mean()),junction_edges=int(np.triu(bonds,1).sum())),candidates


def prepattern(x,d,n,competition):
    for _ in range(400): d,n=regulate(x,d,n,competition,.05)
    return d,n


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); root=args.output
    if (root/'settings.json').exists(): raise RuntimeError('Use fresh output')
    root.mkdir(parents=True,exist_ok=True)
    helpers=['branching_sprouting_assay.py','branching_prerequisites.py','branching_adhesion_assay.py','branching_shared_state.py','branching_baseline.py']
    settings=dict(seeds=[42,43,44],substeps=[16,32],outer_dt=.05,horizon=12,
        initializations=['immediate','prepatterned'],prepattern_duration=20,
        conditions=dict(neither=[False,False],guidance=[True,False],competition=[False,True],both=[True,True]),
        repulsion=40,damping=2,spring=4,break_distance=1.6,contact_cutoff=1.2,
        source_sha256=source_hashes(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256={f:hashlib.sha256((Path(__file__).parent/f).read_bytes()).hexdigest() for f in helpers})
    (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    rows=[]; finals=[]; histories={}; states={}
    for seed in settings['seeds']:
        start,fixed,d0=initial(seed)
        rest=np.linalg.norm(start[:,None,:]-start[None,:,:],axis=2); original_bonds=contact_graph(start)
        for initialization in settings['initializations']:
            for steps in settings['substeps']:
                for name,(guidance,competition) in settings['conditions'].items():
                    x=start.copy(); d=d0.copy(); n=np.zeros(36); bonds=original_bonds.copy()
                    if initialization=='prepatterned': d,n=prepattern(x,d,n,competition)
                    initial_d=d.copy(); initial_n=n.copy()
                    positions=[]; ds=[]; ns=[]; bs=[]; times=[]; candidate_history=[]
                    broken=0; increases=0; min_material=1.; min_geometric=1.; max_overlap=float(np.maximum(1-rest[np.triu_indices(36,1)],0).max())
                    for outer in range(241):
                        r,candidates=observe(x,fixed,d,n,bonds)
                        min_geometric=min(min_geometric,r['attached_mobile_fraction'])
                        if outer%10==0:
                            row=dict(seed=seed,initialization=initialization,substeps=steps,condition=name,time=outer*.05,**r)
                            rows.append(row); positions.append(x.copy()); ds.append(d.copy()); ns.append(n.copy()); bs.append(bonds.copy()); times.append(outer*.05); candidate_history.append(candidates)
                        if outer==240: break
                        d,n=regulate(x,d,n,competition,.05)
                        active=traction(x,d,guidance)
                        x,bonds,diagnostics=mechanical_interval(x,fixed,bonds,rest,active,steps)
                        broken+=diagnostics['broken']; increases+=diagnostics['increases']
                        min_material=min(min_material,diagnostics['min_material_attachment']); max_overlap=max(max_overlap,diagnostics['max_overlap'])
                        assert np.isfinite(x).all() and np.isfinite(d).all() and np.isfinite(n).all()
                        np.testing.assert_array_equal(x[fixed],start[fixed])
                    match=lambda a,b:len(a&b)>=.5*min(len(a),len(b))
                    persistent=sum(int(any(match(a,b) and match(b,c) for a in first for b in second for c in third))
                        for first,second,third in zip(candidate_history,candidate_history[1:],candidate_history[2:]))
                    finals.append(dict(row,broken_junctions=broken,active_work_increases=increases,max_overlap=max_overlap,
                        min_material_attachment=min_material,min_outer_geometric_attachment=min_geometric,
                        persistent_candidate_windows=persistent,peak_candidates=max(map(len,candidate_history))))
                    f=root/f'seed{seed}_{initialization}_steps{steps}_{name}.npz'
                    np.savez_compressed(f,positions=positions,dll4=ds,notch=ns,bonds=bs,time=times,fixed=fixed,
                                        initial_dll4=initial_d,initial_notch=initial_n,rest=rest)
                    histories[seed,initialization,steps,name]=(np.array(positions),np.array(ds))
                    states[seed,initialization,steps,name]=(x.copy(),d.copy(),n.copy())
                print(seed,initialization,steps,'completed',flush=True)
    write(root/'timecourses.csv',rows); write(root/'final_summary.csv',finals)
    differences=[]
    for seed in settings['seeds']:
        for initialization in settings['initializations']:
            for name in settings['conditions']:
                a=states[seed,initialization,16,name]; b=states[seed,initialization,32,name]
                differences.append(dict(seed=seed,initialization=initialization,condition=name,
                    position_rms_difference=float(np.sqrt(np.mean(np.sum((a[0]-b[0])**2,axis=1)))),
                    dll4_rmse=float(np.sqrt(np.mean((a[1]-b[1])**2)))))
    write(root/'resolution_differences.csv',differences)
    fig,axes=plt.subplots(2,4,figsize=(14,8),constrained_layout=True)
    for row,initialization in enumerate(settings['initializations']):
        for col,name in enumerate(settings['conditions']):
            x,d=histories[42,initialization,32,name]; end=x[-1]
            for i in range(36): axes[row,col].plot(x[:,i,0],x[:,i,1],color='gray',alpha=.25,lw=.5)
            im=axes[row,col].scatter(end[:,0],end[:,1],c=d[-1],vmin=0,vmax=1,cmap='viridis',s=25)
            axes[row,col].scatter(end[fixed,0],end[fixed,1],c='red',marker='x',s=12)
            axes[row,col].set_title(initialization+'\n'+name,fontsize=10); axes[row,col].set_aspect('equal')
    fig.colorbar(im,ax=axes,label='Abstract tip activity',shrink=.6)
    fig.savefig(root/'junction_sprouting_states.png',dpi=170); plt.close(fig)
    assert settings['source_sha256']==source_hashes()


if __name__=='__main__': main()
