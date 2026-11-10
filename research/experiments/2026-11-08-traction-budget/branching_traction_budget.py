"""Matched-total-force diagnostic of follower load and activity-weighted traction."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from branching_sprouting_assay import initial,contact_graph,regulate,write
from branching_junction_sprouting import prepattern,mechanical_interval,observe,parent_attachment
from branching_shared_state import source_hashes


def allocate(d,fixed,budget,follower_fraction):
    """Global force budget is a diagnostic intervention, not a local biological rule."""
    weights=follower_fraction+(1-follower_fraction)*np.maximum(d[~fixed],0)**4
    fallback=float(weights.sum())<1e-15
    if fallback: weights=np.ones_like(weights)
    f=np.zeros((len(d),3)); f[~fixed,0]=budget*weights/weights.sum()
    return f,fallback


def candidates(x,fixed,bonds):
    from branching_sprouting_assay import components
    attached=parent_attachment(bonds,fixed); found=[]
    for group in components(bonds,attached&(x[:,0]>0)):
        if len(group)>=3 and x[group,0].max()>=2 and np.ptp(x[group,1])<=2:
            found.append(set(group))
    return found


def persistent(history):
    match=lambda a,b:len(a&b)>=.5*min(len(a),len(b))
    return sum(int(any(match(a,b) and match(b,c) for a in first for b in second for c in third))
               for first,second,third in zip(history,history[1:],history[2:]))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True); args=parser.parse_args(); root=args.output
    if (root/'settings.json').exists(): raise RuntimeError('Use fresh output')
    root.mkdir(parents=True,exist_ok=True)
    conditions={f'F{budget}_b{b}':dict(budget=budget,follower_fraction=b) for budget in [2.5,7.5,15.] for b in [0.,.25,1.]}
    conditions['E018_reference']=dict(budget=None,follower_fraction=None)
    helpers=['branching_junction_sprouting.py','branching_prerequisites.py','branching_sprouting_assay.py','branching_adhesion_assay.py','branching_shared_state.py','branching_baseline.py']
    settings=dict(seeds=[42,43,44],substeps=[16,32],outer_dt=.05,horizon=12,prepattern_duration=20,
        conditions=conditions,weight_rule='b+(1-b)*D**4, normalized across mobile cells to fixed total F; +x direction',
        source_sha256=source_hashes(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256={f:hashlib.sha256((Path(__file__).parent/f).read_bytes()).hexdigest() for f in helpers})
    (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    rows=[]; finals=[]; endpoints={}
    for seed in settings['seeds']:
        start,fixed,d0=initial(seed); rest=np.linalg.norm(start[:,None,:]-start[None,:,:],axis=2)
        d_initial,n_initial=prepattern(start,d0.copy(),np.zeros(36),True)
        for steps in settings['substeps']:
            for name,condition in conditions.items():
                x=start.copy(); d=d_initial.copy(); n=n_initial.copy(); bonds=contact_graph(start)
                positions=[]; ds=[]; ns=[]; bs=[]; times=[]; force_history=[]; geo_history=[]; mat_history=[]
                increases=0; broken=0; min_attachment=1.; max_overlap=float(np.maximum(1-rest[np.triu_indices(36,1)],0).max())
                min_geo=1.; fallbacks=0; max_budget_error=0.; total_trace=[]; max_cell_trace=[]
                last_active=np.zeros_like(x)
                for outer in range(241):
                    r,geo=observe(x,fixed,d,n,bonds); mat=candidates(x,fixed,bonds)
                    attached=parent_attachment(bonds,fixed)
                    r.update(material_attached_max_x=float(max(0,x[attached,0].max())),material_narrow_candidates=len(mat))
                    min_geo=min(min_geo,r['attached_mobile_fraction'])
                    if outer%10==0:
                        row=dict(seed=seed,substeps=steps,condition=name,time=outer*.05,
                            **r,last_interval_total_force=float(last_active[:,0].sum()),last_interval_max_cell_force=float(last_active[:,0].max()))
                        rows.append(row); positions.append(x.copy()); ds.append(d.copy()); ns.append(n.copy()); bs.append(bonds.copy()); times.append(outer*.05)
                        force_history.append(last_active.copy()); geo_history.append(geo); mat_history.append(mat)
                    if outer==240: break
                    d,n=regulate(x,d,n,True,.05)
                    if name=='E018_reference':
                        active=np.zeros_like(x); active[~fixed,0]=.15+2.35*d[~fixed]
                    else:
                        active,fallback=allocate(d,fixed,condition['budget'],condition['follower_fraction']); fallbacks+=int(fallback)
                        max_budget_error=max(max_budget_error,abs(float(active[:,0].sum())-condition['budget']))
                    total_trace.append(float(active[:,0].sum())); max_cell_trace.append(float(active[:,0].max()))
                    x,bonds,stats=mechanical_interval(x,fixed,bonds,rest,active,steps)
                    last_active=active.copy(); increases+=stats['increases']; broken+=stats['broken']
                    min_attachment=min(min_attachment,stats['min_material_attachment']); max_overlap=max(max_overlap,stats['max_overlap'])
                    assert np.isfinite(x).all() and np.isfinite(d).all() and np.isfinite(n).all()
                    np.testing.assert_array_equal(x[fixed],start[fixed])
                finals.append(dict(row,min_material_attachment=min_attachment,min_outer_geometric_attachment=min_geo,
                    max_overlap=max_overlap,broken_junctions=broken,active_work_increases=increases,
                    persistent_geometric_windows=persistent(geo_history),persistent_material_windows=persistent(mat_history),
                    peak_geometric_candidates=max(map(len,geo_history)),peak_material_candidates=max(map(len,mat_history)),
                    allocation_fallbacks=fallbacks,max_budget_error=max_budget_error,
                    mean_total_force=float(np.mean(total_trace)),max_total_force=max(total_trace),max_single_cell_force=max(max_cell_trace)))
                np.savez_compressed(root/f'seed{seed}_steps{steps}_{name}.npz',positions=positions,dll4=ds,notch=ns,bonds=bs,time=times,
                                    fixed=fixed,last_interval_active=force_history,rest=rest)
                endpoints[seed,steps,name]=(x.copy(),d.copy())
                print(seed,steps,name,'attachment',r['material_attachment'],'extension',r['material_attached_max_x'],flush=True)
    write(root/'timecourses.csv',rows); write(root/'final_summary.csv',finals)
    differences=[]
    for seed in settings['seeds']:
        for name in conditions:
            a=endpoints[seed,16,name]; b=endpoints[seed,32,name]
            differences.append(dict(seed=seed,condition=name,position_rms_difference=float(np.sqrt(np.mean(np.sum((a[0]-b[0])**2,axis=1)))),
                                    dll4_rmse=float(np.sqrt(np.mean((a[1]-b[1])**2)))))
    write(root/'resolution_differences.csv',differences)
    assert settings['source_sha256']==source_hashes()


if __name__=='__main__': main()
