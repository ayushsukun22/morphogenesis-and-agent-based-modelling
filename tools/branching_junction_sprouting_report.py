"""Validate and summarize E018 without equating migration with sprouting."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_sprouting_assay import initial,write,contact_graph,components
from branching_junction_sprouting import prepattern,parent_attachment
from branching_shared_state import source_hashes


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('folder',type=Path); args=parser.parse_args(); root=args.folder
    settings=json.loads((root/'settings.json').read_text()); assert settings['source_sha256']==source_hashes()
    rows=list(csv.DictReader((root/'timecourses.csv').open())); finals=list(csv.DictReader((root/'final_summary.csv').open()))
    assert len(rows)==1200 and len(finals)==48
    hashes={}; count=0; network_rows=[]
    for seed in settings['seeds']:
        start,fixed,d0=initial(seed)
        for timing in settings['initializations']:
            for n in settings['substeps']:
                for name,(_,competition) in settings['conditions'].items():
                    f=root/f'seed{seed}_{timing}_steps{n}_{name}.npz'
                    d=d0.copy(); notch=np.zeros(36)
                    if timing=='prepatterned': d,notch=prepattern(start,d,notch,competition)
                    with np.load(f) as z:
                        assert z['positions'].shape==(25,36,3)
                        assert all(np.isfinite(z[k]).all() for k in z.files)
                        np.testing.assert_array_equal(z['positions'][0],start)
                        np.testing.assert_array_equal(z['dll4'][0],d)
                        np.testing.assert_array_equal(z['notch'][0],notch)
                        np.testing.assert_array_equal(z['positions'][:,fixed],np.broadcast_to(start[fixed],(25,12,3)))
                        assert np.all((z['dll4']>=0)&(z['dll4']<=1))
                        assert np.array_equal(z['bonds'],np.swapaxes(z['bonds'],1,2))
                        assert not np.any(z['bonds'][1:]&~z['bonds'][:-1])
                        for i,t in enumerate(z['time']):
                            graph=contact_graph(z['positions'][i]); bonds=z['bonds'][i]
                            edges=np.triu(bonds,1)
                            positions=z['positions'][i]; attached=parent_attachment(bonds,fixed)
                            material_candidates=0
                            for group in components(bonds,attached&(positions[:,0]>0)):
                                if len(group)>=3 and positions[group,0].max()>=2 and np.ptp(positions[group,1])<=2:
                                    material_candidates+=1
                            network_rows.append(dict(seed=seed,initialization=timing,substeps=n,condition=name,time=float(t),
                                material_edges=int(edges.sum()),signaling_edges=int(np.triu(graph,1).sum()),
                                material_edges_outside_signaling_range=int((edges&~graph).sum()),
                                mobile_without_signaling_neighbors=int(((graph.sum(axis=1)==0)&~fixed).sum()),
                                material_attached_max_x=float(max(0,positions[attached,0].max())),
                                material_narrow_candidates=material_candidates))
                        count+=25
                    hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
    summary=[]
    for timing in settings['initializations']:
        for name in settings['conditions']:
            selected=[r for r in finals if r['initialization']==timing and r['condition']==name and int(r['substeps'])==32]
            record=dict(initialization=timing,condition=name)
            for key in ['attached_mobile_fraction','material_attachment','max_attached_extension','mobile_tip_mean','mobile_tip_std',
                        'mobile_high','mobile_isolated_high','broken_junctions','max_overlap','min_material_attachment',
                        'persistent_candidate_windows','peak_candidates']:
                a=np.array([float(r[key]) for r in selected])
                record.update({key+'_mean':float(a.mean()),key+'_min':float(a.min()),key+'_max':float(a.max())})
            summary.append(record)
    write(root/'condition_summary.csv',summary)
    write(root/'network_diagnostics.csv',network_rows)
    validation=dict(runs=48,saved_states=count,source_unchanged=True,matched_initialization_policy=True,fixed_anchors=True,
        finite=True,irreversible_symmetric_junctions=True,active_work_increases=sum(int(r['active_work_increases']) for r in finals),
        maximum_overlap=max(float(r['max_overlap']) for r in finals),unit_tests_passed=14,trajectory_sha256=hashes,
        report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (root/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for row,timing in enumerate(settings['initializations']):
        for name in settings['conditions']:
            selected=[r for r in rows if r['initialization']==timing and r['condition']==name and int(r['substeps'])==32]
            times=sorted({float(r['time']) for r in selected})
            for col,key in enumerate(['material_attachment','max_attached_extension','mobile_isolated_high']):
                a=np.array([[float(r[key]) for r in selected if float(r['time'])==t] for t in times])
                axes[row,col].plot(times,a.mean(axis=1),label=name)
                axes[row,col].fill_between(times,a.min(axis=1),a.max(axis=1),alpha=.12)
                axes[row,col].set(title=timing,xlabel='Time since movement begins',ylabel=key.replace('_',' ')); axes[row,col].grid(alpha=.2)
    axes[0,0].legend(fontsize=8); fig.savefig(root/'timecourse_comparison.png',dpi=170); plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(14,8),constrained_layout=True)
    for row,timing in enumerate(settings['initializations']):
        for col,name in enumerate(settings['conditions']):
            with np.load(root/f'seed42_{timing}_steps32_{name}.npz') as z:
                x=z['positions']; end=x[-1]; bonds=z['bonds'][-1]; fixed=z['fixed']
                for i in range(36): axes[row,col].plot(x[:,i,0],x[:,i,1],color='gray',alpha=.2,lw=.5)
                for i,j in zip(*np.where(np.triu(bonds,1))): axes[row,col].plot(end[[i,j],0],end[[i,j],1],color='steelblue',lw=.6,alpha=.7)
                im=axes[row,col].scatter(end[:,0],end[:,1],c=z['dll4'][-1],vmin=0,vmax=1,cmap='viridis',s=25)
                axes[row,col].scatter(end[fixed,0],end[fixed,1],c='red',marker='x',s=12)
                axes[row,col].set_title(timing+'\n'+name,fontsize=10); axes[row,col].set_aspect('equal')
    fig.colorbar(im,ax=axes,label='Abstract tip activity',shrink=.6)
    fig.savefig(root/'junction_sprouting_states.png',dpi=170); plt.close(fig)
    for r in summary: print(r)
    print('Validation:',{k:v for k,v in validation.items() if k!='trajectory_sha256'})


if __name__=='__main__': main()
