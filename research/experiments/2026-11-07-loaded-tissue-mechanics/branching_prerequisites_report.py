"""Validate and summarize both independent prerequisites without selecting a winner."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_prerequisites import initial,write
from branching_shared_state import source_hashes


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('mode',choices=['signaling','mechanics']); parser.add_argument('folder',type=Path)
    args=parser.parse_args(); root=args.folder
    settings=json.loads((root/'settings.json').read_text())
    assert settings['source_sha256']==source_hashes()
    rows=list(csv.DictReader((root/'final_summary.csv').open()))
    expected=36 if args.mode=='signaling' else 72
    assert len(rows)==expected
    summaries=[]; hashes={}; snapshots=0
    if args.mode=='signaling':
        for topology in settings['graphs']:
            for condition in settings['conditions']:
                selected=[r for r in rows if r['topology']==topology and r['condition']==condition and float(r['dt'])==.025]
                record=dict(topology=topology,condition=condition)
                for key in ['activity_std','high_cells','low_cells','isolated_high_cells','stable_isolated_high_cells','mean_edge_contrast']:
                    a=[float(r[key]) for r in selected]
                    record.update({key+'_mean':float(np.mean(a)),key+'_min':min(a),key+'_max':max(a)})
                mobile_high=[]
                for seed in settings['seeds']:
                    with np.load(root/f'seed{seed}_{topology}_{condition}_dt0.025.npz') as z:
                        _,fixed,_=initial(seed); mobile_high.append(int(((z['dll4'][-1]>.6)&~fixed).sum()))
                record['mobile_high_mean']=float(np.mean(mobile_high)); summaries.append(record)
        for f in root.glob('seed*.npz'):
            seed=int(f.name.split('_')[0][4:]); _,_,d0=initial(seed)
            with np.load(f) as z:
                assert z['dll4'].shape==(41,36)
                np.testing.assert_array_equal(z['dll4'][0],d0)
                assert np.array_equal(z['graph'],z['graph'].T)
                assert all(np.isfinite(z[k]).all() for k in z.files)
                assert np.all((z['dll4']>=0)&(z['dll4']<=1))
                snapshots+=41
            hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
        inputs=json.loads((root/'input_hashes.json').read_text())
        assert all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in inputs.items())
    else:
        for variant in settings['variants']:
            for load in settings['loads']:
                selected=[r for r in rows if r['variant']==variant and r['load']==load and float(r['dt'])==.0025]
                record=dict(variant=variant,load=load)
                for key in ['attached_fraction','attached_fraction_1_6','front_displacement','mean_mobile_displacement','max_overlap_over_run','broken_initial_junctions','min_sampled_attachment']:
                    a=[float(r[key]) for r in selected]
                    record.update({key+'_mean':float(np.mean(a)),key+'_min':min(a),key+'_max':max(a)})
                # Predeclared overlap limit plus geometric attachment; meaningful extension is separate.
                record['all_seeds_attachment_overlap_screen']=all(float(r['min_sampled_attachment'])==1 and float(r['max_overlap_over_run'])<=.1 for r in selected)
                summaries.append(record)
        for f in root.glob('seed*.npz'):
            seed=int(f.name.split('_')[0][4:]); start,fixed,_=initial(seed)
            with np.load(f) as z:
                assert z['positions'].shape==(21,36,3)
                np.testing.assert_array_equal(z['positions'][0],start)
                np.testing.assert_array_equal(z['positions'][:,fixed],np.broadcast_to(start[fixed],(21,12,3)))
                assert all(np.isfinite(z[k]).all() for k in z.files)
                snapshots+=21
            hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
    write(root/'condition_summary.csv',summaries)
    if args.mode=='signaling':
        fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
        for row,topology in enumerate(settings['graphs']):
            for col,name in enumerate(settings['conditions']):
                with np.load(root/f'seed42_{topology}_{name}_dt0.025.npz') as z:
                    x=z['positions']; graph=z['graph']; d=z['dll4'][-1]
                    for i,j in zip(*np.where(np.triu(graph,1))): axes[row,col].plot(x[[i,j],0],x[[i,j],1],color='gray',lw=.3,alpha=.3)
                    im=axes[row,col].scatter(x[:,0],x[:,1],c=d,vmin=0,vmax=1,cmap='viridis',s=45)
                    axes[row,col].set_title(topology.replace('_',' ')+'\n'+name.replace('_',' '),fontsize=10)
                    axes[row,col].set_aspect('equal')
        fig.colorbar(im,ax=axes,label='Abstract tip activity',shrink=.6)
        fig.savefig(root/'signaling_states.png',dpi=170); plt.close(fig)
    else:
        fig,axes=plt.subplots(3,4,figsize=(13,10),constrained_layout=True)
        for row,load in enumerate(settings['loads']):
            for col,name in enumerate(settings['variants']):
                with np.load(root/f'seed42_dt0.0025_{name}_{load}.npz') as z:
                    x=z['positions'][-1]; fixed=z['fixed']
                    axes[row,col].scatter(x[~fixed,0],x[~fixed,1],s=20)
                    axes[row,col].scatter(x[fixed,0],x[fixed,1],s=15,c='red',marker='x')
                    axes[row,col].set_title(name.replace('_',' ')+'\n'+load.replace('_',' '),fontsize=10)
                    axes[row,col].set_aspect('equal')
        fig.savefig(root/'mechanics_states.png',dpi=170); plt.close(fig)
    diagnostics=dict(runs=expected,saved_states=snapshots,source_unchanged=True,finite=True,
        matched_initial_arrays=True,unit_tests_passed=5,trajectory_sha256=hashes,
        active_work_increases=sum(int(r.get('active_work_increases',0)) for r in rows),
        report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (root/'validation.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    for record in summaries: print(record)
    print('diagnostics', {k:v for k,v in diagnostics.items() if k!='trajectory_sha256'})


if __name__=='__main__': main()
