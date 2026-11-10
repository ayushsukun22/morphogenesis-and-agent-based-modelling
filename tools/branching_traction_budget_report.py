"""Validate the fixed-budget intervention and summarize all E019 conditions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_sprouting_assay import initial,write
from branching_junction_sprouting import prepattern
from branching_shared_state import source_hashes


def main():
    p=argparse.ArgumentParser(); p.add_argument('folder',type=Path); args=p.parse_args(); root=args.folder
    settings=json.loads((root/'settings.json').read_text()); assert settings['source_sha256']==source_hashes()
    finals=list(csv.DictReader((root/'final_summary.csv').open())); rows=list(csv.DictReader((root/'timecourses.csv').open()))
    assert len(finals)==60 and len(rows)==1500
    hashes={}; refs={}; count=0
    for seed in settings['seeds']:
        start,fixed,d0=initial(seed); d,n=prepattern(start,d0,np.zeros(36),True)
        for steps in settings['substeps']:
            for name,condition in settings['conditions'].items():
                f=root/f'seed{seed}_steps{steps}_{name}.npz'
                with np.load(f) as z:
                    assert z['positions'].shape==(25,36,3)
                    assert all(np.isfinite(z[k]).all() for k in z.files)
                    np.testing.assert_array_equal(z['positions'][0],start)
                    np.testing.assert_array_equal(z['dll4'][0],d); np.testing.assert_array_equal(z['notch'][0],n)
                    np.testing.assert_array_equal(z['positions'][:,fixed],np.broadcast_to(start[fixed],(25,12,3)))
                    assert np.array_equal(z['bonds'],np.swapaxes(z['bonds'],1,2))
                    assert not np.any(z['bonds'][1:]&~z['bonds'][:-1])
                    active=z['last_interval_active'][1:]
                    np.testing.assert_array_equal(active[:,fixed],0)
                    np.testing.assert_array_equal(active[:,:,1:],0)
                    if condition['budget'] is not None:
                        np.testing.assert_allclose(active[:,:,0].sum(axis=1),condition['budget'],rtol=0,atol=1e-12)
                    else:
                        reference=Path('research/experiments/2026-11-08-junction-sprouting')/f'seed{seed}_prepatterned_steps{steps}_both.npz'
                        with np.load(reference) as old:
                            for key in ['positions','dll4','notch','bonds','time','fixed','rest']:
                                np.testing.assert_array_equal(z[key],old[key])
                        refs[str(reference)]=hashlib.sha256(reference.read_bytes()).hexdigest()
                    count+=25
                hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
    summary=[]
    for name,condition in settings['conditions'].items():
        selected=[r for r in finals if r['condition']==name and int(r['substeps'])==32]
        record=dict(condition=name,**condition)
        for key in ['material_attachment','attached_mobile_fraction','material_attached_max_x','min_material_attachment',
                    'broken_junctions','max_overlap','mobile_high','mobile_tip_std','persistent_geometric_windows',
                    'persistent_material_windows','peak_geometric_candidates','peak_material_candidates',
                    'mean_total_force','max_total_force','max_single_cell_force']:
            a=np.array([float(r[key]) for r in selected])
            record.update({key+'_mean':float(a.mean()),key+'_min':float(a.min()),key+'_max':float(a.max())})
        chosen=[r for r in rows if r['condition']==name and int(r['substeps'])==32]
        record['peak_sampled_material_extension']=max(float(r['material_attached_max_x']) for r in chosen)
        summary.append(record)
    write(root/'condition_summary.csv',summary)
    validation=dict(runs=60,snapshots=count,finite=True,anchors_fixed=True,source_unchanged=True,
        initial_states_matched=True,reference_replay_exact=True,reference_input_sha256=refs,
        active_work_increases=sum(int(r['active_work_increases']) for r in finals),
        maximum_overlap=max(float(r['max_overlap']) for r in finals),
        maximum_budget_error=max(float(r['max_budget_error']) for r in finals),
        allocation_fallbacks=sum(int(r['allocation_fallbacks']) for r in finals),
        trajectory_sha256=hashes,report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (root/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for ax,key,title in zip(axes,['material_attachment_mean','material_attached_max_x_mean','broken_junctions_mean'],
                           ['Final material attachment','Final material-connected max x','Broken junctions']):
        data=np.array([[next(r[key] for r in summary if r['budget']==F and r['follower_fraction']==b)
                        for b in [0,.25,1]] for F in [2.5,7.5,15.]])
        im=ax.imshow(data,aspect='auto',cmap='viridis')
        for i in range(3):
            for j in range(3): ax.text(j,i,f'{data[i,j]:.2f}',ha='center',va='center',color='white' if data[i,j]<data.max()*.5 else 'black')
        ax.set(xticks=range(3),xticklabels=['b=0\nactivity weighted','b=.25\nfollower floor','b=1\nuniform'],
               yticks=range(3),yticklabels=['2.5','7.5','15'],xlabel='Allocation b',ylabel='Total traction budget',title=title)
        ax.tick_params(axis='x',labelsize=8); fig.colorbar(im,ax=ax,shrink=.7)
    fig.savefig(root/'traction_budget_summary.png',dpi=170); plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(10,9),constrained_layout=True)
    for i,F in enumerate([2.5,7.5,15.]):
        for j,b in enumerate([0.,.25,1.]):
            with np.load(root/f'seed42_steps32_F{F}_b{b}.npz') as z:
                x=z['positions'][-1]; bonds=z['bonds'][-1]; fixed=z['fixed']
                for a,c in zip(*np.where(np.triu(bonds,1))): axes[i,j].plot(x[[a,c],0],x[[a,c],1],color='steelblue',lw=.6)
                im=axes[i,j].scatter(x[:,0],x[:,1],c=z['dll4'][-1],vmin=0,vmax=1,cmap='viridis',s=25)
                axes[i,j].scatter(x[fixed,0],x[fixed,1],c='red',marker='x',s=12)
                axes[i,j].set_title(f'F={F}, b={b}',fontsize=10); axes[i,j].set_aspect('equal')
    fig.colorbar(im,ax=axes,label='Abstract tip activity',shrink=.6)
    fig.savefig(root/'traction_budget_states.png',dpi=170); plt.close(fig)
    for r in summary: print(r)
    print('Validation:',{k:v for k,v in validation.items() if 'sha256' not in k})


if __name__=='__main__': main()
