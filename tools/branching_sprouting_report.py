"""Summarize E015, including failed sprout criteria and connectivity sensitivity."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from branching_sprouting_assay import components,contact_graph,write
from branching_shared_state import source_hashes


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('folder',type=Path); args=parser.parse_args()
    p=args.folder
    rows=list(csv.DictReader((p/'timecourses.csv').open()))
    finals=list(csv.DictReader((p/'final_summary.csv').open()))
    settings=json.loads((p/'settings.json').read_text())
    assert len(finals)==24 and len(rows)==600
    assert settings['source_sha256']==source_hashes()
    summary=[]; sensitivity=[]
    for condition in settings['conditions']:
        selected=[r for r in finals if r['condition']==condition and int(r['substeps'])==32]
        record=dict(condition=condition,runs=len(selected))
        for key in ['mobile_mean_x','attached_mobile_fraction','max_attached_extension',
                    'detached_mobile_cells','tip_activity_std','high_tip_activity_cells',
                    'persistent_candidate_windows']:
            a=np.array([float(r[key]) for r in selected])
            record.update({key+'_mean':float(a.mean()),key+'_min':float(a.min()),key+'_max':float(a.max())})
        chosen=[r for r in rows if r['condition']==condition and int(r['substeps'])==32]
        record['peak_sampled_attached_extension']=max(float(r['max_attached_extension']) for r in chosen)
        record['peak_sampled_candidate_count']=max(int(r['narrow_extension_candidates']) for r in chosen)
        mobile_means=[]; mobile_stds=[]
        for seed in settings['seeds']:
            with np.load(p/f'seed{seed}_steps32_{condition}.npz') as z:
                d=z['dll4'][-1,~z['fixed']]
                mobile_means.append(float(d.mean())); mobile_stds.append(float(d.std()))
        record['mobile_tip_activity_mean']=float(np.mean(mobile_means))
        record['mobile_tip_activity_std_mean']=float(np.mean(mobile_stds))
        summary.append(record)
    hashes={}
    for seed in settings['seeds']:
        with np.load(p/f'initial_seed{seed}.npz') as z:
            start={k:z[k].copy() for k in z.files}
        for n in settings['substeps']:
            for condition in settings['conditions']:
                f=p/f'seed{seed}_steps{n}_{condition}.npz'
                hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
                with np.load(f) as z:
                    assert all(np.isfinite(z[k]).all() for k in z.files)
                    np.testing.assert_array_equal(z['positions'][0],start['positions'])
                    np.testing.assert_array_equal(z['dll4'][0],start['dll4'])
                    np.testing.assert_array_equal(z['notch'][0],start['notch'])
                    assert z['positions'].shape==(25,36,3)
                    assert np.all((z['dll4']>=0)&(z['dll4']<=1))
                    assert np.all((z['notch']>=0)&(z['notch']<=1))
                    np.testing.assert_array_equal(z['positions'][:,z['fixed']],np.broadcast_to(start['positions'][z['fixed']],(25,12,3)))
                    end=z['positions'][-1]; fixed=z['fixed']
                for cutoff in [1.1,1.2,1.5,2.5]:
                    graph=contact_graph(end,cutoff); attached=np.zeros(36,bool)
                    for g in components(graph):
                        if fixed[g].any(): attached[g]=True
                    sensitivity.append(dict(seed=seed,substeps=n,condition=condition,cutoff=cutoff,
                                            attached_mobile_fraction=float(attached[~fixed].mean())))
    write(p/'condition_summary.csv',summary); write(p/'connectivity_sensitivity.csv',sensitivity)
    diagnostics=dict(runs=24,snapshots=600,outer_updates_per_run=240,
        inner_updates_total=3*4*240*(16+32),application_source_unchanged=True,
        matched_initial_states=True,fixed_anchors=True,finite=True,
        max_overlap=max(float(r['max_overlap']) for r in finals),
        active_work_increases=sum(int(r['active_work_potential_increases']) for r in finals),
        note='Clock counts follow explicit loop bounds; snapshot shape and initial/fixed arrays independently checked.',
        trajectory_sha256=hashes,
        report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (p/'validation.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    for condition in settings['conditions']:
        selected=[r for r in rows if r['condition']==condition and int(r['substeps'])==32]
        times=sorted({float(r['time']) for r in selected})
        for ax,key in zip(axes,['attached_mobile_fraction','max_attached_extension','tip_activity_std']):
            values=np.array([[float(r[key]) for r in selected if float(r['time'])==t] for t in times])
            ax.plot(times,values.mean(axis=1),label=condition)
            ax.fill_between(times,values.min(axis=1),values.max(axis=1),alpha=.12)
            ax.set(xlabel='Model time',ylabel=key.replace('_',' ')); ax.grid(alpha=.2)
    axes[0].legend(fontsize=7); fig.savefig(p/'timecourse_comparison.png',dpi=170); plt.close(fig)
    for r in summary: print(r)
    print('active-work increases',diagnostics['active_work_increases'],'max overlap',diagnostics['max_overlap'])


if __name__=='__main__': main()
