# E013: coupled mechanical refinement, 8 versus 16 substeps

Date: 2026-11-07 (Europe/Berlin). Protocol recorded before the 16-step run completed.

Compare 8 and 16 inner steps per outer interval from the identical E007 seed-44 hour-10 checkpoint through hour 20. Inner dt is 0.00625 versus 0.003125; outer dt remains 0.05. Reuse the archived E012 eight-step control without modifying it. Original application code remains unchanged.

Use exactly the E012 splitting policy: active regulation/fields, arrested growth/division, zero future regulatory noise, polarity and propulsion gate held within each outer interval, first-inner force buffered for the next GRN update. No active-work potential diagnostic is introduced.

Primary comparison: matched-cell position RMS difference and six-node GRN RMSE, at hour 20 and averaged over early (10,15] and late (15,20] observations. Compare these differences with E012 4-versus-8 differences, without treating a smaller difference as proof of convergence. Also record spread, high-SNAIL counts, YAP and force-input histories. High-SNAIL is not a functional tip label.

Validation: identical initial arrays/source/helper hashes; 256 fixed IDs/radii; finite snapshots; equal 200 biological/chemical/polarity outer updates and 51,200 GRN calls/draws/polarity-cell evaluations. Force calls should be 1,600 versus 3,200. Preserve raw snapshots, scripts, input hashes, observations and limitations.

## Results

The run completed to hour 20. All planned validation checks passed, including
42 finite snapshots, identical initial arrays, fixed IDs/radii, source/helper
hashes and update counts. The 6 coupled-substep and 19 checkpoint/mechanics
tests were rerun and passed. Application code was not changed.

| Final readout, hour 20 | 8 steps | 16 steps |
| --- | --- | --- |
| Inner dt | 0.00625 | 0.003125 |
| RMS radius (model length units) | 1.37768 | 1.41362 |
| High-SNAIL cells (SNAIL > 0.8) | 2 | 6 |
| Mean SNAIL | 0.09606 | 0.11174 |
| Mean YAP | 0.29571 | 0.29069 |
| Mean force-proxy input actually read by final GRN update | 1.14215 | 0.64397 |
| Connected components | 5 | 4 |
| Largest-component fraction | 0.96094 | 0.97266 |

RMS spread differs by 2.61% relative to the eight-step result (versus 8.49%
for 4 versus 8). Mean YAP differs by 0.00501, but matched-cell states and
positions still differ. Cell count is not an independent replicate count.

| Matched-cell disagreement | 4 vs 8, E012 | 8 vs 16, E013 |
| --- | --- | --- |
| Final position RMS difference | 0.47745 | 0.49745 |
| Final difference / mean cell radius | 0.63474 | 0.66133 |
| Final six-node GRN RMSE | 0.18859 | 0.08028 |
| Early mean position RMS difference, (10,15] | 0.02665 | 0.02068 |
| Late mean position RMS difference, (15,20] | 0.34018 | 0.42306 |
| Early mean GRN RMSE | 0.00498 | 0.00370 |
| Late mean GRN RMSE | 0.04809 | 0.03535 |

Final GRN RMSE is about 57% smaller, and late mean GRN RMSE about 26%
smaller, with this refinement. Final position disagreement is about 4%
larger, and late mean position disagreement about 24% larger. Therefore
refinement improves some regulatory/collective readouts, but not uniformly.
The early high-SNAIL pulse persists (sampled peaks 27 and 26); the late
endpoint counts of 2 and 6 do not demonstrate branches or stable tips.

## Interpretation and next question

Both fine runs show expansion followed by strong contraction, low late
mean SNAIL, and similar final mean YAP. The large late force-proxy elevation
seen with four steps is absent from both finer sampled trajectories.
This is evidence for greater consistency of some coarse observables, not
proof that the 16-step trajectory is correct. Position-level agreement
does not improve monotonically, and the final force inputs still differ.
No production timestep or model mechanism is adopted from this test.

The comparison cannot distinguish numerical integration error from
sensitive trajectories driven by nonlinear contact/neighbor changes. The
fixed outer clock, held polarity/eligibility, force-proxy definition and
feedback lag remain untested sources of error. One seed/checkpoint with
growth arrested and noise removed does not validate biological branching.

A useful next diagnostic is the within-outer active-work-adjusted potential
check proposed in E012. It could test whether individual inner updates
respect the fixed-force interval's dissipation property, without demanding
that an actively regulated system lose energy globally. This has not been
implemented or run; additional blind timestep halving is not treated as a
substitute for that check.

## Reproduction and artifacts

```powershell
rtk proxy python tools/branching_coupled_substeps.py --checkpoint research/experiments/2026-11-06-shared-state-timesteps/checkpoint --outer-dt 0.05 --substeps 16 --horizon 20 --sample-every 0.5 --output _branching_coupled_substeps16
rtk proxy python tools/branching_coupled_substeps_report.py research/experiments/2026-11-06-coupled-substeps/steps4 research/experiments/2026-11-07-coupled-refinement/steps8 research/experiments/2026-11-07-coupled-refinement/steps16 --output research/experiments/2026-11-07-coupled-refinement/analysis
```

Use fresh output folders to rerun. `steps8` is an unchanged copy of E012;
`steps16` contains the new run. Analysis also references the archived E012
four-step arm, retained there rather than duplicated. `analysis` contains
the figure, final summary, half-hour pairwise differences, window summary
and report provenance. Archived scripts, `validation.json` and `sha256.json`
preserve implementation, checks and hashes. Pressure/velocity timing follows
E012: GRN inputs are lagged outer-start buffers; last-inner pressure and
velocity are separate, and initial force inputs are inherited checkpoint
values. Force proxies are not biological stresses. The figure's pressure
axis uses symlog scaling.
