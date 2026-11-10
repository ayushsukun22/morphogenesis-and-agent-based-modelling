# E012: coupled regulation with mechanical substeps

Date: 2026-11-06 (Europe/Berlin). Internal protocol recorded before execution.

## Question and design

Do finer mechanical steps improve agreement when regulation and chemical
feedback are active, with their update timing held fixed? Restore the SAME
complete E007 seed-44 hour-10 checkpoint, not E011's reduced final state.
Run to hour 20, keeping outer regulatory/field dt=0.05 model hours and
growth/division arrested. Zero GRN noise as in E007. Compare 1, 4 and 8
mechanical steps per outer interval (inner dt=0.05, 0.0125 and 0.00625).

## Explicit operator-splitting policy

- The original simulation loop samples fields, integrates each GRN, synchronizes
  phenotype, handles age/fate and chemical consumption/deposition once per
  OUTER interval. Field diffusion/decay advances once afterward at outer dt.
- At the outer mechanics start, evaluate the original neighbor signals and
  polarity turn ONCE, with original per-update memory 0.2.
- Hold polarity direction, propulsion eligibility and current GRN-derived
  strengths fixed during the inner force steps. Eligibility uses the original
  no-neighbor/symmetric-direction guards at outer start; it is not recomputed
  until the next outer interval. Pair forces are recalculated at every inner
  step using updated positions and the original force arithmetic/order.
- Preserve original GRN input timing: cell.pressure stores the first inner
  force norm (outer-start force), while neighbor_signal also comes from outer
  start. The next GRN update sees these buffers with the original outer lag.
  The last inner force norm and velocity are recorded separately and do not
  silently replace the GRN buffer with a newer force sample.
- sim.time advances only on the original outer clock. Inner steps do not
  cause extra regulatory, field or polarity updates.

This is an explicit split diagnostic, not a claim that it uniquely defines
continuous polarity dynamics. It avoids confounding refinement with more
frequent GRN/field updates or a changed feedback lag. Application code and
pair-force laws remain unchanged. Copied mechanics phases must be checked
against the original class, including an EXACT one-substep replay of E007.

## Planned observations and checks

Compare half-hour positions, regulatory states, high-SNAIL counts, fields,
spread, connectivity and YAP/force-input histories. Primary resolution
comparison is 4 versus 8 substeps; retain 1-step control to diagnose the
effect relative to the original coupled scheme. Examine full timecourses,
early (10,15] and late (15,20], not a desired final branching appearance.

All arms must share the checkpoint, initial arrays and source hashes; keep
256 fixed IDs/radii. Each arm must have 200 outer/field updates, 51,200
GRN calls, midpoint draws and polarity-cell evaluations; mechanical calls
are 200/800/1600. The 1-step arm must match all original E007 saved arrays
exactly for shared keys. Phase and pressure-buffer tests must pass.

No potential-decrease requirement is applied to active, changing phenotypes.
E011's passive criterion is not valid as a full-model energy constraint.
No universal production timestep, biological validation or branch claim
is selected from this single checkpoint and artificial growth arrest.

## Implementation validation

Six new tests passed: one-substep full-update equality with the original;
held polarity/neighbor signals during inner steps; restoration of the FIRST
force to the regulatory buffer instead of the last; fixed GRN/field/polarity
clock counts; original symmetric/isolated propulsion guards; and invalid
substep-count rejection. The nineteen earlier checkpoint/mechanics tests
were rerun and passed. A separate 256-cell checkpoint one-substep check
matched source positions, velocities, polarity, pressure and neighbor signals
exactly. The full 1-step control matches all 21 E007 observations exactly
for every array E007 stored, including both chemical grids and GRN nodes.

The helper preserves pair-force arithmetic and cell accumulation order;
it separates the once-per-outer polarity decision from the force integration.
Its 1-step replay validates this split at the reference resolution, not the
physical accuracy of the finer operator-split model.

## Readout conventions

Raw pressure is the outer-start regulatory buffer. grn_pressure_inputs is
what the CURRENT GRN update actually read from the PREVIOUS buffer. The
separate last_inner_pressure and velocity describe the latest inner force
evaluation, not an outer-interval average velocity. At the initial observation
these inherited force values come from the pre-intervention checkpoint.
propulsion_active is -1 initially (prior gate not reconstructed), then 0/1
for the outer-start decision held through the inner interval.

The pressure plot uses symlog scaling (linear near zero, logarithmic above
one) to show large spikes without hiding lower values. Mean input force is
not a biological stress. For interpretation the report also computes the
fraction of final GRN pressure inputs above 1 and mean clip(pressure-1,0,1),
matching the existing YAP driving term. A mean-force threshold is not used
to classify individual cells, and these readouts do not establish a unique
cause of any regulatory-state difference.

Commands from repository root with fresh output folders:

```powershell
rtk proxy python tools/branching_coupled_substeps.py --checkpoint research/experiments/2026-11-06-shared-state-timesteps/checkpoint --outer-dt 0.05 --substeps 1 --horizon 20 --sample-every 0.5 --output _branching_coupled_substeps1
rtk proxy python tools/branching_coupled_substeps.py --checkpoint research/experiments/2026-11-06-shared-state-timesteps/checkpoint --outer-dt 0.05 --substeps 4 --horizon 20 --sample-every 0.5 --output _branching_coupled_substeps4
rtk proxy python tools/branching_coupled_substeps.py --checkpoint research/experiments/2026-11-06-shared-state-timesteps/checkpoint --outer-dt 0.05 --substeps 8 --horizon 20 --sample-every 0.5 --output _branching_coupled_substeps8
rtk proxy python tools/test_branching_coupled_substeps.py
rtk proxy python tools/test_branching_shared_state.py
```

## Results

All three runs completed from hour 10 to 20. Each retained the same 256
cell IDs/radii and performed 200 outer updates, 200 updates of each field,
51,200 GRN calls, 51,200 midpoint draws and 51,200 polarity-cell evaluations.
Only inner force calls differed: 200, 800 and 1,600. All 63 saved snapshots
were finite; initial arrays matched exactly across arms. Application-source
hashes remained unchanged. The 25 tests described above passed.

Final observations at model hour 20:

| Inner steps | Inner dt | RMS radius | High-SNAIL cells | Mean SNAIL | Mean YAP | Mean GRN force-proxy input |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.05 | 3.05241 | 34 | 0.51174 | 0.67007 | 4.02542 |
| 4 | 0.0125 | 1.50543 | 4 | 0.26103 | 0.66394 | 38.32835 |
| 8 | 0.00625 | 1.37768 | 2 | 0.09606 | 0.29571 | 1.14215 |

High-SNAIL means SNAIL > 0.8; it does not identify a functional branch tip.
RMS radius measures spread about the collective centroid in model length
units. Force-proxy inputs are the per-cell values actually supplied to the
final GRN update, not biological stress measurements or final-inner forces.

The primary 4-versus-8 comparison agrees closely early, then separates.
Mean per-cell position RMS difference is 0.02665 over (10,15] and 0.34018
over (15,20]; mean six-node GRN RMSE is 0.00498 and 0.04809 respectively.
At hour 20, position RMS difference is 0.47745 model units (0.63474 mean
cell radii), and six-node GRN RMSE is 0.18859. Final RMS spread differs by
about 8.5% relative to the 4-step result. Agreement is not yet established.

The final fractions of cells with GRN force input above 1 are 0.71094,
0.95703 and 0.24219 for 1, 4 and 8 steps. Mean clip(input-1,0,1), the
existing YAP force-driving term, is 0.61413, 0.94462 and 0.10750. Regulatory
states depend on their histories and other inputs; these endpoint measures
do not isolate a single causal pathway. The early high-SNAIL pulse survives
refinement: sampled peaks before hour 14 are 23, 21 and 27 cells.

## Interpretation and limits

Mechanical resolution changes the coupled regulatory outcomes even when
GRN, chemical and polarity clocks are held fixed. The finer arm has much
smaller late force inputs than the 4-step arm, but this does not prove that
eight steps are accurate or converged. Large active forces alone do not
prove numerical overshoot. No timestep or source-code change is adopted.

This is one checkpoint, one seed, no future GRN noise and arrested growth.
It does not establish behavior with division or larger populations, nor
biological validity or sustained branching. Holding inner polarity and
propulsion eligibility is an explicit splitting choice; outer-step errors
and the original feedback lag remain. Endpoint marker counts cannot turn
the earlier activation into a demonstrated numerical artifact.

A next proposed test is 8 versus 16 inner steps. A complementary diagnostic
could account for active work within each outer interval: with phenotype
and propulsion vector f held fixed, W = U_pair - sum(f_i dot x_i) has
overdamped derivative -sum(|F_pair + f_i|^2)/damping, away from the source's
near-zero separation guard. Check W within an interval and rebase it when
regulation or polarity changes; do not demand global monotonicity across
those changes. This diagnostic has NOT been implemented or run here.

## Artifacts

- `steps1`, `steps4`, `steps8`: settings, summary, timecourse and 21 snapshots each.
- `analysis/coupled_substeps.png`: timecourses and matched-cell comparisons.
- `analysis/final_summary.csv`, `pairwise_differences.csv`, `window_summary.csv`:
  numeric readouts; the latter two are inside `analysis` as well.
- Archived runner, report and test scripts alongside baseline/checkpoint helpers.
- `analysis/provenance.json`: report input hashes; `sha256.json`: archive hashes.
