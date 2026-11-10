# E018: recombined junction-based sprouting assay

Date: 2026-11-08 (Europe/Berlin). Protocol recorded before execution.

Question: do the reference inhibitory signaling circuit, guidance and E017 material junction mechanics together generate persistent connected narrow extensions without a preselected leader? Reuse E015 initial 36-cell sheet (12 fixed parent, 24 mobile), seeds 42/43/44. Four factorial conditions: neither guidance nor inhibition, guidance alone, inhibition alone, both. The cue magnitude affects signaling in all arms; guidance changes traction direction only.

Use E015 regulate/traction functions unchanged: prescribed sigmoid cue increasing in x; D/N relaxation tau 0.5; inhibitory Hill exponent 4 and half inhibition 0.3; active magnitude 0.15+2.35D; guided direction +x or unguided direction opposite local unit-neighbor sum. All cells can change activity; no tip location is assigned. Forces on fixed parent cells are suppressed. No growth, ECM, lumen or flow.

Use E017 material mechanics: all-pair overlap repulsion 40, damping 2, initial contacts <=1.2 receive spring constant 4/rest length equal initial distance, irreversible rupture above 1.6, no new junctions or broad pair attraction. Contact-based signaling still uses current geometric distances <=1.2, not spring connectivity. This distinction is explicit: junctions may persist beyond the signaling contact range.

Compare two timings. Immediate: same seeded D around 0.5 and N=0, signaling/motion begin together. Prepatterned: run the SAME condition-specific signaling on the stationary initial positions for 20 model units at dt 0.05 before time-zero migration. No mechanical relaxation occurs during prepatterning. Thus initial regulatory states differ deliberately between timing/competition arms; geometry and seed are matched. The prepatterned no-inhibition control also receives this waiting period. This is a sequencing test, not a new circuit.

For each timing run 12 time units of coupled dynamics at outer dt 0.05 with 16/32 mechanical substeps. Total: 48 runs. Hold traction within each outer interval; recompute elastic forces and check rupture each inner step. Monitor active-work-adjusted potential before rupture/rebase afterward, finite states, fixed anchors, maximum overlap and source hashes. No global energy-decrease requirement across outer updates.

Measure geometric parent attachment, material-junction attachment, broken junctions, mobile tip levels/heterogeneity/isolated high cells and E015 narrow-extension proxy. Candidates require >=3 parent-attached cells with x>0, max x>=2 and y-span<=2; persistence is ID overlap across three half-unit observations (one time unit). This is NOT a validated skeleton/bifurcation detector. Compare all arms, seeds and resolutions; report negative results and inherited 10%-diameter overlap screen. No application edits or automatic production adoption.

## Results

Completed 48 runs and validated 1,200 saved observations. At the finer
mechanical resolution (32 substeps), model time 12:

| Timing | Condition | Parent attachment through material junctions | High-activity mobile cells, seed range | Broken junctions, seed range | Persistent geometric narrow-extension windows |
| --- | --- | --- | --- | --- | --- |
| Immediate | Neither | 100%, all seeds | 16-17 | 5-8 | 0 |
| Immediate | Guidance | 0%, all seeds | 24 | 6 | 0 |
| Immediate | Competition | 100%, all seeds | 7-13 | 0-2 | 0 |
| Immediate | Both | 0%, all seeds | 10-12 | 6 | 0 |
| Prepatterned | Neither | 0% in seed 42; 100% in seeds 43/44 | 17-18 | 8-9 | 0 |
| Prepatterned | Guidance | 0%, all seeds | 24 | 6 | 0 |
| Prepatterned | Competition | 100%, all seeds | 6-12 | 1-2 | 0 |
| Prepatterned | Both | 0%, all seeds | 10-11 | 6 | 0 |

Attachment here describes paths through intact material junctions, not
cell-body contact or a continuous vessel wall. Geometric attachment at
distance 1.2 is much lower in the unguided arms as cells spread apart while
some spring paths remain. High activity alone is not functional leadership.

For BOTH, the six junctions crossing from fixed parent to mobile tissue
break in every seed at BOTH mechanical resolutions and BOTH timings.
At finer resolution the first saved completely disconnected state is time
5.0 in all six timing/seed combinations; exact rupture timing within the
half-unit output interval is not claimed. Prepatterning gives 1-3 high
mobile cells initially, versus zero initially in the immediate arm, but
does not rescue parent attachment.

By time 12, the combined arm has pronounced activity contrast: mean
within-mobile activity SD across seeds is 0.3828 immediate and 0.3767
prepatterned. This differs from the near-uniform E015 mobile state, but
is observed in already detached tissue. There are about 10 isolated high
cells among 24 mobile cells, not a single sparse migrating tip. Guidance
alone drives all 24 mobile cells high. The effect of inhibition is visible,
but does not establish a controlled branch-forming leader/follower system.

The predefined geometric sprout proxy finds zero candidates at every saved
time in every run. A supplemental material-path version of the SAME size,
length and width thresholds also finds zero candidates in BOTH. Before
detachment, peak sampled material-attached max x is 1.608 immediate and
1.923 prepatterned (finer resolution, max across seeds), below the length-2
candidate threshold. This shows that the negative BOTH result is not
solely caused by using geometric attachment instead of spring paths.

The supplemental material-only screen does flag one transient candidate
at times 8.0, 8.5 and 9.0 in the prepatterned NEITHER control, seed 42, at
each resolution. This was discovered during post-run diagnostic analysis,
not a predeclared success criterion. It is not detected by the geometric
screen; the control later loses all parent material attachment. It is not
evidence that guidance plus competition produces a stable connected vessel.
Candidate-ID persistence was not separately computed for this supplemental
screen, and no bifurcation/skeleton measure was evaluated.

## Mechanical and numerical checks

Maximum overlap is 0.08981 model units (8.98% nominal diameter), inherited
from the initial jittered geometry; the new mechanics avoids E015's severe
compaction in these trials. No active-work-adjusted potential increase
exceeds tolerance across 276,480 inner steps. Each run has 240 regulatory
updates during movement; prepatterned runs additionally have 400 stationary
regulatory updates. No global monotonic energy requirement is imposed.

Four new coupling tests and ten existing signaling/mechanics tests passed
(14 total). Arrays are finite; initial states match the prescribed seed,
timing and condition policy exactly; anchors remain exactly fixed; junction
graphs stay symmetric and only lose edges. Application-source hashes match.

Across all arms, final position RMS disagreement between 16/32 substeps
reaches 0.45596 model units and activity RMSE reaches 0.12638. Therefore
not every detailed trajectory is resolution-independent. The BOTH outcome
of complete parent detachment and no persistent geometric extensions
repeats across all seeds and resolutions; neither this repetition nor a
passing energy check validates the biological force/signaling model.

## Interpretation and next question

The prerequisite successes do not automatically combine into sprouting.
This recombination achieves low overlap and activity contrast, but not a
parent-connected narrow extension. Prepatterning alone is insufficient.
It would be inaccurate to conclude either that inhibition has no effect
or that biological VEGF/Notch cannot produce sprouts.

The traction mapping still gives every mobile cell a positive drive, with
strength 0.15+2.35D; this is not the isolated single-cell load that passed
E017. Elastic junctions have finite strength, no renewal and no new-cell
supply. A next controlled test could vary basal/follower traction and
total active load while keeping junction mechanics fixed, with a
matched-total-force localization control. This would distinguish load
capacity from load distribution rather than silently strengthening springs.
Growth and junction remodeling are separate future hypotheses, not fixes
demonstrated by these results.

An additional modeling issue is that signaling uses geometric contacts
<=1.2 whereas material junctions survive to 1.6. At the combined finer-step
endpoints, 4-6 of the remaining 54 material edges lie outside signaling
range. All mobile cells nevertheless have at least one signaling neighbor
at those endpoints. In unguided controls, a high cell can be marked
isolated simply because geometric neighbors are lost, even without any
inhibitory circuit. Thus counts of isolated high activity must not be
interpreted as evidence of lateral inhibition on their own. A junction-based
signaling graph is another distinct, untested model choice.

Limits: prescribed cue, planar centers, fixed parent boundary, 3 seeds,
one parameter set, no growth/division, no explicit ECM or lumen/flow,
irreversible material edges, uncalibrated regulatory dynamics and forces.
No production-model change or branch claim is adopted.

## Reproduction and artifacts

```powershell
rtk proxy python tools/test_branching_junction_sprouting.py
rtk proxy python tools/test_branching_prerequisites.py
rtk proxy python tools/test_branching_sprouting_assay.py
rtk proxy python tools/branching_junction_sprouting.py --output research/experiments/2026-11-08-junction-sprouting
rtk proxy python tools/branching_junction_sprouting_report.py research/experiments/2026-11-08-junction-sprouting
```

Choose a fresh output folder to rerun. Settings preserve source, runner and
helper hashes; archived scripts accompany trajectory NPZ arrays and CSVs.
`network_diagnostics.csv` includes the supplemental material-path screen;
`validation.json` and `sha256.json` preserve checks and provenance. In the
seed-42 state figure, gray tracks show motion, blue segments intact final
junctions, red crosses fixed parent cells and color abstract tip activity.
Panels have independent axis limits. Timecourse bands show seed min/max,
not confidence intervals. Times are model units since movement begins,
excluding the separately recorded stationary prepatterning period.
