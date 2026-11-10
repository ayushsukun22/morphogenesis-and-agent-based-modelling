# E019: follower traction and matched force budgets

Date: 2026-11-08 (Europe/Berlin). Protocol recorded before execution.

Question: with E018 material mechanics fixed, does reducing total traction or redistributing the same total load preserve attachment and permit narrow extensions? Use reference contact-based inhibitory signaling and +x guidance in ALL arms. Seeded 36-cell sheet, 12 fixed/24 mobile cells, prepattern signaling for 20 units, then simulate 12 units. Three seeds (42/43/44), outer dt 0.05 and inner substeps 16/32. No cell is assigned leadership.

Factorial diagnostic: total force budget F=2.5, 7.5 or 15 model force units; allocation weights w_i=b+(1-b)D_i^4 for mobile cells, b=0, 0.25 or 1. Force is +x F*w_i/sum(w). Thus b=1 is uniform, b=0 strongly favors high activity, and b=0.25 adds a follower floor. b is a weight parameter, NOT the fraction of total force going to followers, and b=0 does not make every low-activity cell force-free. If all weights vanish, use uniform allocation and log the event. Reference arm repeats original E018 prepatterned BOTH force .15+2.35D without normalization. Total 60 runs. All settings are exploratory, not measured biological constants.

Global normalization is deliberately a nonlocal diagnostic intervention, not a proposed cell-autonomous biological mechanism. It controls TOTAL traction exactly despite regulatory changes. Matched-F comparisons isolate allocation rules as interventions, but evolving state feedback means their cell identities need not remain matched. Budget comparisons isolate force scale within each allocation rule. The reference does not have a fixed total budget.

Mechanical rule remains repulsion 40/damping 2/spring 4/initial material contacts <=1.2/rest length initial distance/rupture >1.6/no new bonds. Regulatory rule remains geometric contacts <=1.2 with original E015 inhibition, prescribed cue and prepattern policy. Hold allocated forces during each outer interval. No growth, ECM, lumen or blood flow.

Primary readouts: minimum material-parent attachment over every inner update, rupture count, overlap maximum, material-connected extension and geometric tissue attachment. Apply E018 narrow-extension screening separately to geometric and material paths, including three-sample ID-overlap persistence. Material-only links are not a continuous vessel wall. A low-force connected state with no extension is NOT successful sprouting. Keep all arms, no post-hoc parameter selection. Check force-budget errors/fallbacks, finite states, anchors, every-inner active-work potential, resolution differences and exact shared-array replay of E018 reference trajectories. Application code stays unchanged.

## Results

Completed 60 runs and validated 1,500 saved states. At the finer resolution,
final material-connected max x (mean over seeds, in nominal cell diameters
because diameter=1) is:

| Total force F | Activity weighted, b=0 | Follower floor, b=0.25 | Uniform, b=1 |
| --- | --- | --- | --- |
| 2.5 | 0.532 | 0.300 | 0.283 |
| 7.5 | 1.489 | 0.874 | 0.772 |
| 15 | 2.945 | 1.451 | 0 |

This is max x among cells with intact material paths to the fixed parent,
not path length, branch length, displacement from each cell's initial
position, or a guarantee of continuous cell-body contact. At F=15 the
averages conceal different failure outcomes and must be read with the
attachment results below.

At F=2.5 and 7.5, ALL allocations retain all 24 mobile cells through intact
material paths, with ZERO ruptures across all three seeds and BOTH inner
resolutions. Minimum material attachment stays 100% at every inner step.
At F=7.5 activity-weighted traction yields mean final extent 1.489 versus
0.772 for uniform traction (1.93x), at the SAME total applied force.
These effects cannot be explained solely by a lower total force budget.

At F=15, allocation affects which failure occurs:

| Allocation | Fine-step final attachment, seeds 42/43/44 | Broken edges, seeds 42/43/44 |
| --- | --- | --- |
| b=0 | 24/24, 23/24, 23/24 | 6, 3, 8 |
| b=0.25 | 0/24, 24/24, 24/24 | 6, 5, 5 |
| b=1 | 0/24, 0/24, 0/24 | 6, 6, 6 |

These attachment classifications repeat at both steps. Activity weighting
reduces catastrophic parent separation relative to uniform loading here,
but does NOT prevent junction damage or occasional single-cell escape.
The follower-floor arm is seed-sensitive, not a robust fix. At b=0 maximum
single-cell force reaches 3.724, 7.276 and 5.252 across fine-step seeds,
versus 0.625 for every cell under the uniform F=15 condition. Concentrating
load therefore trades broad stress for potentially large local loads.

The E018 reference exactly reproduces its earlier detached trajectories.
Its mean total force over the movement phase is 23.13-24.06 across seeds,
with peak total 26.36-27.90, above the three diagnostic budget values.
This supports testing lower load, but the reference also uses a different
allocation rule; fixed-budget factorial comparisons provide the cleaner
separation of total scale from allocation effects.

## Shape and continuity outcomes

NO condition produces a persistent narrow extension by either geometric
or material-path screening. The geometric screen detects no candidates
at all. The material screen detects a single candidate at time 3.5 in
seed 43, F=15/b=0, at both resolutions; it does not persist across three
saved observations. That run ultimately loses a mobile cell. Neither
large max x nor the transient candidate demonstrates a stable sprout.

Material and geometric continuity differ. At F=2.5 all fine-step endpoints
are geometrically attached at cutoff 1.2. At F=7.5/b=0 geometric attachment
is 0/24, 24/24 and 24/24 across seeds; at b=0.25 and b=1 it is 0/24 in
every seed despite intact material paths. Thus elastic links can stretch
outside the geometric contact range without rupturing. Calling every
material-connected configuration a continuous tissue would be incorrect.
The pre-existing signaling/junction range mismatch remains a model issue.

## Conclusions and limits

1. Lowering load from the reference regime to F<=7.5 prevents material
   rupture in the tested parameter/geometry/time window.
2. At matched total force, activity-weighted allocation produces more
   forward extension than uniform allocation. Load DISTRIBUTION matters.
3. High total load remains failure-prone; high selective force can tear
   individual cells out, whereas uniform traction can detach the whole
   mobile tissue. No high-budget allocation succeeds robustly on all
   integrity/shape criteria.
4. Preserving connections and extending a sheet is insufficient for
   branching. The next question is how local cellular rearrangement,
   junction remodeling or cell supply could support a narrow connected
   extension. Those mechanisms were not tested here; simply increasing
   force is not justified by these results.

This is an intervention assay with global knowledge of all mobile cells,
not a proposed biological force controller. Allocation includes detached
mobile cells; an escaped high-activity cell can receive a substantial
part of the global budget and change loading of the remaining tissue.
An attached remainder must therefore not be interpreted as biologically
regulated protection. No cell is manually selected as leader, but D**4
weighting is itself an uncalibrated hypothesis. No claim of physiological
VEGF, Notch, cadherin or traction parameters is made.

Further limits: three seeds, a planar small sheet, fixed parent positions,
prescribed cue, fixed radii/population, immutable rest lengths, irreversible
rupture without rebinding, no ECM, lumen, perfusion or proliferation. The
screening geometry and half-unit output cadence limit branch detection.

## Validation and artifacts

Four allocation tests plus nine prior coupling/mechanics tests passed
(13 total). Fixed-budget force errors never exceed 7.11e-15; no zero-weight
fallback occurs. Across 345,600 mechanical updates no active-work potential
increase exceeds tolerance. Maximum overlap is 0.08981, inherited from the
initial geometry. All saved states are finite and anchors remain exactly
fixed. Initial regulatory states and geometry match between arms. The six
E018 reference runs match all shared saved arrays EXACTLY, including
positions, signaling, junctions and times. Application-source hashes match.

Maximum final matched-cell position RMS difference between 16/32 substeps
is 0.00903 model units; activity RMSE is at most 0.00330. This improves
numerical agreement in these comparisons but is not proof of universal
convergence. Findings are exploratory computational evidence.

```powershell
rtk proxy python tools/test_branching_traction_budget.py
rtk proxy python tools/test_branching_junction_sprouting.py
rtk proxy python tools/test_branching_prerequisites.py
rtk proxy python tools/branching_traction_budget.py --output research/experiments/2026-11-08-traction-budget
rtk proxy python tools/branching_traction_budget_report.py research/experiments/2026-11-08-traction-budget
```

Use a fresh output folder for reruns. Settings preserve runner/helper/source
hashes; archived scripts accompany all 60 NPZ histories, timecourse/final/
condition summaries and resolution differences. `validation.json` records
exact-reference input hashes and trajectory hashes; `sha256.json` covers
the archive. Saved forces are those held during the preceding outer
interval, with zeros at the initial observation before movement.

`traction_budget_summary.png` shows finer-step seed means, not uncertainty
estimates; consult individual CSVs for seed failures. In the seed-42 state
figure, blue segments are intact junctions, red crosses fixed cells and
colors abstract tip activity. Panels have independent axis limits. Model
time begins when movement starts, excluding 20 units of prepatterning.
No production code or parameter set has been adopted.
