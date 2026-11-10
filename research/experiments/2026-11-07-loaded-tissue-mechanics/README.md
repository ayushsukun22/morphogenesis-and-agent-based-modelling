# E017: loaded tissue mechanics

Date: 2026-11-07 (Europe/Berlin). Protocol recorded before execution.

Separate mechanical ability to maintain a parent connection under imposed load from any signaling. Use E015 36-cell initial planar sheet, 12 fixed parent cells, 24 mobile cells and seeds 42/43/44. No regulation, division, growth, ECM or state-dependent forces. Compare original repulsion 4/cutoff 2.5; repulsion 40/cutoff 2.5; repulsion 40/short cutoff 1.2; and explicit initial material junctions with repulsion 40. Pair attraction is 0.49 in the first three variants. Damping remains 2.

Material junctions are an abstract alternative, not measured VE-cadherin: place springs on initial contacts <=1.2 with initial distance as rest length, spring constant 4, and irreversible removal when distance>1.6 after an integration step. All overlapping pairs repel. No uniform pair attraction or new junction formation in this variant. This changes more than one mechanism and cannot isolate spring effects from removal of broad attraction.

Loads: none; one front cell closest to y=0 with +x force 2.5; all mobile cells each with +x force 2.5. The singled-out load site is imposed for a mechanical test, NOT an emergent leader. Uniform load therefore has 24 times total applied force and is not a matched-total-force localization comparison. Run 5 model time units at dt 0.005/0.0025. All variants use identical starting arrays and load sites.

Primary readouts: geometric parent attachment at cutoff 1.2, alternative attachment at 1.6, maximum overlap over run, front displacement, broken junctions and timestep sensitivity. Record no-load stability as well as traction response. Monitor fixed-load W=U-sum(f dot x) every step, count increases before junction deletion, then rebase after energy-reducing deletions. Finite states and exact fixed anchors required. Operational screening threshold for acceptable overlap is <=0.1 model units (10% nominal diameter), not biological validation. No branch claim from imposed extension. Application code unchanged.

## Results

Completed 72 runs, saving 1,512 position states. Under the one-front-cell
load, at dt=0.0025 and time 5:

| Mechanics | Final parent-attached mobile cells | Mean loaded-cell displacement | Worst overlap over run, across seeds | Attachment/overlap screen |
| --- | --- | --- | --- | --- |
| Original | 24/24 | -3.315 | 0.873 | Fail |
| Repulsion x10 | 24/24 | +0.400 | 0.251 | Fail |
| Repulsion x10, short attraction range | 23/24 | +6.085 | 0.090 | Fail: pulled cell escapes |
| Material junction springs + repulsion x10 | 24/24 | +1.165 | 0.090 | Pass |

The material-junction loaded-cell displacement spans 1.159-1.171 across
seeds. Mean mobile-cell displacement is 0.194, so the entire population
does not simply translate by the loaded-cell displacement. No material
junction breaks under the single-cell load; all cells remain geometrically
attached at saved observations using both cutoff 1.2 and 1.6. The maximum
overlap in these arms is already present in the jittered initial geometry,
not increased by loading beyond that initial maximum.

No-load controls matter: original mechanics contracts by mean loaded-site
displacement -3.553 with severe overlap (worst 0.877); stiffer repulsion
alone reduces this to -0.825 but still exceeds the overlap screen. Short
range and material junctions preserve the sheet without severe overlap,
with mean site displacement -0.080 and +0.037 respectively. Thus low
overlap alone does not establish the ability to transmit pulling forces.

At the higher TOTAL load of 2.5 per mobile cell, all four variants detach
all 24 mobile cells. The spring arm breaks six initial junctions per seed,
disconnecting from the fixed parent. This does not isolate load localization
from total load: that would require a matched-total-force comparison.

The original no-load and single-load arms briefly lose geometric attachment
at saved intermediate times before rejoining the fixed-parent region. Their
final 24/24 attachment must not be mistaken for persistent tissue integrity.
The screen requires full attachment at every saved time and overlap <=0.1,
not merely a connected endpoint. Between saved times, geometric attachment
is not independently measured; intact spring topology is supported by the
zero break count in the successful arm.

## Numerical checks and limitations

All states were finite and fixed anchors stayed exactly fixed. Across
108,000 integration updates, no fixed-load work-adjusted potential increase
exceeded tolerance. Maximum final position RMS disagreement between the
two tested steps is 0.02152 (original), 0.00439 (stiffer), 0.000175 (short
range) and 0.000598 (material junctions). Qualitative success/failure
patterns repeat at both steps. This is not proof of convergence at all
loads or populations. Five prerequisite tests passed; see E016 for the
recorded floating-roundoff correction in the reference-force test.

The successful arm is a candidate mechanical model, not a measured or
validated VE-cadherin law. It changes the force description and interaction
graph together; a single spring coefficient and break distance were tested.
Junctions are initialized, retain rest lengths, break irreversibly and do
not reform. This may constrain rearrangement and cannot establish that
sprouting will occur when coupled to regulatory selection. The 10% overlap
screen is an explicit modeling tolerance, not a biological standard.
The loaded cell is selected by us, so its displacement is NOT evidence
of self-organized tip selection or a branch. No growth or ECM is present.

## Combined interpretation with E016

Two prerequisites now have partial support: the reference signaling rule
can form persistent high/low patterns on an uncontracted contact network,
and spring junctions with stiffer repulsion can transmit one localized
load without collapse or detachment under this assay. Neither demonstrates
a complete sprout. A next coupled test could use these mechanics with the
REFERENCE signaling circuit and guidance, preserving factorial controls.
It should track mobile leader density, evolving contact/junction graphs,
connectivity, overlap and persistent extensions. It must allow a negative
result: the static pattern has many high cells and the material graph may
restrict elongation. This follow-up has not been run or adopted in the app.

## Reproduction and artifacts

```powershell
rtk proxy python tools/test_branching_prerequisites.py
rtk proxy python tools/branching_prerequisites.py mechanics --output research/experiments/2026-11-07-loaded-tissue-mechanics
rtk proxy python tools/branching_prerequisites_report.py mechanics research/experiments/2026-11-07-loaded-tissue-mechanics
```

Use a fresh output folder to rerun. Settings, archived scripts, initial and
trajectory NPZ files, timecourse/final/condition CSVs, resolution differences,
validation JSON and SHA-256 archive manifest preserve the experiment.
`mechanics_states.png` shows finer-step seed-42 endpoints (mobile cells blue,
fixed cells red). Panels have independent axis limits. Application code
remains unchanged.
