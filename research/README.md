# Research notebook

Project: agent based, self-organized branching in cell collectives, using
vessel-like growth as the primary reference.

| Record | Question | Status |
| --- | --- | --- |
| [Initial baseline](../tools/BRANCHING_BASELINE.md) | What does the unchanged model do through hour 10? | Completed exploratory runs and sensitivity checks |
| [E002](experiments/2026-11-06-long-baseline/README.md) | Does activation emerge with a longer observation window? | Completed, 3 seeds through hour 12 |
| [E003](experiments/2026-11-06-trajectory-analysis/README.md) | Are activated cells peripheral and moving like leaders? | Completed descriptive reanalysis |
| [E004](experiments/2026-11-06-protrusion-ablation/README.md) | Does expansion depend on protrusive force after hour 10? | Protocol and results recorded in experiment folder |
| [E005](experiments/2026-11-06-growth-arrest/README.md) | What changes when growth/division, protrusion, or both stop? | Completed four-condition comparison, 3 matched seeds |
| [E006](experiments/2026-11-06-fixed-population-long/README.md) | Does more time yield sustained peripheral activity in a fixed population? | Completed to hour 20; exploratory timestep sensitivity found |
| [E007](experiments/2026-11-06-shared-state-timesteps/README.md) | Does timestep sensitivity persist from an identical state with no future noise? | Completed 3 timesteps; substantial late differences remain |
| [E008](experiments/2026-11-06-polarity-memory-scaling/README.md) | Does time-scaled polarity memory improve timestep agreement? | Completed; no consistent improvement from this scalar rescaling |
| [E009](experiments/2026-11-06-frozen-regulation/README.md) | Does timestep sensitivity persist with regulatory outputs frozen? | Completed; substantial mechanical disagreement persists |
| [E010](experiments/2026-11-06-passive-mechanics/README.md) | Does sensitivity persist with only adhesion/repulsion moving cells? | Completed; large-step passive potential increases observed |
| [E011](experiments/2026-11-06-passive-refinement/README.md) | Do finer passive steps improve consistency and agreement? | Completed; every-step potential check passes, individual arrangements remain sensitive |
| [E012](experiments/2026-11-06-coupled-substeps/README.md) | What changes with finer mechanics while biological update timing stays fixed? | Completed 1/4/8 substeps; early agreement improves, late regulatory differences remain |
| [E013](experiments/2026-11-07-coupled-refinement/README.md) | Do remaining differences shrink with 8 versus 16 mechanical substeps? | Completed; regulatory/spread agreement improves, matched-cell position agreement does not |
| [E014](experiments/2026-11-07-pair-specific-adhesion/README.md) | Does pair-specific adhesion alter sorting beyond reduced overall attraction? | Completed passive 3-seed assay; selective adhesion reduces mixing and changes shape, no branching established |
| [E015](experiments/2026-11-07-guidance-tip-competition/README.md) | Do guidance, tip competition, or both create connected sprouts? | Completed 24 diagnostic runs; contraction or detached cluster migration, no persistent connected extensions |
| [E016](experiments/2026-11-07-static-tip-selection/README.md) | Can the inhibition circuit pattern a fixed network? | Completed 36 runs; reference circuit patterns sparse sheet, contracted network has no high mobile cells |
| [E017](experiments/2026-11-07-loaded-tissue-mechanics/README.md) | Can tissue resist overlap and stay attached under imposed traction? | Completed 72 runs; spring junctions pass localized-load screen, all variants detach at high total load |
| [E018](experiments/2026-11-08-junction-sprouting/README.md) | Do reference signaling, guidance and material junctions combine into sprouts? | Completed 48 runs; low overlap and activity contrast, combined arm still detaches despite prepatterning |
| [E019](experiments/2026-11-08-traction-budget/README.md) | Do total traction and follower allocation separately affect attachment and extension? | Completed 60 runs; lower loads preserve junctions, activity weighting extends farther, no persistent narrow sprouts |

These records distinguish observations from interpretations and proposed
mechanisms. They are exploratory computational evidence, not biological
validation. Preserve settings, source and input hashes, seeds, analysis
definitions, failed hypotheses and limitations alongside each experiment.
Cell counts are not counts of independent statistical replicates.
