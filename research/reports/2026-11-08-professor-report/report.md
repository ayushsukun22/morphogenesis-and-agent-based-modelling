# Self organized branching in a cell based simulation

Research progress report for supervisory discussion

Project: embryo_twin.digital.kit | Reporting date: 8 November 2026

## Main finding and purpose of this meeting

The project investigates whether vessel-like tissue branches can emerge from interacting cell agents without assigning branch locations, branch points or tip-cell identities in advance. The current work has demonstrated growth, regulatory activation, cell sorting, signaling contrast and transmission of localized pulling forces. It has not demonstrated robust, persistent branching.

The experimental record is useful because it separates these ingredients. In the original growing model, delayed SNAIL activation did not establish functional tip behavior. Numerical diagnostics revealed substantial sensitivity to mechanical stepping. Later synthetic assays showed that selective adhesion sorts cells and that explicit junctions can transmit a localized load. Recombining guidance, inhibition and junctions still caused parent separation. Lower force budgets preserved junctions and activity-weighted pulling increased extension, but did not yield persistent narrow sprouts.

The decision requested from the supervisor is how to turn these findings into a focused research question. In particular, should the next stage investigate cell rearrangement and junction remodeling, growth and cell supply, or a different cell representation? The requirement remains a cell-based, self-organizing model; explicitly programming branches to split would answer a different question.

## How to read the report

Sections 1-5 explain the project and executable rules. Section 6 defines the measurements. Sections 7-8 catalogue every experimental record. Sections 9-17 present methods, quantitative results and limits. Sections 18-20 distinguish supported conclusions from hypotheses, identify decisions for the meeting, and explain reproducibility. The final appendices locate the full notebooks and biological references.

Evidence status: exploratory computational experiments, primarily three initialization seeds; no experimental biological validation or inferential statistical claims. The original application source was not changed by the diagnostic experiments. Later assays are separate models in research scripts, not implemented application upgrades.

<!-- PAGE -->
# 1 Research aim and biological motivation

## What the project seeks to explain

The long-term aim is to understand how local cell interactions produce connected, branching tissue. Vessel-like branching is the main reference, while reproducing every difference between blood vessels and nerves is outside the immediate scope. A useful first target is a persistent cellular sprout attached to a parent tissue; a perfused blood vessel is a much later target.

An agent is an individual cell with position, size, lineage and internal regulatory state. Self-organization means that tissue-scale geometry emerges from interactions, environmental sensing, movement and population dynamics. Explicit local rules are necessary in any agent-based model; the research issue is whether those rules prescribe the branching pattern itself. In this project, branch coordinates and splitting events should not be assigned directly.

## Why angiogenesis is relevant

Experimental studies support VEGF-guided tip migration and Dll4-Notch restriction of neighboring tip behavior [B1, B2]. Tip and stalk roles are dynamic: cells can compete and exchange positions [B3]. Regulated VE-cadherin junction dynamics contribute to endothelial rearrangement [B4]. These findings motivate hypotheses about guidance, local competition and cohesion; they do not validate the equations or coefficients used here.

The original application is not yet a mechanistic angiogenesis model. Its morphogen is generic, its neighbor input is mean neighboring SNAIL, and its polarity is based on neighbor geometry. The later VEGF-like assays introduce different, explicitly documented abstract rules. They do not contain measured VEGF concentrations, explicit filopodia, integrin traction, ECM degradation, lumen formation or blood flow.

## What would count as progress

A stronger result would combine localized leader behavior, continuous attachment, narrow advancing tissue and persistence across seeds and numerical resolutions. A branch point would additionally require a topology or geometry measure capable of distinguishing bifurcation from a stretched sheet. High marker values, elongation or a large displacement are intermediate observations, not sufficient evidence of branching.

Model time and length are computational units. Earlier scripts label time as model hours, but no mapping to measured developmental hours or micrometers has been calibrated. The phrase digital twin is the project name and aspiration; the present system is a research prototype rather than a validated twin of a specific biological specimen.

<!-- PAGE -->
# 2 Software architecture and operation

## Application and model layers

The NVIDIA Omniverse Kit application embryo_twin.digital.kit enables the embryo_twin.sim extension. The extension constructs an EmbryoSimulation, exposes Start, Stop, Reset and simulation-speed controls, and synchronizes simulated cells to USD spheres. Sphere color changes from blue toward red as SNAIL rises. This is a visualization of an internal state, not a validated vessel-tip classification.

| File in embryo_twin.sim | Responsibility |
| --- | --- |
| simulation.py | Cell registry, environment sampling, update ordering, growth and division gates, fate labels and field evolution |
| cell.py | Cell state, phenotype synchronization, radius growth, daughter construction and lineage |
| grn.py | Six regulatory nodes and their mechanical phenotype readouts |
| mechanics.py | Pair forces, neighbor signals, polarity, protrusion and overdamped movement |
| chemical_field.py | Grid sampling, consumption, deposition, diffusion and decay |
| extension.py | Kit UI, event callback, USD geometry, metadata and rendering synchronization |

The mechanics are custom NumPy calculations. The presence of PhysX-related Kit dependencies does not mean that the cell force law is a PhysX collision model. The research runners load the pure Python model without launching Kit, enabling fixed timesteps, controlled interventions and data export.

## One application update in execution order

1. Advance model time and obtain the current living-cell list.
2. For each cell, sample oxygen and morphogen, copy the previous mechanical force proxy into the GRN, integrate the GRN and synchronize its phenotype.
3. If oxygen permits, increase radius and age. Consume oxygen, deposit morphogen and update the discrete fate label.
4. If division criteria are met, construct daughters and mark the parent dead. The current-cell mechanics list is processed before daughters are added.
5. Compute pair forces, update neighbor signals and polarity, apply protrusion, move cells and record the force proxy for the next GRN update.
6. Add daughters, advance both chemical grids and prune dead parents.

The interactive timestep is min(real-frame dt, 0.1) multiplied by the speed setting. Thus changing frame rate or speed can change numerical behavior. Headless results at fixed dt are not automatically equivalent to an interactive demonstration. The UI stops at the 1,500-cell cap; the diagnostic runner also records stopping conditions.

<!-- PAGE -->
# 3 Cell state growth fields and fate labels

## Cell initialization and population rules

The application starts with two cells, radius 0.5, with independent position coordinates sampled from a normal distribution with standard deviation 0.65. Reset uses random seed 42. Each cell stores a unique ID, position, velocity, polarity, radius, age, parent ID, generation, lineage string, oxygen, morphogen, fate and GRN. Daughter cells inherit deep copies of the parent's GRN and lineage information.

Radius growth is dr/dt = 0.18 times the cell-cycle readout, provided oxygen > 0.2. Age advances when that growth method is called. A cell divides when radius >= 0.9, oxygen > 0.2 and the population gate permits. Each daughter radius is parent radius divided by 2^(1/3), approximating equal total volume. Daughters are placed on opposite sides along a normalized mixture of parental polarity and a random unit vector, with offset 1.15 times daughter radius. Parent removal follows daughter creation.

## Environmental fields

Both fields are 32 x 32 x 32 grids spanning coordinates -10 to 10 in each dimension. Oxygen starts at 1, with grid diffusion coefficient 0.8 and no decay. Morphogen starts at 0, with coefficient 0.3 and decay 0.05. Cells remove oxygen at 0.05 times radius times dt and deposit morphogen at 0.02 times dt. These amounts are distributed over eight nearby voxels using trilinear weights; sampling uses the corresponding trilinear interpolation.

The explicit field update uses a six-neighbor grid Laplacian with edge padding for zero-flux boundaries. Its code does not divide that Laplacian by physical grid spacing squared, so these coefficients must not be reported as calibrated physical diffusivities. Concentrations are clipped at zero. There is no configured oxygen replenishment in the application update, although the field class supports optional sources.

Cell positions are not confined by the chemical grid. Out-of-grid sampling clamps coordinates to the edge. Three centers left the grid in the seed-42 E002 endpoint, limiting interpretation of that run's local chemical exposure.

## Two different classification systems

The application assigns TYPE_A if morphogen > 0.75, otherwise TYPE_B if oxygen < 0.30, otherwise PROGENITOR, after a minimum time-in-fate gate of 0.6. These labels do not directly determine adhesion. A separate GRN classify_fate method returns MESENCHYMAL, PROGENITOR or TRANSITIONAL from SNAIL/SOX2; simulation.update_fate does not use that method. These labels must not be conflated with each other or with endothelial tip and stalk cells.

<!-- PAGE -->
# 4 Regulatory equations and phenotype readout

## Implemented regulatory drift rules

Let O be sampled oxygen, M sampled morphogen, P the previous net-force magnitude, and Q the previous neighbor signal. Let H, C, S, T, Y and N denote HIF1a, cell cycle, SOX2, SNAIL, YAP and Notch. The code applies the following drift expressions multiplied by dt:

| Node | Drift expression |
| --- | --- |
| HIF1a H | (1 - O) - 0.4 H |
| Cell cycle C | 0.4 O - 0.3 H - 0.4 C |
| SOX2 S | 0.3 + 0.6 M - 0.9 T - 0.3 S |
| SNAIL T | 1.2 Y - 0.4 S - 0.6 N - 0.4 T + 0.05 (u - 0.5) |
| YAP Y | 0.5 clip(P - 1, 0, 1) - 0.5 Y |
| Notch N | Q - 0.4 N |

Here u is a uniform random draw. Updates occur sequentially in the table's order: later expressions can see nodes already updated earlier in the same call. Nodes are then clipped to [0, 2]. The noise is multiplied by dt as part of the drift; it is not a calibrated stochastic differential equation with square-root-dt forcing. Changing dt also changes the number of random draws.

Initial values are S=1, T=0, Y=0, N=0, C=0.5 and H=0. The regulatory names express intended roles; their coefficients and interactions are exploratory hypotheses, not fitted biochemical kinetics.

## Mechanical phenotype

Define sigmoid_k,z(a) = 1 / (1 + exp(-k(a-z))). After each GRN update:

- Adhesion A = sigmoid_5,0.5(SOX2).
- Protrusion strength B = sigmoid_5,0.4(0.6 SNAIL + 0.4 YAP).
- Polarity bias q = sigmoid_5,0.3(SNAIL).
- Division-rate readout = clip(cell cycle, 0, 2).

cell.sync_phenotype copies adhesion and protrusion to the fields used by mechanics. Despite the configuration's nominal global adhesion 1.8 and protrusion 0.12, the current force implementation uses these per-cell values. Fate labels such as TYPE_A do not introduce a separate binding-compatibility rule.

The neighbor input Q is mean neighboring SNAIL within the mechanics cutoff, not explicitly Dll4 ligand binding to a Notch receptor. P is a net-force norm, not contact pressure: balanced forces can cancel even under substantial compression. Consequently a high YAP/SNAIL value is a model state, not proof of biological stress or functional tip behavior.

<!-- PAGE -->
# 5 Mechanical forces movement and feedback

## Pair interactions in the original application

For cells i and j, let d be center distance, e the unit direction from i to j and R = 2.5 the interaction cutoff. Pairs with d < 1e-8 or d > R are skipped. When d < r_i+r_j, the overlap produces equal-and-opposite repulsion of magnitude 4(r_i+r_j-d). Adhesion produces attraction of magnitude A_i A_j min(1, R/d).

Because attraction is evaluated only at d <= R, the factor min(1, R/d) is 1 for every evaluated pair. Thus attraction does not taper within the interaction range; it stops abruptly at the cutoff. It also acts on overlapping pairs. Multiple attractive neighbors can therefore compress a cell population strongly against a soft repulsive law.

The same force law applies to progenitor-progenitor, TYPE_A-progenitor and other pairs. Strengths can differ indirectly through SOX2-derived A values, but the model contains no explicit surface-protein repertoire or pair-compatibility matrix. E014 tested that missing distinction in a separate assay.

## Polarity and active motion

Each cell sums unit directions toward neighbors within 2.5. Its proposed polarity points opposite that sum, scaled by its GRN polarity bias. With memory m=0.2, it blends m times previous polarity plus (1-m) times the proposal, then normalizes the vector. Cells without neighbors or with a directionally balanced sum do not receive protrusion in that update. This is geometric outward steering, not chemical-gradient sensing.

Normalization means polarity bias is not simply a multiplicative reduction of the final propulsion magnitude. It influences directional blending. The active force is B_i times the resulting unit polarity for eligible cells. Total pair and active forces give velocity v_i = F_i / 2 and position x_i <- x_i + dt v_i. There is no inertia or resolved cell shape.

The final force norm is saved as cell.pressure. Neighbor SNAIL averages are also buffered for later regulatory updates. Feedback is therefore lagged, and the fixed memory coefficient acts per update rather than per physical hour. These implementation choices motivated E007-E013.

## Numerical consistency diagnostic

For frozen passive interactions, a pair potential U can be constructed whose negative gradient equals the pair force. Continuous overdamped motion should dissipate U away from skipped zero-distance configurations. During an interval with fixed active force f, W = U - sum(f_i dot x_i) should dissipate instead. These are computational integration checks, not measured biological energies; active regulation changes the reference between outer updates.

<!-- PAGE -->
# 6 Measurements controls and interpretation

| Measurement | Definition and interpretation |
| --- | --- |
| Random seed | Integer initializing the pseudorandom sequence. Seeds 42, 43 and 44 vary initialization and subsequent draws; their numerical values are labels, not biological quantities. |
| RMS radius | sqrt(mean_i norm(x_i-centroid)^2). Describes population spread, not an individual cell radius, outer tissue boundary or branch length. |
| Elongation | A ratio based on spatial covariance eigenvalues. Near 1 is roughly isotropic. Baseline reports regularize by cell size; later planar assays use two-dimensional covariance. Values from these definitions are not interchangeable. |
| High SNAIL | SNAIL > 0.8 in the original-model studies. This is a marker threshold, not a validated tip classifier. |
| Geometric connectivity | Baseline graph edges use d <= 1.25(r_i+r_j); later radius-0.5 assays usually use d <= 1.2. Component counts depend on the threshold. |
| Material connectivity | A path through explicitly retained spring junctions. A stretched spring can remain even when cell centers no longer satisfy geometric contact. |
| Position RMS difference | RMS Euclidean separation of matching cell IDs between resolutions, without rigid alignment. Distinct from within-population RMS radius. |
| GRN RMSE | Root mean squared difference across matched cells and regulatory nodes. Original six-node GRNs and later D/N circuits are different state spaces. |
| Unlike contact fraction | Fraction of geometric contact edges joining distinct fixed adhesion identities in E014. Lower values indicate sorting only when interpreted with connectivity. |

## Extension screening in the synthetic sprouting assays

A candidate is a connected component of parent-attached cells with x>0, at least three cells, max x>=2, and y-span<=2. Persistence requires overlapping candidate IDs across three half-unit observations, covering one model time unit; consecutive overlaps must include at least half the smaller ID set. Some reports additionally evaluate the same screen using material paths. This is a screening proxy, not a validated vessel skeleton or bifurcation measure.

Length 2, width 2, activity thresholds and the 10%-diameter overlap tolerance are explicit operational choices, not biological standards. A negative screen cannot exclude all brief or differently shaped protrusions. Conversely, passing one screen does not demonstrate a perfused or even continuous vessel.

## Experimental design principles

Seeds are exploratory realizations. Interacting cells and related daughters are not independent replicates. No significance tests or confidence intervals are claimed. Later multi-condition comparisons reuse three seeds; increasing the number of parameter combinations does not create an equivalent number of independent biological samples. Comparisons use matched initial states, logged interventions and timestep checks where available. Complete checkpoints, selected observations and cached plotting data are distinguished.

<!-- PAGE -->
# 7 Experiment catalogue for the original model

The initial baseline and E002-E013 concern the original model or explicit diagnostic reductions of it. E003 reanalyzes observations; it is not a new simulation experiment. E007 onward uses a complete seed-44 hour-10 checkpoint for the main numerical comparisons.

| Record | Question and intervention | Main result |
| --- | --- | --- |
| Initial baseline | Existing growing model to hour 10; seeds 42-44, dt=0.05 | 256 cells per seed, no high SNAIL; connectivity sensitive to proximity threshold and dt. |
| E002 | Extend unchanged growing model to hour 12 | Delayed high-SNAIL activation, roughly rounded collectives; no demonstrated sustained branches. |
| E003 | Analyze E002 location, lineage motion and neighbor retention | Activated groups often more internal; only five start-high cells available for prospective motion analysis. |
| E004 | Turn protrusive force off from hour 10 to 12 | Outward-motion summary falls, but residual expansion remains. |
| E005 | Factorial growth/division and protrusion arrest at hour 10 | Growth strongly supports spread; protrusion contributes motion without proving branches. |
| E006 | Growth-arrested populations to hour 20; added dt sensitivity run | Contraction/clustering and variable activation; confounded dt repeat motivates an identical-state test. |
| E007 | Same full checkpoint, no future GRN noise, three timesteps | Large late differences persist; not explained solely by different initial tissues or random draws. |
| E008 | Scale polarity memory with dt | No consistent improvement; not adopted. |
| E009 | Freeze regulatory nodes and outputs | Mechanical timestep sensitivity persists without evolving regulation. |
| E010 | Additionally disable protrusion | Passive large-step potential increases indicate non-dissipative finite-step behavior in this test. |
| E011 | Refine passive mechanics and monitor every step | Potential checks pass at finer steps; bulk spread becomes close, individual arrangements remain sensitive. |
| E012 | Fixed outer dt=0.05, 1/4/8 mechanical substeps | Mechanical resolution changes regulatory outcomes despite equal biological clocks. |
| E013 | Compare 8/16 substeps using the same protocol | Some aggregate/regulatory agreement improves; position disagreement does not uniformly shrink. |

No numerical experiment establishes a universally converged production timestep. Frozen regulation, arrested growth, zero future noise and operator splitting are diagnostic interventions, not the unchanged biological model. Exact replay checks establish implementation fidelity for the comparisons tested, not biological accuracy.

<!-- PAGE -->
# 8 Experiment catalogue for synthetic assays

These assays introduce new model assumptions in separate tools. Their results should not be presented as if they occurred in the original Kit application. Most remove growth and division so that a specific mechanism can be isolated.

| Record | Design and number of runs | Main result |
| --- | --- | --- |
| E014 | Passive mixed aggregate, 64 cells; four adhesion conditions, three seeds, two steps; 24 runs | Pair compatibility produces sorting and anisotropy beyond a uniform strength reduction; no branching established. |
| E015 | 36-cell planar sheet with fixed parent; neither/guidance/competition/both, three seeds, two steps; 24 runs | Unguided tissue contracts; guided tissue detaches. Severe overlap limits biological interpretation. |
| E016 | Static initial or contracted network; no inhibition/reference/steeper circuit, three seeds, two steps; 36 runs | Reference circuit patterns the sparse sheet; contracted-network high cells occur only at former fixed-parent sites. |
| E017 | Four mechanics variants, no/single-cell/all-mobile load, three seeds, two steps; 72 runs | Spring junctions plus stiff repulsion transmit one imposed pull without gross overlap; large total load still detaches tissue. |
| E018 | Recombine guidance/inhibition/junctions; immediate or prepatterned movement, three seeds, two steps; 48 runs | Low overlap and activity contrast, but the guided combined arm breaks parent links in all trials. |
| E019 | Three force budgets, three allocations plus reference, three seeds, two steps; 60 runs | Lower budgets retain junctions; selective allocation extends farther at matched force; no persistent narrow sprouts. |

## How these assays relate to self organization

E014 assigns fixed abstract identities to test mechanical compatibility. It does not show differentiation. E017 deliberately chooses a load site to test force transmission; that cell is not an emergent tip. E015 and E018 allow activity to evolve without assigning a leader, but use a prescribed environmental direction. E019 further imposes a global force budget, so it is a diagnostic intervention rather than a fully local biological model.

The later assays preserve the research goal by not assigning branch points or splitting events, but they have not achieved it. Their value is to identify which proposed ingredients can work in isolation and which fail when combined. The 264 runs across E014-E019 are parameterized computational trials, not 264 independent biological replicates.

The following sections describe the tested rules, quantitative outcomes, failures and caveats. Complete per-condition records remain in the experiment notebooks listed in Appendix A.

<!-- PAGE -->
# 9 Baseline activation growth and protrusion

## Initial baseline and delayed activation

At hour 10, every baseline seed has 256 cells and zero cells above SNAIL 0.8. At hour 12 [E002]:

| Seed | Cells | High SNAIL | Elongation | Mean radial rank of high group in E003 |
| --- | --- | --- | --- | --- |
| 42 | 919 | 17 | 1.084 | 0.092 |
| 43 | 849 | 73 | 1.139 | 0.421 |
| 44 | 938 | 24 | 1.125 | 0.199 |

Radial rank runs from central 0 to peripheral 1. Activation is therefore not equivalent to a peripheral leader. E003 has only five cells already high at the start of an available forward-motion interval; four divide during that interval. Their outward motion and neighbor retention do not establish persistent causal leadership.

## Delayed ablations

Turning protrusive force off during hours 10-12 lowers the mean outward-motion summary by 0.05868, 0.11206 and 0.11443 model length units per model hour for seeds 42, 43 and 44 [E004]. Positive net outward movement remains. The relative reductions of 33%, 30% and 44% depend on the interval aggregation and are not fractions of expansion causally attributable to protrusion.

The factorial growth/division arrest [E005] gives the following RMS spread changes from hour 10 to 12:

| Seed | Normal | Protrusion off | Growth and division off | Both off |
| --- | --- | --- | --- | --- |
| 42 | +7.08% | +5.57% | -0.06% | -1.60% |
| 43 | +28.72% | +23.16% | +1.46% | -4.75% |
| 44 | +13.59% | +8.63% | -1.23% | -5.11% |

Growth and division together support increased spread; these joint interventions do not separate their individual contributions. Protrusion also contributes to motion when growth is arrested. Marker counts can rise without protrusive force, reinforcing that marker activation and functional migration are different outcomes.

Conclusion: expansion and delayed regulatory activation occur, but neither demonstrates organized branching. Early connectivity dips are threshold-sensitive and can coincide with division, so they should not automatically be called permanent fragmentation.

<!-- PAGE -->
# 10 Longer observation and controlled timestep tests

## Growth arrest to hour 20

In E006, all three populations remain at 256 cells after arrest at hour 10. By hour 20, RMS spread changes are -10.47%, -44.19% and -44.30%; high-SNAIL counts are 0, 16 and 70. Seed 42 clusters, while the others become compact. No cell passes the full high-SNAIL/peripheral/connected two-hour persistence screen. Late seed-44 activity is right-censored, however, so the finite horizon cannot establish that it is short-lived.

An added dt=0.025 run changes both the earlier tissue history and future random draws. It is a sensitivity warning, not a controlled convergence experiment. E007 addresses those confounds with a complete seed-44 hour-10 checkpoint, verified by exact round-trip and future stochastic-update tests. All continuation arms arrest growth/division and set subsequent GRN draws to 0.5, canceling the noise term.

| E007 dt | Final RMS radius | High SNAIL | Mean SNAIL |
| --- | --- | --- | --- |
| 0.05 | 3.0524 | 34 | 0.5117 |
| 0.025 | 2.1869 | 71 | 0.6285 |
| 0.0125 | 1.5454 | 23 | 0.4357 |

The final position RMS difference between 0.025 and 0.0125 is 1.93948 model units. Identical initial state and absent future noise do not remove the late disagreement. The finest tested arm is not established ground truth.

## Isolating proposed contributors

E008 applies m(dt)=0.2^(dt/0.05), preserving scalar memory retention under subdivision. Because polarity is normalized, this does not prove equivalent directional dynamics. For the middle/fine pair, late mean position disagreement increases by 17.46% and GRN disagreement by 26.73%; the candidate rescaling is not adopted.

E009 freezes all six nodes and four outputs while mechanics remains active. Final RMS radii at the three dt values are 6.4869, 2.5493 and 2.6882; substantial sensitivity remains. E010 additionally removes protrusion. Large passive potential increases remain at dt=0.05 and 0.025, with 6 and 8 increases between saved observations. The finest arm has none at that sampling cadence. These results locate a mechanical numerical issue without proving that all remaining disagreement has one cause.

<!-- PAGE -->
# 11 Passive refinement and coupled mechanical substeps

## Passive refinement E011

The reduced passive helper retains original pair arithmetic/order, fixed radii and frozen cell-specific strengths, omitting state that cannot affect translation in this regime. The dt=0.0125 control exactly reproduces E010 positions and force proxies at all 21 observations. An optional compiled prototype failed the exact-replay criterion at roundoff scale and was not used in any scientific arm.

| Passive dt | Updates | Potential increases at every update | Final RMS radius |
| --- | --- | --- | --- |
| 0.0125 | 800 | 0 | 2.74950 |
| 0.00625 | 1,600 | 0 | 2.74031 |
| 0.003125 | 3,200 | 0 | 2.73437 |

The spread range is only 0.554% of the finest spread, but the middle/fine final position difference is 0.47558 model units, about 0.632 mean cell radii. Passing a dissipation check and agreeing in bulk spread do not establish identical cell arrangements.

## Coupled substeps E012 and E013

These diagnostics restore regulation and fields while holding the outer clock at dt=0.05. Polarity and propulsion eligibility update once per outer interval; forces update at each inner step. The first inner force norm is retained for the next GRN update, preserving the original lag. One substep exactly replays E007. All arms retain 256 cells and 200 outer updates.

| Substeps | Inner dt | Final RMS radius | High SNAIL | Mean YAP |
| --- | --- | --- | --- | --- |
| 1 | 0.05 | 3.0524 | 34 | 0.6701 |
| 4 | 0.0125 | 1.5054 | 4 | 0.6639 |
| 8 | 0.00625 | 1.3777 | 2 | 0.2957 |
| 16 | 0.003125 | 1.4136 | 6 | 0.2907 |

The 8/16 final spread differs by 2.61% relative to the eight-step result. Final GRN RMSE falls from 0.18859 for 4/8 to 0.08028 for 8/16, but final position disagreement changes from 0.47745 to 0.49745. Some regulatory and aggregate measures improve; cell-level convergence is not established. No universal timestep is adopted, and the outer splitting policy itself remains an assumption.

<!-- PAGE -->
# 12 Pair specific adhesion and sorting

## Hypothesis and controlled assay E014

The original scalar adhesion product cannot independently express strong same-type binding and weak cross-type binding. E014 tests a symmetric pair multiplier in a passive, jittered 4 x 4 x 4 aggregate of 64 cells. Two fixed abstract identities are assigned randomly and equally. There is no GRN evolution, growth or propulsion. Three seeds, four conditions and dt=0.003125/0.0015625 run for five model units.

Equal adhesion uses multiplier 1 for all pairs. Selective adhesion uses 1 for like identities and 0.2 for unlike identities. A uniformly reduced control matches the selective arm's initial summed attraction magnitudes over interacting pairs. A compensated selective arm instead rescales both coefficients to match the equal arm's initial summed magnitude. This matching does not continue automatically as geometry changes and does not match vector force sums.

| Condition | Unlike contact fraction | RMS radius | Elongation |
| --- | --- | --- | --- |
| Equal adhesion | 0.512 | 0.656 | 1.003 |
| Uniformly reduced | 0.513 | 0.779 | 1.008 |
| Selective | 0.250 | 0.854 | 1.636 |
| Selective compensated | 0.320 | 0.724 | 1.721 |

Values are final means across seeds at the finer step. All final collectives remain one component at geometric cutoff 1.2. The selective unlike-contact fraction ranges 0.246-0.253, versus 0.510-0.518 for the matched uniform control. Pair compatibility causes sorting and shape change beyond simply lowering general stickiness in this model.

![Adhesion comparison](figures/adhesion.png)

Figure 1. Finer-step seed means and min/max bands, not confidence intervals. Sorting remains connected. Soft repulsion permits strong compaction, so geometry is not a calibrated tissue prediction. The identities are not measured protein expression or actual progenitor/TYPE_A differentiation. Sorting and elongation do not establish branching.

<!-- PAGE -->
# 13 Guidance and competition in a synthetic sheet

## New assay rules E015

The later sprouting assays use a planar 6 x 6 sheet embedded in three-dimensional coordinates: 36 cells of radius 0.5, with 12 cells in the two left columns fixed as parent tissue and 24 mobile cells. Initial spacing is 1 with positional jitter 0.025. The prescribed VEGF-like cue is V(x)=sigmoid((x+2)/3), increasing in +x without lateral branch paths. It is not generated by simulated hypoxic tissue.

Abstract ligand/tip activity D relaxes toward V/(1+(N/0.3)^4) when competition is enabled, otherwise toward V. N relaxes toward mean D of geometric neighbors within 1.2. Both timescales are 0.5. Guided traction points +x; unguided traction follows the original outward geometric direction. Magnitude is 0.15+2.35D. Initial D is perturbed around 0.5, N=0, with no future noise. There is no growth, lumen, ECM or flow.

The four conditions are neither mechanism, guidance alone, competition alone and both. Cue magnitude affects D in every arm; the guidance factor changes force direction only. E015 initially retains the original soft pair mechanics. Three seeds and 16/32 substeps per outer dt=0.05 give 24 trials to time 12.

| Condition | Outcome repeated across seeds and steps |
| --- | --- |
| Neither | Contracts into the parent region |
| Guidance | Mobile cells detach and migrate as a cluster |
| Competition | Activity is reduced, but tissue contracts |
| Both | Migration is slower; mobile cluster still detaches |

No persistent connected narrow extension is detected. Maximum overlap reaches 0.876 for a nominal diameter of 1, a severe mechanical limitation. The circuit does not create strong mobile leader/follower contrast in this moving, contracting assay. That result motivates testing signaling on a fixed network and testing mechanics under imposed loads separately.

This is a negative result for this implementation and parameter set. It does not refute VEGF-guided angiogenesis or Dll4-Notch biology. It also does not test explicit surface proteins, cell deformation, matrix resistance or replacement of trailing cells by proliferation.

<!-- PAGE -->
# 14 Static signaling reveals a geometry dependence

## Separating pattern formation from movement E016

Hold positions fixed and use either the initial sheet network or the contracted E015 control network. Give all cells a uniform cue V=0.7 to remove positional cue variation. Compare no inhibition, the reference Hill exponent 4/half inhibition 0.3, and a steeper exponent 8/half inhibition 0.2. Three seeds, two dt values and two networks produce 36 runs to time 20.

An isolated high cell has D>0.6 with no high contact neighbor. A stable high cell is high in at least 80% of observations during times 15-20. These are signaling proxies; no cell moves or leads tissue in this test.

| Network | Rule | Stable isolated high cells per seed | High cells in formerly mobile group |
| --- | --- | --- | --- |
| Initial sheet | No inhibition | 0 | 24 |
| Initial sheet | Reference | 15 | Mean 9.33 |
| Initial sheet | Steeper | 15 | Mean 9.33 |
| Contracted | No inhibition | 0 | 24 |
| Contracted | Reference | 2 | 0 |
| Contracted | Steeper | 4 | 0 |

Without inhibition all cells approach 0.7; none is an isolated high cell. The reference circuit can therefore produce contrast on a sparse fixed sheet, revising the interpretation of E015. On the contracted graph, high cells occur only at former fixed-parent sites. Strengthening inhibition does not restore high activity among the formerly mobile cells in these tests.

![Static signaling states](figures/static_signaling.png)

Figure 2. Reference circuit, seed 42, dt=0.025, final abstract activity on two networks. The initial-sheet pattern is dense: 15 high cells among 36 is not evidence of one sparse sprout leader. The static horizon, uniform cue and fixed graph differ from E015; this is not a causal replay of its failure. Final step-pair activity RMSE is at most 0.000979, supporting numerical consistency for this assay rather than biological validation.

<!-- PAGE -->
# 15 Mechanical integrity under imposed traction

## Mechanical variants E017

The same initial sheet is tested without signaling under three loads: none, a +x force of 2.5 on one front cell, or 2.5 on every mobile cell. The chosen front cell is closest to y=0 in the front column. This is an imposed mechanical probe, not self-organized leadership. Uniform loading has 24 times the total force of the single-cell probe.

Four variants compare original repulsion 4/cutoff 2.5, repulsion 40/cutoff 2.5, repulsion 40/short cutoff 1.2, and repulsion 40 with explicit material junctions. The first three retain attraction 0.49. Material junctions instead connect initial geometric contacts with spring constant 4 and rest lengths equal to initial separation; edges break irreversibly beyond 1.6. No new edges form. All overlapping pairs repel, and damping remains 2. There are 72 runs over three seeds, two dt values and five model time units.

| Variant under one cell load | Final attached cells | Mean loaded cell displacement | Worst overlap across seeds |
| --- | --- | --- | --- |
| Original | 24/24 | -3.315 | 0.873 |
| Stiffer repulsion | 24/24 | +0.400 | 0.251 |
| Short range | 23/24 | +6.085 | 0.090 |
| Material junctions | 24/24 | +1.165 | 0.090 |

The material arm passes the declared screen: all cells remain attached at saved geometric observations, no junction breaks under the single-cell load, and overlap stays below 0.1 at both steps and all seeds. Mean mobile displacement is 0.194, so the loaded-cell movement is not simply whole-sheet translation. The worst overlap is already present in the initial jittered sheet.

No-load controls show that short-range attraction and material junctions preserve the sheet better than broad attraction. But short-range attraction lets the pulled cell escape. At the much larger all-mobile load, every variant detaches the mobile tissue; the material arm breaks six parent junctions.

Conclusion: reduced overlap alone is insufficient. Explicit junctions can transmit one localized load under these conditions, but are not calibrated VE-cadherin mechanics. The fixed material graph may restrict cell rearrangement, and an imposed successful pull is not a branch. These tests justify a coupled trial, not production adoption.

<!-- PAGE -->
# 16 Recombining signaling guidance and junctions

## Coupled test and timing control E018

Combine the reference D/N circuit, guidance and E017 material mechanics without assigning a leader. Compare neither/guidance/competition/both. Start movement either immediately or after 20 stationary units of the same condition-specific signaling. All geometries remain matched; prepatterning deliberately changes initial regulatory state. Three seeds and 16/32 mechanical substeps give 48 runs to time 12, with outer dt=0.05.

The original traction mapping gives every mobile cell a positive force, 0.15+2.35D. In the combined arm, prepatterning produces 1-3 initially high mobile cells instead of zero, but does not prevent later separation. All six links from parent to mobile tissue break in every seed and both resolutions. At the finer step the first sampled fully disconnected state is time 5.0 for both timings.

The combined arm now shows activity contrast, with final within-mobile SD about 0.383 immediate and 0.377 prepatterned. However, that contrast is observed in detached tissue and includes approximately 10-12 high cells among 24 mobile cells, not a single narrow advancing tip. Guidance alone drives all mobile cells high. Maximum overlap remains below 0.09, so gross collapse is improved without solving attachment.

![Junction sprouting states](figures/combined.png)

Figure 3. Combined guidance/inhibition arm, seed 42, finer step, comparing the two timings. Blue lines show intact final junctions; red crosses mark fixed parent cells. Panels have independent axes. Both lose parent attachment. Separately, unguided controls can retain material paths while losing geometric contact.

No persistent geometric narrow extension is detected. The material-path variant also detects no candidate in the combined arm. A transient material-only detection in one prepatterned neither control does not establish a stable sprout. Geometric signaling ends at distance 1.2 while junctions survive to 1.6, an explicit model mismatch. Final position disagreement reaches 0.456 across resolutions, despite repeated qualitative combined-arm failure. Stronger contrast and low overlap are therefore necessary progress, but insufficient results.

<!-- PAGE -->
# 17 Matched force budgets distinguish scale and allocation

## Fixed mechanics with different traction E019

Keep the reference inhibitory circuit, prepatterning, +x guidance and material mechanics fixed. Test total force F=2.5, 7.5 or 15 with allocation weights w_i=b+(1-b)D_i^4 and force F w_i/sum(w). b=0 favors high activity, b=0.25 adds a follower floor, and b=1 is uniform. b is not the fraction of total force delivered to followers. A reference retains the E018 force rule. Three seeds and two resolutions give 60 runs.

Global normalization is a diagnostic intervention using all mobile cells, not a local cell-autonomous biological controller. Detached cells remain included in that allocation. No leader ID is selected, but activity weighting itself is an uncalibrated model assumption.

| F | Activity weighted extent | Follower floor extent | Uniform extent |
| --- | --- | --- | --- |
| 2.5 | 0.532 | 0.300 | 0.283 |
| 7.5 | 1.489 | 0.874 | 0.772 |
| 15 | 2.945 | 1.451 | 0 |

Extents are mean final material-connected max x at the finer step, not branch lengths. At F=2.5 and 7.5, all 36 seed/allocation/resolution trials retain every mobile cell through material paths, with zero ruptures. At F=7.5 activity weighting yields 1.93 times the mean extent of uniform traction at the same force.

At F=15, b=0 retains 24/24, 23/24 and 23/24 mobile cells across seeds; b=0.25 retains 0/24, 24/24 and 24/24; uniform loading retains none. Every high-budget allocation causes junction damage. The reference mean total force is 23.13-24.06, and it detaches in all runs.

![Force budget comparison](figures/force_budget.png)

Figure 4. Finer-step seed means with min/max error bars, not confidence intervals. Low load protects junctions; allocation matters at matched load. No persistent narrow extension is found by either screen. A material-only candidate in one high-selective seed at time 3.5 does not persist. Some intact spring configurations are geometrically disconnected. Numerical checks improve, with maximum resolution position RMS 0.00903 and exact six-run E018 reference replay, but the shape objective remains unmet.

<!-- PAGE -->
# 18 Supported findings and unresolved explanations

## Findings supported within the tested protocols

The original model exhibits delayed SNAIL activation, growth-related expansion and a contribution of protrusion to outward motion. It does not establish sustained functional tip behavior. Large-step passive mechanics can violate its expected dissipation property; finer stepping improves some diagnostics without proving full trajectory convergence.

Synthetic pair-specific adhesion changes sorting beyond a uniform reduction in attraction. A reference inhibitory circuit creates static high/low patterns on an uncontracted graph. Explicit material junctions with stronger repulsion transmit a single imposed load while avoiding gross overlap. Under coupled movement, activity contrast and low overlap can coexist with parent rupture. At matched force budgets, allocation changes extension and failure mode; lower loads consistently retain junctions over the tested finite horizon.

These are positive mechanistic and numerical findings, even though they do not achieve branching. They narrow the question from whether movement or activation occurs to how a connected population can continually reorganize into a narrow extension.

## Claims not supported by this record

- That the original rules can never branch, or that a single missing ingredient has been identified conclusively.
- That high SNAIL, high D, elongation, graph degree or large max x establishes a functional vessel branch.
- That TYPE_A has a measured protein repertoire or that the tested adhesion coefficients represent actual cadherin affinities.
- That the later assays are validated angiogenesis, or that E019's global budget is a biological feedback mechanism.
- That any timestep, force allocation or junction coefficient is universally correct or ready for adoption.

## Structural questions that the latest assays leave open

Most late assays remove growth and cell supply, use cell centers rather than deformable shapes, and retain initial springs without junction formation or rest-length remodeling. Such a network can deform a sheet while resisting the neighbor rearrangements that might support a narrow extension. This is a plausible limitation, not an experimentally isolated cause of failure.

Other unresolved issues include geometric versus material signaling contacts, force-proxy cancellation, load sharing, substrate/ECM coupling and how leader density should be regulated. The current data justify choosing a specific structural hypothesis with supervisor guidance; they do not justify silently adding all missing biology or programming branching events directly.

<!-- PAGE -->
# 19 Questions and decisions for the supervisor

## Proposed opening statement for the meeting

The project aims to produce vessel-like branching from cell-level self-organization. So far, I have separated growth, active motion, signaling, adhesion and junction mechanics, and tested numerical sensitivity. I can produce sorting, signaling contrast and connected extension under some conditions, but not a robust persistent sprout. I would like guidance on which structural hypothesis to test next and what biological evidence should constrain it.

## Decisions requested

1. What is the appropriate first biological target: an early connected sprout, a bifurcating cellular network, or a lumenized vessel? What model simplifications are acceptable for that target?
2. Is a center-based cell representation sufficient, or should deformable cells, contact areas or an alternative tissue representation be considered before further parameter sweeps?
3. Should the next isolated mechanism be local junction turnover and neighbor exchange, growth and cell supply, or ECM-dependent traction? Which can be anchored to accessible experimental data?
4. Should lateral inhibition act through actual retained junctions, geometric proximity or contact area? What is a defensible model of tip/stalk plasticity and leader density?
5. What quantitative evidence should define success: branch width, parent continuity, persistence, skeleton bifurcation, division/rearrangement histories and sensitivity across seeds?
6. What numerical tolerances, number of independent seeds and biological calibration are needed before reporting a mechanism as a research result?

## Candidate next studies for discussion rather than adoption

| Candidate | Isolated question | Evidence needed before combining |
| --- | --- | --- |
| Local junction remodeling | Can cells exchange neighbors while maintaining tissue continuity under bounded load? | Logged bond/contact changes, limited overlap, persistent parent attachment and reversible/irreversible rules justified independently |
| Growth and cell supply | Does local proliferation support extension rather than broad expansion or buckling? | Matched no-growth controls, division histories, density and narrow-extension measures |
| Substrate or ECM coupling | Can local traction and resistance stabilize directed migration? | Explicit substrate rule, controlled resistance contrast and calibrated or justified force scales |

A useful scope would choose one of these questions, retain a stable mechanical baseline and specify the success/failure criteria before execution. These studies have not been performed. Branch positions and splitting events should remain emergent; using an explicitly branching reference algorithm could be a separate engineering comparison, but would not answer the central cell-level research question.

<!-- PAGE -->
# 20 Reproducibility validation and study limitations

## What has been preserved

The research notebook index and experiment folders retain protocols, assumptions, observed failures, run settings, scripts, CSV measurements, NPZ observations and figures. SHA-256 manifests identify archived artifacts. E007 additionally preserves a complete compressed JSON checkpoint including cell/GRN/mechanics/field attributes, array shapes and dtypes, dictionary ordering and PCG64 RNG state. Its canonical payload digest is:

2c846bb447e3193e79998844de7b19ee15a65557cf132a1369c327f367f30c95

Selected NPZ observations should not be confused with complete restart checkpoints. Early observations omit some fields and generator state. The report is a synthesis of archived computational evidence, not a new simulation run. Source snapshots and an evidence manifest accompanying this report identify the files inspected for its creation.

## Validation examples and their scope

Checkpoint tests verify exact restoration and future updates; intervention tests check timing and restoration; force/potential tests compare negative gradients and implementation behavior; fixed-clock tests check equal biological updates under substepping. E011 and E012/E019 historical replays preserve specified observables exactly. Different test counts in successive notebooks are cumulative or overlapping suites, not independent replications to be summed.

In E019, 13 relevant tests passed, force-budget error stayed below 7.11e-15, no zero-weight fallback occurred, all 1,500 saved states were finite, anchors were exactly fixed and no work-adjusted potential increase exceeded tolerance across 345,600 mechanical updates. Maximum overlap was 0.08981, already present initially. These checks substantiate the implementation and reported comparisons, not the biological adequacy of the assumptions.

## Practical reproduction

Run commands from the repository root with the documented Python/NumPy environment; shell commands follow the repository's rtk convention. Consult the specific README before selecting a run: the later scripts are separate models with different states and clocks. Use fresh output folders to avoid overwriting evidence. A fixed-dt headless run should accompany any interactive demonstration.

The records are local working-tree artifacts; no external backup, publication or peer review is implied. Some controls are referenced from earlier folders rather than duplicated. A source hash establishes identity, not correctness. There is no calibrated mapping to physical tissue, no biological dataset validation, no inferential statistical analysis and no exhaustive parameter search. These limits should accompany any presentation of the results.

<!-- PAGE -->
# Appendix A Experiment records and evidence locations

All paths below are relative to the project root C:/Users/ayush/omniverse/kit-app-template. Under research/experiments, each listed folder contains its README and available settings, CSVs, figures, archived scripts and provenance files. The initial baseline is documented separately.

| Record | Notebook location |
| --- | --- |
| Initial baseline | tools/BRANCHING_BASELINE.md |
| E002 | 2026-11-06-long-baseline |
| E003 | 2026-11-06-trajectory-analysis |
| E004 | 2026-11-06-protrusion-ablation |
| E005 | 2026-11-06-growth-arrest |
| E006 | 2026-11-06-fixed-population-long |
| E007 | 2026-11-06-shared-state-timesteps |
| E008 | 2026-11-06-polarity-memory-scaling |
| E009 | 2026-11-06-frozen-regulation |
| E010 | 2026-11-06-passive-mechanics |
| E011 | 2026-11-06-passive-refinement |
| E012 | 2026-11-06-coupled-substeps |
| E013 | 2026-11-07-coupled-refinement |
| E014 | 2026-11-07-pair-specific-adhesion |
| E015 | 2026-11-07-guidance-tip-competition |
| E016 | 2026-11-07-static-tip-selection |
| E017 | 2026-11-07-loaded-tissue-mechanics |
| E018 | 2026-11-08-junction-sprouting |
| E019 | 2026-11-08-traction-budget |

The notebook index is research/README.md. Production configuration is source/apps/embryo_twin.digital.kit. Python application rules are under source/extensions/embryo_twin.sim/embryo_twin/sim. Research runners are tools/branching_*.py, with corresponding test scripts. Later synthetic assays are not available by pressing Start in the Kit UI unless integrated in a future explicitly approved implementation.

Appendix locations are intended to make the report auditable in the shared project; a professor receiving only the PDF will not automatically receive the raw NPZ data. The accompanying evidence package includes notebooks, application source snapshots and selected summary data, while full raw data remain in the repository.

<!-- PAGE -->
# Appendix B Biological references and notation

## Biological motivation

[B1] Gerhardt H et al. VEGF guides angiogenic sprouting utilizing endothelial tip cell filopodia. Journal of Cell Biology 161, 1163-1177 (2003). DOI: 10.1083/jcb.200302047. https://pmc.ncbi.nlm.nih.gov/articles/PMC2172999/

[B2] Hellstrom M et al. Dll4 signalling through Notch1 regulates formation of tip cells during angiogenesis. Nature 445, 776-780 (2007). DOI: 10.1038/nature05571. https://doi.org/10.1038/nature05571

[B3] Jakobsson L et al. Endothelial cells dynamically compete for the tip cell position during angiogenic sprouting. Nature Cell Biology 12, 943-953 (2010). DOI: 10.1038/ncb2103. https://doi.org/10.1038/ncb2103

[B4] Bentley K et al. The role of differential VE-cadherin dynamics in cell rearrangement during angiogenesis. Nature Cell Biology 16, 309-321 (2014). DOI: 10.1038/ncb2926. https://doi.org/10.1038/ncb2926

These studies motivate mechanism choices; none supplies a calibration for the present coefficients. They are not an exhaustive literature review. The computational claims in this report are supported by the experiment records in Appendix A and the inspected code, not by inferring successful simulation behavior from biological papers.

## Notation and distinctions

| Symbol or term | Meaning in this report |
| --- | --- |
| dt or h | Outer timestep or mechanical inner step, as specified for each assay |
| D and N | Abstract ligand/tip activity and inhibitory state in synthetic assays |
| SNAIL and Notch | Original model GRN nodes; not numerically interchangeable with D and N |
| A and B | Original per-cell adhesion and protrusion phenotype strengths |
| F and b | E019 total force budget and follower-floor weight parameter |
| U and W | Model pair potential and held-active-force-adjusted potential |
| Checkpoint | Complete restart state when explicitly described as complete |
| Observation | Selected saved state used for analysis, not necessarily sufficient for restart |
| Stable pattern | Persistent regulatory contrast on a network; not a vessel |
| Connected extension | Advance with a specified geometric or material connection; not automatically a branch |

The central unresolved question is how local cell-level behavior can generate a persistent, narrow, growing tissue extension while preserving continuity. Supervisor input is requested before choosing the next structural hypothesis and success criteria.
