# E016: static-network tip selection

Date: 2026-11-07 (Europe/Berlin). Protocol recorded before execution.

Hold all cells fixed and test whether the abstract inhibitory circuit can produce leader/follower contrast. Use the E015 initial sheet and its final contracted control topology for seeds 42/43/44; graph edges are geometric contacts <=1.2. Use uniform cue V=0.7 in all cells to remove positional cue variation. Initial D is the identical seeded E015 perturbation around 0.5, N=0. Compare no inhibition, E015 Hill exponent 4/half inhibition 0.3, and steeper exponent 8/half inhibition 0.2. Both states relax with tau=0.5 for 20 model units at dt 0.05 and 0.025. No future noise.

Primary screen: isolated high-activity cells D>0.6 with no high-activity contact neighbor, alongside low cells D<0.2 and edge contrast. Persistent high cells are high for >=80% of saved observations over times 15-20; persistent isolated cells apply the same no-high-neighbor criterion to that persistent set. These are operational signaling proxies, not physical tip cells. Report all arms and both networks: do not infer that a circuit works in dense tissue from a sparse grid alone.

No mechanics, migration or branching is tested here. Varying both Hill parameters is an exploratory circuit contrast, not an isolated effect of either parameter. The contracted topology is a deliberately confounded geometric/topological comparison, not solely a degree manipulation. Application code unchanged.

## Results and interpretation

Completed 36 runs, saving 1,476 signaling states. At dt=0.025 and time 20:

| Fixed network | Inhibition | Stable isolated high cells, each seed | High cells in previously mobile group, seed mean | Mean edge contrast |
| --- | --- | --- | --- | --- |
| Initial sheet | None | 0 | 24 | 0 |
| Initial sheet | E015 reference | 15 | 9.33 | 0.554 |
| Initial sheet | Steeper | 15 | 9.33 | 0.606 |
| Contracted E015 | None | 0 | 24 | 0 |
| Contracted E015 | E015 reference | 2 | 0 | 0.029 |
| Contracted E015 | Steeper | 4 | 0 | 0.171 |

Without inhibition, all cells approach the uniform cue value 0.7; all are
high and none is an isolated high cell. The reference circuit DOES create
stable high/low contrast on the sparse sheet. Thus E015 did not demonstrate
that the circuit is intrinsically unable to pattern. The contracted graph
changes the outcome: high cells occur only at the former fixed-parent
locations, with none among the formerly mobile population. Increasing the
Hill steepness/altering the half-inhibition parameter does not fix that
mobile-group outcome in these tests.

The initial-sheet reference state has 15 high cells among 36, roughly 42%.
This is a spatial pattern, not a demonstration of a sparse number of vessel
leaders, and it does not select cells at a moving sprout front. No cells
move here; the original parent/mobile labels are only readout partitions.
Topology is held fixed for each run, all cells receive the same cue, and
the horizon is longer than E015. Consequently this is not a causal replay
of E015's failure; evolving contacts, cue exposure and timing remain coupled
in that earlier experiment. Geometry/contact preservation is a plausible
priority for recombination, not a uniquely proven cause of failed sprouts.

Final activity RMSE between dt=0.05/0.025 is at most 0.000979, and matched
half-unit trajectory RMSE at most 0.005742. The high/low-pattern conclusions
repeat across all tested seeds and both steps. These are numerical
consistency observations, not biological validation or an exhaustive
parameter sweep.

## Validation and artifacts

Five prerequisite tests passed. An initial bitwise comparison of reference
forces failed at floating roundoff (max absolute difference 1.8e-15): this
diagnostic uses decimal 0.49, while the prior helper computes 0.7 squared.
The force-equivalence test was corrected to explicit tight floating-point
tolerances (rtol 1e-13, atol 2e-14), rather than claiming bitwise replay.
No force implementation was changed after execution began.

All arrays are finite, initial D values match their seeds, graphs are
symmetric, and original application source hashes match. Archived scripts,
settings, input hashes for E015 topologies, NPZ timecourses, summary CSVs,
resolution differences, validation JSON and `sha256.json` preserve provenance.
`signaling_states.png` shows seed 42 at the finer step; positions remain
fixed, and colors indicate abstract activity rather than measured proteins.

```powershell
rtk proxy python tools/test_branching_prerequisites.py
rtk proxy python tools/branching_prerequisites.py signaling --output research/experiments/2026-11-07-static-tip-selection
rtk proxy python tools/branching_prerequisites_report.py signaling research/experiments/2026-11-07-static-tip-selection
```

Use a fresh output folder to rerun. No production-model choice is adopted.
