# Branching baseline

Research log: [E002: longer unchanged-model experiment](../research/experiments/2026-11-06-long-baseline/README.md).
The [research notebook index](../research/README.md) links subsequent analyses
and the delayed protrusion ablation. The runner's optional
`--protrusion-off-at` switch is an explicit intervention; omitting it keeps
the normal force rules. The onset must align with dt and precede the horizon.
An optional `--growth-off-at` switch arrests radius growth and division while
allowing GRN, age and field updates to continue. It can be combined with
protrusion ablation; E005 documents this explicit experimental protocol.
That folder preserves raw observations, source copies, settings and figures.
The runner now also records source hashes, regulatory diagnostics, parent
division counts and cell centers outside the field. It streams observations
to disk, refuses nonempty output folders, and stops at the application cell
limit rather than silently advancing beyond it.

Run from the repository root with a Python environment containing NumPy:

```powershell
rtk proxy python tools/branching_baseline.py --hours 10 --dt 0.05 --seeds 42 43 44
```

This loads the pure simulation without starting Kit or changing its rules.
Each seed initializes the existing two-cell setup with its own random generator.
Output goes into the ignored `_branching_baseline/` folder. Use a distinct
`--output` folder for each experiment to avoid overwriting earlier results.

Outputs include run settings, per-seed CSV measurements, and numbered NPZ
snapshots containing time, cell IDs, positions, radii and SNAIL values.
Times are model hours, not wall-clock hours. Results describe this headless
timestep and need not match an interactive run using a different timestep.

Connectivity uses center distance <= 1.25 times the sum of cell radii.
Change `--contact-factor` to assess sensitivity. This is a proximity proxy,
not a measurement of actual junctions. Gaps immediately after division may
be temporary. Elongation is the square root of the largest to smallest
spatial covariance eigenvalue ratio, regularized by cell size. A sphere-like
cloud is near 1; a large value alone does not establish branching.
High-SNAIL means SNAIL > 0.8, matching the existing model threshold.
No branch count is reported: contact graph degree would mislabel compact
tissue as branching. Persistence and branch geometry need separate analysis.

## Initial observations

With the existing rules, dt=0.05 and 10 model hours:

| Seed | Cells | Components | Largest component fraction | Elongation | High-SNAIL cells |
| --- | --- | --- | --- | --- | --- |
| 42 | 256 | 52 | 0.3555 | 1.0923 | 0 |
| 43 | 256 | 2 | 0.9961 | 1.1595 | 0 |
| 44 | 256 | 2 | 0.9961 | 1.1437 | 0 |

These short runs do not demonstrate sustained branching. They indicate a
roughly isotropic distribution, no above-threshold SNAIL cells at the final
sample, and seed-sensitive proximity connectivity. They do not prove that
longer runs cannot branch. Inspect temporal snapshots, contact-threshold
sensitivity and timestep sensitivity before selecting changes to the model.

## Follow-up sensitivity checks

Reanalyzing the same final snapshots with different proximity tolerances:

| Seed | Factor | Components | Largest component fraction |
| --- | --- | --- | --- |
| 42 | 1.00 | 108 | 0.2344 |
| 42 | 1.25 | 52 | 0.3555 |
| 42 | 1.50 | 8 | 0.9609 |
| 43 | 1.00 | 6 | 0.9375 |
| 43 | 1.25 | 2 | 0.9961 |
| 43 | 1.50 | 1 | 1.0000 |
| 44 | 1.00 | 26 | 0.7383 |
| 44 | 1.25 | 2 | 0.9961 |
| 44 | 1.50 | 1 | 1.0000 |

Connectivity is strongly threshold-sensitive. Do not call these component
counts evidence of permanent fragmentation.

Repeating seed 42 at dt=0.025 through 10 hours yielded 256 cells, one
component at factor 1.25, elongation 1.1276, and zero high-SNAIL cells.
Compare with 52 components at dt=0.05. This is evidence that connectivity
is not yet robust to experimental settings, not a numerical convergence
study: changing dt also changes the number of random draws and stochastic
forcing. Multiple seeds and finer steps are needed to separate these effects.

The saved hourly baseline snapshots have maximum SNAIL below 0.002 through
hour 8, then approximately 0.14-0.15 at hour 9 and 0.42-0.44 at hour 10.
Thus SNAIL activation is beginning rather than demonstrably failing.

The existing `grn_debug.py --pair --no-plot` experiment at dt=0.1, seed 1,
pressure 3 and morphogen 0.5 separates the two cells by about hour 15:
SNAIL is approximately 0 and 1.67. It reaches 0 and 2 by hour 20.
This demonstrates symmetry breaking under sustained, prescribed inputs;
it does not establish that the coupled tissue supplies those inputs long
enough or maintains connected branches.

Next hypothesis: the coupled tissue needs sufficient sustained regulatory
activation time for spatially separated motile states to emerge. Test a
longer observation window before changing regulatory or mechanical rules.
Track growth and division events alongside connectivity, because synchronized
division may produce transient proximity gaps. The current 1500-cell limit
and quadratic mechanics cost must be accounted for in longer experiments.
