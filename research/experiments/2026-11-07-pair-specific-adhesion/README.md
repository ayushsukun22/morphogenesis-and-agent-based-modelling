# E014: pair-specific adhesion assay

Date: 2026-11-07 (Europe/Berlin). Protocol recorded before execution.

Question: does lower heterotypic adhesion alter mixing and collective structure independently of simply reducing overall attraction? This is a controlled passive assay, not a full coupled branching simulation.

Use 64 equal-radius cells in a jittered 4x4x4 aggregate, with balanced randomly assigned fixed abstract A/progenitor identities. Seeds 42/43/44 vary geometry and identity placement. These identities do not assert actual fate or measured protein expression. Freeze regulatory/chemical dynamics, growth and propulsion. Keep original cutoff 2.5, repulsion 4, per-cell adhesion 0.7, damping 2. Add a symmetric pair multiplier to adhesion only.

Four matched conditions: equal (all multipliers 1); selective (same-type 1, cross-type 0.2); uniform_reduced (all equal to the selective mean multiplier over initial interacting pairs); selective_compensated (selective multipliers divided by that mean). Thus selective versus uniform_reduced have equal initial summed attraction magnitudes, as do compensated versus equal. This matching applies initially only and does not match vector force sums or energy. Cross-type reduction of 80% is an exploratory parameter, not a measured biological value.

Run 5 model time units at dt 0.003125 and 0.0015625. Save half-unit observations and final positions. Primary readout: fraction of geometric contact edges (center distance <=1.2) joining unlike identities; secondary readouts: contact count, RMS radius, connected components, largest component, elongation and type-centroid separation. Contact threshold is fixed independently of the 2.5 force cutoff. Lower unlike-contact fraction may reflect sorting or fragmentation: interpret together with connectivity. No branch claim from elongation alone.

Validate equal-arm force equivalence to original passive helper, pair selectivity, force/potential gradient and reproducible balanced initialization. Check finite states, centroid drift, every-step passive potential increases and timestep disagreement. Application source hashes must remain unchanged. This assay cannot establish coupled GRN effects, genuine protein compatibility, robust branching, or behavior in growing tissue.

## Results

All 24 runs completed (3 seeds x 4 conditions x 2 timesteps). Four meaningful
implementation tests passed. Across 57,600 integration updates, no passive
potential increases above the stated relative tolerance were detected.
States remained finite, centroid drift remained below 1e-14 model units,
and application-source hashes were unchanged.

Final values at the finer dt, means across the three seeds:

| Condition | Unlike-type contact fraction | RMS radius | Elongation |
| --- | --- | --- | --- |
| Equal adhesion | 0.512 | 0.656 | 1.003 |
| Uniformly reduced adhesion | 0.513 | 0.779 | 1.008 |
| Selective adhesion | 0.250 | 0.854 | 1.636 |
| Selective, initial strength compensated | 0.320 | 0.724 | 1.721 |

All final states are single connected components at the specified contact
threshold. Selective unlike-contact fractions span 0.246-0.253 across seeds,
versus 0.510-0.518 for the strength-matched uniform control. Type-centroid
separation is 1.165-1.168 with selective adhesion, versus 0.124-0.188 for
the uniform-reduced arm. Identity-dependent attraction therefore produces
sorting and anisotropy in this assay; merely lowering attraction does not
reproduce it. The initial-strength-compensated arm also sorts, supporting
an effect beyond reduced initial total attraction magnitude.

Refining dt changes final matched-cell positions by RMS 0.00154-0.04237
model units across arms/seeds (maximum about 0.085 cell radii). These are
not identical trajectories. The sorting result persists at both tested
steps; this does not establish universal timestep convergence.

## Interpretation and limits

This supports the user's hypothesis that pair compatibility can matter,
within an explicitly abstract mechanical model. It does not identify real
adhesion proteins, validate the coefficient ratio, or demonstrate that the
existing GRN should assign these identities. The labels are random, fixed,
balanced groups, not current fate classifications. The tiny synthetic
aggregate and soft repulsion permit strong compaction/overlap, so lengths
and shapes must not be interpreted as calibrated tissue predictions.

No branching metric was measured and no branching is established. Reduced
mixing and elongation are evidence for sorting/shape change, not vessel
formation. A logical next experiment is to couple pair compatibility to
the existing continuous regulatory states and restore propulsion, keeping
matched controls; that has NOT been performed here. Protein expression
and pair compatibility would remain separate assumptions to document.

## Reproduction and files

```powershell
rtk proxy python tools/test_branching_adhesion_assay.py
rtk proxy python tools/branching_adhesion_assay.py
```

The runner uses this fixed experiment output folder and refuses to overwrite
completed results. Preserve the archive and select a fresh folder in a copy
before rerunning. `settings.json` records source and runner hashes;
`timecourses.csv`, `final_summary.csv`, `resolution_differences.csv` and
initial/final NPZ arrays preserve observations. `adhesion_comparison.png`
shows finer-step seed means and min/max bands (not confidence intervals).
Archived scripts and `sha256.json` preserve implementation and provenance.
