# E015: guidance and tip competition factorial assay

Date: 2026-11-07 (Europe/Berlin). Protocol recorded before execution.

Question: do directional guidance, local tip competition, or both generate persistent connected narrow extensions? Four matched conditions: neither, guidance only, competition only, both. Three initialization seeds (42/43/44); outer regulatory step 0.05, mechanical 16/32 substeps (0.003125/0.0015625), duration 12 model time units.

This is a separate abstract endothelial assay, not a patch of the original GRN or a calibrated angiogenesis model. Start with a planar 6x6 cell aggregate, center spacing 1, jitter 0.025, radius 0.5. Two left columns (12 cells) are fixed as a parent tissue reservoir; 24 cells move. All cell pairs retain equal adhesion (0.7 squared), original cutoff 2.5, repulsion 4, damping 2. No growth, division, ECM, lumen, perfusion, or evolving oxygen/chemical field is modeled.

The prescribed VEGF-like cue is V(x)=sigmoid((x+2)/3), with uniform +x gradient and no lateral paths or preselected tip cells. Cue magnitude influences tip activity in ALL arms. The guidance factor changes movement DIRECTION only; unguided direction is opposite the sum of unit directions toward neighbors within 2.5. Direction and active force are held within each outer interval.

Abstract tip/Dll4 activity D relaxes toward V/(1+(N/0.3)^4) when competition is enabled, and toward V otherwise. N relaxes toward mean contact-neighbor D (center distance <=1.2), including fixed-parent cells. Both timescales are 0.5. Initial D=0.5 plus seed-specific 0.03 noise; N=0; no future noise. Traction magnitude is 0.15+2.35D for all mobile cells, including followers. These uncalibrated constants are hypothesis choices, not measured proteins or kinetics. This replaces the original SOX2/SNAIL circuit only in the diagnostic, not the application.

Readouts: extension in +x still connected to parent; fraction of mobile cells attached; detached cells; components; tip heterogeneity; planar elongation. A narrow-extension candidate is a connected component of parent-attached cells with x>0, at least 3 cells, max x>=2, and y-span <=2. Persistence proxy: overlapping candidate IDs across three half-unit observations (one time unit), with intersection >=half the smaller set for consecutive pairs. This is a screening proxy, not a branching skeleton or bifurcation measure. Detached migrating cells must not count as sprouts.

Check fixed anchors, finite states, matched initial states, equal outer clocks, source hashes and timestep differences. Because propulsion is held within each outer interval, monitor W=U_pair-sum(f_active dot x) each inner step, rebasing each outer interval. Flag increases above 1e-9 times max(1,abs(old W),abs(new W)). Do not require global monotonicity across regulatory/polarity changes. Strong overlap and connectivity cutoff dependence are limitations to report.

## Results

All 24 runs completed: 4 conditions x 3 seeds x 2 inner resolutions. Nine
tests passed (5 for this assay plus 4 for its pair-force helper). All 600
saved states were finite, initial arrays matched between conditions and
resolutions, the 12 fixed cells remained exactly fixed, and source hashes
remained unchanged. Each run performed 240 outer updates; 138,240 total
inner updates were checked, with zero active-work-adjusted potential
increases above tolerance. This is an integration diagnostic, not evidence
that the chosen force law is biologically adequate.

At model time 12, means over 3 seeds at 32 substeps:

| Condition | Mobile cells attached to parent | Mean mobile x | Mean mobile tip activity | Mobile tip activity SD | Persistent narrow extensions |
| --- | --- | --- | --- | --- | --- |
| Neither | 100% | -4.344 | 0.315 | 0.0329 | 0 |
| Guidance only | 0% | 10.127 | 0.978 | 0.00343 | 0 |
| Competition only | 100% | -4.372 | 0.231 | 0.0188 | 0 |
| Both | 0% | 3.357 | 0.335 | 0.00484 | 0 |

The unguided arms retract into the fixed parent. Both guided arms compact
into a mobile cluster and detach: all 24 mobile cells are already fully
detached at the first half-unit saved observation for every seed. The
guided cluster travels farther without inhibition. No narrow-extension
candidate appears at any sampled time in any condition; therefore no
persistent candidate appears either. There is no evidence for bifurcation.
The maximum attached +x position at saved times is only 0.0502 model units,
inherited from initial jitter, not a newly growing extension.

Final attachment classifications are unchanged when the geometric cutoff
is changed from 1.2 to 1.1, 1.5 or 2.5, at both resolutions. Thus final
detachment is not an artifact of the chosen 1.2 contact threshold alone.
The exact unsampled breakage time remains unmeasured.

Final matched-cell position RMS differences between 16 and 32 steps span
0.00008-0.09858 model units; DLL4 RMSE spans approximately 0.000005-0.00603.
Some individual arrangements remain sensitive, but the qualitative
contraction/detachment failure repeats for all tested seeds/resolutions.
This does not prove general convergence or exclude brief unsampled events.

## Conclusions

1. Directional guidance changes migration, but the implemented forces do
   not maintain a parent-connected extension. Migration is not sprouting.
2. This particular contact-averaged inhibition circuit reduces activity
   but does not establish a strong tip-versus-follower pattern among mobile
   cells. In the combined arm their final activity is nearly uniform,
   mean 0.335 and within-mobile SD about 0.00484. The global SD also includes
   immobile parent cells and must not be mistaken for leader selection.
3. Adding both abstract mechanisms is insufficient under these parameter,
   geometry and boundary assumptions. This is a negative result for this
   implementation, not evidence against VEGF/Dll4/Notch biology.
4. Mechanical adequacy is unresolved even though inner integration passes
   the dissipation check. Maximum pair overlap reaches 0.876 model units
   for a nominal diameter of 1: severe interpenetration. Adhesion acts on
   every nearby pair up to 2.5, rather than explicit deformable junctions.
   This strongly compacting soft-cell regime limits biological conclusions.

We cannot attribute failure uniquely to one missing mechanism: the assay
changes the original regulatory model, supplies a prescribed field, fixes
one parent boundary, uses planar cell centers, omits proliferation and
matrix mechanics, and tests only one set of coefficients. Retaining
equal adhesion avoids promoting type sorting as a substitute for cohesive
tip-stalk junctions; this experiment does not test state-dependent junction
turnover or pair-compatibility together with guidance.

## Next research decisions

Before another combined sprouting trial, test two prerequisites separately:
(a) a static contact network should develop localized leader/follower
contrast under the proposed signaling circuit; (b) a small mechanically
loaded tissue should resist gross overlap and maintain appropriately
extensible junctions under directional traction. Vary repulsion, interaction
range and junction dynamics in explicit controlled comparisons rather than
silently tuning parameters until a branch-shaped image appears. Then
recombine and test growth/proliferation separately. These follow-ups have
not been run here.

Full ECM degradation, filopodia, lumen formation and blood flow are outside
this early-sprouting assay. No application source changes were made.

## Reproduction and provenance

```powershell
rtk proxy python tools/test_branching_sprouting_assay.py
rtk proxy python tools/test_branching_adhesion_assay.py
rtk proxy python tools/branching_sprouting_assay.py --output research/experiments/2026-11-07-guidance-tip-competition
rtk proxy python tools/branching_sprouting_report.py research/experiments/2026-11-07-guidance-tip-competition
```

Use a fresh output folder for reruns. Settings preserve equations' constants
and source/helper hashes. NPZ files preserve each initialization and all 25
half-unit observations per run. CSV files preserve individual timecourses,
final results, seed summaries, resolution differences and connectivity
sensitivity. `validation.json` records checks and trajectory hashes;
`sha256.json` covers this archive. Archived diagnostic scripts accompany
the data. The trajectory figure shows seed 42 at the finer resolution;
gray tracks/open circles show motion/start, filled colors show final tip
activity and red crosses identify fixed cells. Panels have independent
axis limits. Timecourse bands are seed min/max, not confidence intervals.

Biological motivation (these references do not validate our abstract rules):
[VEGF guidance, Gerhardt et al. 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC2172999/),
[dynamic tip competition, Jakobsson et al. 2010](https://doi.org/10.1038/ncb2103),
[VE-cadherin junction dynamics, Bentley et al. 2014](https://doi.org/10.1038/ncb2926).
