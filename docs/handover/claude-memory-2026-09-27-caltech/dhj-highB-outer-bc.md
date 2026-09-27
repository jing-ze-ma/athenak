---
name: dhj-highB-outer-bc
description: "FIXED (bc6b5774 then 767fff43): the outer-x1 Maxwell term scales as B^2/(rho g dr); it is now applied as a clamped effective gravity inside the hydrostatic solve, dfloor-independent 3-300 G. Includes the measured STS cost/accuracy A/B and the general-EOS cost factor."
metadata:
  type: project
---

Found 2026-08-16 while answering "make the timestep not diffusion-limited and let the code
run at 10 G". The user pushed back on my first (wrong) diagnosis and was right.

## Issue 1 -- diffusion-limited dt: ALREADY SOLVED, no new code

`<mhd>/use_rkg_sts = true` (RKG super-time-stepping) is fully wired and WORKS. Measured on
the dhj 64x64x128 polar grid with `ohmic_resistivity = constant`:

| eta | STS off | STS on |
|---|---|---|
| 1e13 | dt = 4.68 s | dt = 13.41 s (= CFL) |
| 1e15 | (0.047 s) | dt = 13.41 s (= CFL) |

### STS costs MORE per step — measured properly 2026-08-16 (the ~75x figure I first quoted
### was from a 2-cycle login-node sample and is WRONG; it is 45x)

RKG2, alpha = 0.5, `s = ceil(-alpha + sqrt((alpha+1)^2 + tau(3+2alpha)))` with
tau = dt_hydro/dt_diff, **floored at s = 3**. Benchmark, 100 cycles, 16x7 on p.shared,
dhj 64x64x128, constant eta, general+table EOS. Cost = wall seconds per SIMULATED second
(baseline with no resistivity at all = 0.00635):

| eta | tau | s | STS off | STS on | gain |
|---|---|---|---|---|---|
| 1e11 | 0.03 | 3 | 0.0060 | 0.0084 | **0.71x (STS 41% SLOWER)** |
| 1e12 | 0.29 | 3 | 0.0060 | 0.0087 | **0.69x** |
| 1e13 | 2.87 | 4 | 0.0174 | 0.0097 | 1.79x |
| 1e14 | 28.7 | 11 | 0.169 | 0.0159 | 10.7x |
| 1e15 | 287 | 34 | 1.674 | 0.0372 | 45x |

**Break-even is tau ~ 1.5-2.** One STS sub-stage costs ~0.14 of a full RK2 hydro step here,
so `cost_on/cost_off = (1 + 0.14 s)/max(1,tau)` reproduces every row above to <10%.

**ACCURACY COST, and it is not small:** `driver.cpp:406-421` runs the RKG super-step AFTER
the complete RK2 ideal step, once per cycle — **Lie (first-order) operator splitting**. With
STS OFF the resistive fluxes and EMFs are tasks INSIDE the RK2 stages, i.e. second-order
coupled. So turning STS on demotes the flow-field coupling from 2nd to 1st order in dt. Not
measured; do an A/B on a short run before trusting it for production.

Practical rule: watch dt. At the hydro CFL value -> STS is a pure ~40% loss. More than ~2x
below it -> STS wins. At the dhj `max_eta = 1e12` diffusion never binds, so **leave STS
off**; only turn it on together with max_eta = 1e13/1e14 if the tabulated x_e (rather than
the cap) is wanted over the cold part of the domain — see [[resistivity-perna-uhj]].

## Issue 2 -- the high-B crash is the OUTER BOUNDARY, not the Alfven speed

**Reproduction** (400 cycles, `problem/bbot` and `mhd/dfloor` swept):

| bbot | dfloor 7.26e-12 | dfloor 7.26e-15 |
|---|---|---|
| 10 G | dt 13.4 -> 13.0 | identical |
| 30 G | 13.4 -> 3.9 | identical |
| 100 G | 6.3 -> 3.3 | **dt collapses to 9.4e-7 s in 50 cycles** |

**Localisation.** Dump the broken case: at t=25.0 everything is benign (min rho 3.9e-10,
max v_A 1.9e6, max|v| 3.3e5). ONE step later min rho = 7.5e-14 and max|v| = 1.7e13 cm/s.
97% of the cells with |v| > 1e7 sit at i >= 62 of 0..63 — the last two ACTIVE radial
shells — and only 4% are near the poles. Shell i=62 alone drops to rho = 7.5e-14 while
i=61 and i=63 stay at ~3e-10. That is a boundary signature, not an evacuation.

**Mechanism.** `deep_hot_jupiter_rt.cpp` outer-x1 user BC (~line 1100) extrapolates
hydrostatically and then applies a MAGNETIC correction:
```
e0_ip = e0_hyd - e_i*dM1mag/rho_i/grav_acc;   if (e0_ip < 0.0) e0_ip = e_i;
rho0_ip = rho0_hyd*(e0_ip/e0_hyd);
```
`dM1mag` is div(Maxwell stress) + curvature, so the correction relative to the hydrostatic
term is `(B^2/dr)/(rho g)`. At the outermost shell (rho = 4.02e-10, dr = 5.625e7, g = 942):

| bbot | 3 G | 10 G | 30 G | 100 G | 300 G |
|---|---|---|---|---|---|
| ratio | 0.030 | 0.333 | 3.00 | 33.3 | 300 |

which reproduces the sweep exactly. Only `e0_ip < 0` is guarded; an O(1) or larger positive
swing goes straight into the ghost density. **`dfloor` was limiting how far rho could
crash, hence how big m/rho became — that is the whole reason a low floor "breaks" it.**

**FIXED 2026-08-16, commit `bc6b5774`:** `problem/bc_outer_maxwell`, DEFAULT FALSE, drops
the term (the cell-centred ghost field sharing that branch is still filled -- the rest of
the boundary needs it). Verified over 400 cycles: the result is then INDEPENDENT of dfloor
to all printed digits at bbot = 3, 10, 30, 100 and 300 G against dfloor 7.26e-12 and
7.26e-15. 10 G costs almost nothing (dt 13.33 vs 13.41 at 3 G); 300 G runs at dt 1.9 s,
which is the honest Alfven CFL of that field. At 3 G over 1000 cycles, dropping it moves
mass by 7.6e-6 and total energy by 9.4e-6, with the wave energies decorrelating at the
percent level (chaotic, not systematic).

## The Alfven-speed ceiling: built, worked, then DELETED at the user's request

`<mhd>/va_ceiling` floored density at `B^2/va^2` inside `SingleC2P_{Ideal,General}MHD`. It
did cure the collapse (va = 1e7 with dfloor = 7.26e-15 reproduced the 7.26e-12 reference to
6 digits, and 300 G ran), and it does NOT break div(B) — it never touches B, only the
density floor, so the face fields and CT are bitwise unaffected; what it violates is mass
conservation, exactly as `dfloor` already does. But it treats the symptom, so it was
reverted once the BC was identified. Do not re-add it without a reason the BC fix does not
cover.

**Method note worth keeping:** the thing that settled this was dumping the state one step
before and one step after the collapse and histogramming WHERE the bad cells are by index.
"97% at i >= 62 of 0..63, 4% near the poles" is what turned a plausible bulk explanation
into a boundary one.


## 2026-08-16 (later): the BC done PROPERLY, and the STS A/B

**`767fff43` supersedes the drop.** `bc_outer_maxwell` now DEFAULTS TRUE and means "include
the magnetic force, bounded": `dphi -> dphi*(1 - dM1mag/(rho|g|))` applied INSIDE the
hydrostatic solve (so inside WBAdvance under a general EOS), with the factor clamped to
[-1, 0.9]. Same to first order, but it changes the ghost's SCALE HEIGHT instead of swinging
its density linearly through zero, and rho/e now come from one solve instead of rho being
rescaled by e's relative shift. Upper clamp 0.9 = the force-free limit; above 1 the force is
net outward, which an outward-decaying ghost cannot represent.

400-cycle sweep, dfloor 7.26e-12 vs 7.26e-15, identical to every printed digit at all five
field strengths. dt at cycle 400, dropped vs effective-gravity: 3 G 13.40571/13.40571
(identical), 10 G 13.32741/13.22936, 30 G 5.792189/5.206927, 100 G 3.567449/3.718875,
300 G 1.917274/**2.138728**. At 300 G the magnetically supported ghost gives a LARGER dt
than pretending the field is absent.

## STS A/B: the splitting error is small (measured, not argued)

Constant eta so the cap binds everywhere, both to the same simulated time t=1500, bbot=3.
STS off = resistive fluxes inside the RK2 stages; STS on = Lie-split super-step after it.

| | eta=1e13 | eta=1e14 |
|---|---|---|
| 1-ME (dominant) | identical to the 6 digits the .hst carries | identical |
| 2-ME | 1.1e-5 | 1.5e-5 |
| 3-ME (1e-5 of 1-ME) | 2.1e-3 | 1.1e-2 |
| kinetic energies | 2.7e-4 .. 5.1e-3 | 3.8e-4 .. 7.6e-3 |
| mass, tot-E | identical | identical |

The DECAY of 1-ME (4.76524e10 -> 4.70907e10 at 1e14) is the same with and without. So the
thing the resistivity exists to do is captured; the percent-level differences are in
components that are 1e-5 of the total. **Verdict: at eta <= 1e14, STS on is worth 10.7x.**
Re-check if the resistivity ever shapes the field rather than just damping it.

**The viper input now runs `max_eta = 1.0e14` + `use_rkg_sts = true` + `bc_outer_maxwell =
true`.** Revert to `1e12` / off to reproduce the old runs.

## Cost of the general EOS on this problem (2000 cycles, 16x7)

ideal+perna vs general+table+eos, wall seconds per SIMULATED second:
bbot=3: 2.48e-3 -> 8.46e-3 (**3.4x**); bbot=10: 3.85e-3 -> 1.11e-2 (**2.9x**). Per CYCLE it
is 4.4x, partly clawed back because the table EOS gives a larger dt (7.64 vs 6.33 s at 3 G).
