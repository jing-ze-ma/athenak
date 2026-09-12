---
name: solar-convection-general-eos
description: "Ideal-vs-tabulated-EOS comparison of solar_convection: the ideal-gas p->e conversions that broke the first table run, where the runs live, and the offline EOS analysis tooling"
metadata: 
  node_type: memory
  type: project
  originSessionId: b25dc958-761c-420e-9c05-1aabb20efebb
  modified: 2026-08-15T16:16:23.268Z
---

Started 2026-08-14. Goal: run `solar_convection` with `general_eos = table` (H2 + H/He
Saha, **radiation OFF**) and compare against the ideal-gas control. See
[[general-eos-project]] and [[solar-convection-test]].

## RESUME HERE (2026-08-16)

**The atmospheric runaway is STILL OPEN.** Both the corrected sponge AND the taller box are
done and neither gives a steady atmosphere. The tall box buys a quasi-steady t~2000-15000
and then resumes. Next candidate is real radiative damping — see the next section.

Both 192^2 jobs COMPLETED (189934 ideal192 5h01m, 189932 table192 7h16m), 20 dumps each,
final dumps intact. The 64^2 and 192^2 comparisons are both DONE — see "Results".

**THE OPEN QUESTION: the table run's atmosphere is not in steady state.** Fraction of the
box with Mach > 1, versus time:

| t | table (64^2) | table192 | ideal / ideal192 |
|---|---|---|---|
| 5000 | 0.02% | 0.52% | 0 |
| 10000 | 0.01% | 1.51% | 0 |
| 15000 | 0.39% | 9.27% | 0 |
| 20000 | 4.80% | **28.0%** | 0 |

Monotonic growth at BOTH resolutions, ~5000 s faster at 192^2 (better-resolved waves =
less numerical damping). Total KE tracks it (table192 1.21e30 -> 1.56e30 over the last
third; ideal192 flat at ~1.4e30). So the 64^2 headline "v_z,rms 3.44 km/s at the top,
peak Mach 0.89" was a growing transient caught early, NOT a converged EOS effect. The
sub-photospheric comparison is unaffected — base and photosphere quantities are steady.

**The evidence points at the top boundary**, not at convectively driven acoustic waves:
1. The supersonic region appears FIRST AT THE LID and spreads DOWNWARD (t=10000: only
   z/zmax 0.82-0.95; t=20000: down to 0.445). Convective wave driving would start just
   above the photosphere and move up.
2. The horizontally averaged T INVERTS near the lid in every table run (192^2 at t=20000:
   3159 K at z/zmax 0.885 rising to 3663 K at the top), while ideal192 falls monotonically
   to 3525 K ~ `T_top_fix` = 3500. The table atmosphere wants to be far colder (min T
   ~2000-2400 K), so the 3500 K fixed-T lid is a hot, reflecting boundary sitting on an
   atmosphere ~100x more rarefied than the ideal run's.
NOT PROVEN — could still be genuine acoustic amplification needing a taller box or a
damping layer.

Note both ICs start with a nearly isothermal ~4860 K top, so `T_top_fix` = 3500 K is a
FORCED cooling of the top in both runs; the ideal atmosphere relaxes down to 3525 K and
sits consistently with it, the table one does not.

**THE LID TEST IS DONE, AND IT IS NEGATIVE — the lid TEMPERATURE is not the driver.**
Four 96x64^2 table runs, identical except `problem/T_top_fix` (jobs 189970-2, ~90 min each
on p.shared, all 21 dumps intact; binary pinned at `soleos/lidtest/athena`). Fraction of
the box with Mach > 1:

| t | 3500 K (`table`) | 4860 K (`lidtest/lid4860`) | 3000 K (`lid3000`) | 2500 K (`lid2500`) |
|---|---|---|---|---|
| 10000 | 0.01% | 0.00% | 0.24% | 0.12% |
| 14000 | 0.13% | 0.03% | 0.56% | 1.06% |
| 17000 | 1.61% | 1.66% | 3.75% | 1.91% |
| 20000 | 4.80% | 8.57% | 6.71% | 9.16% |

All four are quiet until t ~ 12000-15000 and then take off together; the final spread
4.8-9.2% has NO ordering in `T_top_fix`, and the onset time is the same in all four. The
4860 K run — where the lid is ~750 K HOTTER than the gas below it (lid 4681 K, min 4089 K)
and the ideal-gas-tuned mismatch is largest — grows if anything faster than the 3500 K
baseline. So the T inversion near the lid is a SYMPTOM, not the cause; my earlier
"hot reflecting lid" hypothesis is dead in its temperature form.

## TALLER BOX: DONE, AND IT DELAYS THE RUNAWAY WITHOUT REMOVING IT (2026-08-16)

Jobs 190139/190140 finished, 21 intact dumps each, t=20000.
`soleos/tall/{nosponge,sponge}`, binary `soleos/tall/athena` (md5 `88b24b76`, commit
`cabf4922`). Outputs `soleos/{lidtest,acoustic}_table_nosponge_nosponge_sponge.{txt,png}`.
144 x 64^2 over [0,3e8], dz unchanged, `z_ph_frac = 1/3` so the convection zone is identical
and the atmosphere goes from 8.7 to 16.7 pressure scale heights.

| E_atm t=11000 -> 20000 | dE/dt / F+ at t=20000 |
|---|---|
| table (2e8, old layer): 7.95e10 -> 3.97e11 (x5.0) | +0.007 |
| sponge2/nosponge (2e8, undamped): 1.54e11 -> 8.40e11 (x5.5) | +0.012 |
| **tall/nosponge (3e8): 8.13e10 -> 3.50e11 (x4.3)** | **+0.014** |
| **tall/sponge (3e8): 7.08e10 -> 1.91e11 (x2.7)** | **+0.007** |

**NOT STEADY.** dE/dt is still positive and tall/nosponge has the LARGEST value of all four
at t=20000. The prediction written before the runs held.

**What the tall box DOES buy: a long quasi-steady interval.** From t~2000 to t~15000
dE/dt/F+ oscillates about zero with repeated sign flips (-0.002,+0.001,-0.002,...), which
the short box never does after t=12000; the runaway restarts at t~16000. f(M>1) likewise
saturates at 18-24% from t=2000 to 18000 instead of growing monotonically — but that metric
is diluted by the added rarefied region and is weak evidence.

**The growth rate at a fixed height is UNCHANGED.** At z - z_ph = +0.5e8, Mach_rms goes
0.103 -> 0.325 (short) and 0.080 -> 0.260 (tall) between t=10000 and 20000: x3.2 vs x3.25.
The tall box starts lower and has more room before the lid; it does not slow amplification.

**CAVEAT that weakens the comparison: the tall runs' atmosphere is ~1.9x LESS DENSE than the
short runs' at every height above the photosphere** (rho0 at +0.25e8: 1.64e-7 vs 8.56e-8;
at +0.5e8: 1.34e-8 vs 7.06e-9; at +0.9e8: 2.6e-10 vs 1.39e-10). The sub-photospheric IC is
verified identical, so this is relaxation with the lid 1e8 further away, and/or the IC
integrator's N=10000 spanning 1.1x3e8 instead of 1.1x2e8. **Understand this before quoting
the tall runs quantitatively** — "more room" is currently confounded with "different
atmosphere".

Best of the four is `tall/sponge` (x2.7, +0.007, Mach 0.69 at +0.9e8 against 1.22
unsponged), and it is still growing.

**Remaining candidate: a real energy sink, i.e. radiative damping in the atmosphere.**
Geometry and velocity damping have both now been tried and neither closes the budget.

## SUPERSEDED: the taller box while it was running

### (original setup notes)
## TALLER BOX: RUNNING (2026-08-16, jobs 190139 / 190140)

`/orion/u/jinma/ATHENAK/soleos/tall/{nosponge,sponge}`, binary pinned at
`soleos/tall/athena` (md5 `88b24b76`, = commit `cabf4922`). p.exclusive, 1 node each,
32 MPI x 3 OMP, submitted 2026-08-16, ~80-90 min expected (1.5x the sponge2 runs).
**Analysis to run when they finish:** `lidtest.py` and `acoustic.py` with the run set
`table sponge2/nosponge tall/nosponge tall/sponge`.

Setup: 144 x 64^2 over x1 = [0, 3e8], dz UNCHANGED at 20.8 km. `problem/z_ph_frac = 1/3`
(the new parameter, commit `cabf4922`) puts the photosphere back at z = 1e8, so the
convection zone is identical to the 96 x 64^2 runs and **all** the extra height is
atmosphere: 16.7 pressure scale heights above the photosphere instead of 8.7.
`tall/sponge` uses `sponge_zbot = 0.8, sponge_c = 0.3`, which damps the top 600 km — the
same PHYSICAL thickness as `sponge2/sponge07` did in the short box — so box height and
sponge strength stay separated.

**Why the pgen needed changing first:** the IC shoots on p_top until tau = 2/3 lands at
`z_ph_frac` of the box, so before this parameter existed, raising `x1max` deepened the
convection zone as well and would have been two experiments at once.

**The exponential constraint, which limits how tall is worth trying:** above the
photosphere this atmosphere is nearly isothermal (T ~ 3900-4900 K) with H_p ~ 95-117 km,
so every extra ~100 km of box costs a factor e in density. The 2e8 reference run ends with
rho ~ 1.6e-10 at its top against `dfloor = 1e-11` — only ~2.8 scale heights of headroom.
The tall run therefore also needs `dfloor = 1e-18` and `eos_logd_min = -20` (table 401x551,
26 MB). Its IC top is rho = 5.95e-14 and it will evolve to ~1e-15/1e-16. Going much beyond
3e8 runs into the floor no matter what.

IC verified before submitting (the check from "The bug that invalidated the first run"):
base state T = 14711.1, rho = 1.69982e-5, p = 1.81474e7, grad_ad = 0.155278 and initial
dt = 5.395427e-01 all reproduce the reference to four digits; photosphere rho = 9.760e-7,
T = 5628 K against the reference's 9.756e-7 / 5627; HSE residual 3.8e-3 against the
reference's 4.3e-3.

**Prediction to test, and the reason for doubt:** the diagnosis is that the atmosphere has
no energy SINK — it is radiatively undamped (t_rad/t_ac ~ 11-17) and the shock heating that
inflates it has nowhere to go. A taller box adds room for the shocks to dissipate before
they reach the lid, but it does not add a sink, and the sponge (which IS a sink, and a
perfect no-reflection boundary) did not stop the accumulation. So the honest expectation is
that this delays the runaway rather than removing it. If so, the remaining candidate is
real radiative damping in the atmosphere, not geometry.

## THE CORRECTED 64^2 SPONGE COMPARISON IS DONE — AND THE SPONGE DOES NOT FIX IT (2026-08-16)

Jobs 190049-51 in `soleos/sponge2/{nosponge,sponge08,sponge07}` all finished to t=20000
with 21 intact dumps. Analysed with the working sponge (commits `3d1640e8` + `bfa566c2`).
Outputs: `soleos/{lidtest,acoustic}_table_nosponge_sponge08_sponge07.{txt,png}` (the tools
now tag by run set, commit `ae33271d`). `table` = the OLD always-on fixed tau = 50 s layer,
`nosponge` = a genuinely undamped top, which no earlier run had.

| at t=20000 | nosponge | sponge08 (zbot .8, c .1) | table (old 50 s layer) | sponge07 (zbot .7, c .3) |
|---|---|---|---|---|
| f(M>1) | 20.50% | 8.34% | 4.80% | 1.05% |
| max Mach_rms | 1.77 | 1.17 | 0.89 | 0.57 |
| E_atm t=11000 -> 20000 | 1.54e11 -> 8.40e11 (x5.5) | 1.14e11 -> 5.91e11 (x5.2) | 7.95e10 -> 3.97e11 (x5.0) | 1.05e11 -> 2.92e11 (**x2.8**) |
| dE/dt / F+ at t=20000 | +0.012 | +0.012 | +0.007 | +0.002 |
| mid-atm Mach, t=10000 -> 20000 | 0.154 -> 0.613 | 0.131 -> 0.413 | 0.109 -> 0.336 | 0.110 -> 0.375 |

**NOT SOLVED. Every configuration is still accumulating wave energy at t=20000** — dE/dt is
positive in all four and E_atm nearly triples even in the strongest sponge. The corrected
layer attenuates the runaway; it does not arrest it. This KILLS the earlier "sponge08 works,
use it" conclusion, which came from the garbage damping mask.

**f(M>1) is a compromised diagnostic for sponge runs and must not be quoted alone.**
sponge07 has the lowest supersonic fraction (1.05%) yet a LARGER mid-atmosphere Mach (0.375)
than the old-layer control (0.336): the sponge damps precisely the top cells where M>1 is
counted, while the mid-atmosphere keeps building underneath it. Use E_atm and dE/dt/F+.

Ordering is monotonic in damping strength for f(M>1) and max Mach, so the layer is doing
what it says. sponge08 being WORSE than the old fixed 50 s layer (8.3% vs 4.8%) is either
run-to-run scatter in an exponentially growing quantity — a ~1000 s shift in onset explains
a factor 2 — or the earlier estimate that the defaults are "~2x stronger at the topmost
cell" was wrong. Not resolved.

**Next, if this thread is picked up:** the remaining candidates are a genuinely taller box
(the sponge was meant to make that affordable and is not a substitute for it), a much
stronger/deeper sponge than zbot 0.7 / c 0.3, or real radiative damping in the atmosphere.
The sub-photospheric ideal-vs-table comparison is unaffected by any of this and stands.

## THE SPONGE WAS BROKEN — the 64^2 sponge results below are INVALID (2026-08-15)

See [[sponge-inert-at-192]], now resolved: the kernel read `pcoord->x1v`, a 1x1 placeholder
View on Cartesian meshes, so its `x1v > zs` test compared against out-of-bounds heap
memory. Inert at 192^2 (hence `table192sp` = `table192` bit-for-bit, and job 189980's
result is just a `table192` repeat); at 64^2 it damped whatever cells the neighbouring heap
happened to select, NOT the top 20%. **So the sponge08 / sponge07 table below characterises
a garbage damping mask, not the intended layer, and has to be re-derived.**

Fixed and verified at 96x192^2 (`soleos/spongefix/`): rms v_z ratio on/off is exactly 1.000
below z/zmax = 0.8, ramping quadratically to 0.393 at the top cell in 100 s.

**Also folded in (2026-08-15):** this pgen had a SECOND, unconditional damping layer in
`SourceFunc` (top 20%, fixed tau = 50 s, same ff^2 ramp, KE likewise discarded) that had
been in EVERY run including all the controls. It is now deleted and `problem/sponge`
defaults to **true**, so the parameterised layer is the only one. At the defaults
(zbot 0.8, c 0.1) it is ~2x stronger at the topmost cell than the old fixed 50 s layer.
`sponge = false` gives a genuinely undamped top, which NO earlier run had.

## SUPERSEDED (see above): the 64^2 sponge result, 2026-08-15

`problem/sponge` (commit `1f0cbc9a`, off by default) damps velocity over the top
(1 - `sponge_zbot`) of the box at rate `sponge_c`*cs/dz with a quadratic ramp, DISCARDING
the kinetic energy rather than thermalising it. Two 96x64^2 table runs vs the `table`
control (jobs 189974/189975, ~1h28m each, 21 intact dumps):

| | control | sponge08 (zbot .8, c .1) | sponge07 (zbot .7, c .3) |
|---|---|---|---|
| f(M>1) at t=20000 | 4.80% | ~0.1%, no trend | ~0.1%, no trend |
| E_atm t=11000 -> 20000 | 7.95e10 -> **3.97e11** (x5.0) | 1.88e11 -> 1.51e11 | 1.37e11 -> 1.37e11 |
| dE/dt / F+ after t=12000 | persistently +0.002..+0.009 | oscillates about 0 | oscillates about 0 |
| mid-atm Mach, t=10000 -> 20000 | 0.109 -> **0.336** | 0.085 -> 0.066 | 0.070 -> 0.061 |
| mid-atm F+ | 2.7e7 -> **2.8e8** | 1.3e7 -> 1.1e7 | 8.7e6 -> 7.8e6 |
| R_mid, late | **1.28-1.71** (standing waves) | < 0.70 | < 0.70 |
| totE drift | +1.88% | +1.37% | +1.14% |

(INVALID — the layer was damping garbage-selected cells; kept only as a record.)
**Both work; USE `sponge08`** — as effective while damping only the top 20%. The energy sink
is modest and IMPROVES the energy drift, because much of the control's "drift" was the
accumulating wave energy. Mass drift unchanged at ~+2.1% (the old boundary imbalance).

**This does not contradict the lid test.** Reflection never distinguished ideal from table
(R ~ 0.3-0.7 in both); the AMPLITUDE does. But once the table run goes nonlinear nothing
removes the energy, so it accumulates — the sponge supplies the missing sink. The control's
R_mid > 1 appears only LATE (early it is 0.18-0.80, like the sponge runs), so it is a
symptom of the shocked state, not its cause.

**Caveat:** the sponge is a numerical device. Its steady state is not proven to match a
genuinely taller box — the sponge runs' mid-atmosphere is ~40% less dense than the
control's. A taller-box run is the real validation, and the sponge makes it affordable.

## The diagnosis behind it: physical amplification, not a boundary artifact (2026-08-15)

`tools/acoustic.py` (new) settles it, and it rules out BOTH of the hypotheses the lid test
left standing. It splits each plane into mean + fluctuation and forms the up/down-going
acoustic fluxes `F+- = (rho0 cs0/4)<(w' +- p'/(rho0 cs0))^2>` (whose difference is exactly
`<p'w'>`, a built-in check), the stored wave energy above the photosphere, and the Spiegel
radiative time `t_rad = cv/(16 kappa sigma T^3)` against `t_ac = H_p/cs0`. `cv` comes from
the run's OWN eos by finite differencing `temperature(rho, e)`. Run `python3 acoustic.py
[run ...]` -> `soleos/acoustic.txt`, `plots/acoustic.png`.

- **Reflection is NOT the differentiator.** `R = F-/F+` is ~0.3-0.7 at the top in BOTH
  runs, noisy, with no systematic difference. Same BC, same partial reflection.
- **Radiative damping is NOT the differentiator either.** `t_rad/t_ac` is ~11-17 at
  mid-atmosphere in BOTH — both atmospheres are radiatively undamped.
- **The acoustic INPUT is the same.** `F+` crossing the photosphere is ~5-10e9 erg/cm^2/s
  in both, and steady in time in both.
- **What differs is the AMPLITUDE that same luminosity reaches.** At mid-atmosphere the
  table run is ~10x thinner (rho0 1.09e-8 vs 9.73e-8) with a 34% smaller `cs`, so at
  t=10000 it carries a SMALLER flux (2.65e7 vs 7.37e7) at a 1.8x LARGER `w'_rms`
  (7.3e4 vs 4.0e4) and 2.7x larger Mach (0.109 vs 0.040) — that is `F ~ rho0 cs0 <w'^2>`
  read backwards. The ideal run stays linear forever; the table run reaches shock
  amplitude, shock heating inflates the atmosphere, and more flux gets through: by
  t=20000 `F+` at mid-atmosphere has grown 10x to 2.84e8 (3x the ideal run) at Mach 0.34.
- The rarefaction itself is real EOS physics — H2 formation raises mu 1.15 -> 1.29 and
  halves the scale height.

**Conclusion: the runaway is physical acoustic amplification, and the box is too short and
too undamped to host it.** Not a bug, but the atmospheric statistics of every run so far
(v_z,rms at the top, peak Mach) must NOT be quoted as converged. To get a steady
atmosphere: a taller box with a sponge/damping layer near the top, or real radiative
damping. The sub-photospheric comparison stands unaffected.

Energy budget backing it: `dE_atm/dt` oscillates around zero for the ideal run at all
times, and is persistently POSITIVE for the table run from t=12000 on (+0.5 to +0.9% of
the photospheric acoustic luminosity per unit time), i.e. steady accumulation.

Analysis: `cd soleos/tools && python3 lidtest.py` -> `soleos/lidtest.txt` and
`plots/lidtest.png` (f(M>1)(t), max Mach_rms(t), near-lid `<T>(z)`). Takes run names.

**COMMITTED 2026-08-15** on `general-eos`: `88ace084` "Make solar_convection EOS-agnostic"
(the whole fix below) and `6ccac245` "make the top lid temperature a runtime parameter"
(`problem/T_top_fix`, default 3500 = unchanged behaviour). Still unfixed: the same
ideal-gas p<->e bug in `cooling_convection.cpp` (10 sites).

**Truncated-dump guard: DONE.** New `tools/dumps.py` — a dump more than 64 KB below the
run's LARGEST dump is truncated; `dump_files()` returns the full sorted list with those
entries replaced by `None` so frame numbers keep matching file numbers, and prints what it
skipped. `compare.py` / `slices.py` / `tau_surface.py` all use it. It flags exactly the
three known-corrupt files (`ideal192` 00005+00006, `table192` 00004). In `dumps.resolve`,
a NEGATIVE index counts back through the intact dumps only, so -1 is always usable.
`compare.py` is now `compare.py [idx] [ideal_run] [table_run]`, and tags its outputs
`_<ideal>_<table>` unless the pair is the default `ideal`/`table`, so the 64^2 figures are
not overwritten.

## Results: the 192^2 comparison (t=20000, `summary_ideal192_table192.txt`)

Below the photosphere it reproduces the 64^2 result closely: base T 13560/14630 K, base
Gamma_1 1.667/1.253, grad_ad 0.400/0.156, photosphere z/zmax 0.469 vs 0.504, 42 vs 45
superadiabatic cells. Photosphere T 5472 (ideal) vs 4879 K (table). Mass/totE drift:
ideal192 +7.87%/+6.64% (WORSE than 64^2's +5.31%), table192 +2.60%/+1.92% (same as 64^2).
The atmospheric numbers (v_z,rms at top 3.55 km/s, peak Mach 1.98) are NOT converged —
see "RESUME HERE".

## Results: the 64^2 ideal-vs-table comparison (VALID, t=20000)

| | ideal (mu=0.602 fixed) | table |
|---|---|---|
| base T | 13520 K | 14640 K |
| base rho | 2.354e-6 | 1.698e-5 |
| base p | 4.389e6 | 1.799e7 |
| base mu | 0.6025 | 1.149 |
| base Gamma_1 | 1.667 | 1.253 |
| base grad_ad | 0.400 | 0.156 |
| photosphere z/zmax | 0.457 | 0.480 |
| photosphere T | 5486 K | 5278 K |
| photosphere Gamma_1 / grad_ad | 1.667 / 0.400 | 1.549 / 0.340 |
| superadiabatic cells (of 96) | 42 | 44 |
| v_z,rms base / photosphere / top | 1.88 / 1.69 / 0.56 km/s | 0.81 / 1.53 / **3.44** km/s |
| peak Mach | 0.26 | **0.89** |
| mass / totE drift | +5.31% / +2.72% | +2.59% / +1.88% |

Headlines: the two runs now agree STRUCTURALLY (photosphere at nearly the same height,
same number of unstable cells, tau profiles crossing at the same place) — the differences
are real EOS physics, not a broken setup. Gamma_1 dips to 1.25 and grad_ad to 0.156 in the
H-ionization zone. The striking difference is the ATMOSPHERE: higher mu (H2 forms, mu
1.15 -> 1.29) gives a ~2x smaller scale height, so the top is ~100x more rarefied
(1e-10 vs 1e-8), waves amplify as rho^-1/2, and v_z,rms reaches 3.4 km/s at Mach 0.89
versus 0.56 km/s / Mach 0.26 for the ideal gas. Note the table run is BETTER conserved
than the ideal control; the residual linear ~+0.13%/1000 s mass drift is the pgen's
long-standing boundary imbalance (see [[solar-convection-test]]), not EOS related.

## The bug that invalidated the first run (FIXED, uncommitted)

The table run lost **56% of its mass** and reached Mach 1.8. Root cause: the pgen's IC
stored the internal energy as `p*igm1`, the **ideal-gas** conversion, while `get_wb_eos`
correctly returned the EOS-consistent `(rho,p)`. Under the table EOS `e` also carries H2
dissociation + ionization latent heat, so `(gamma-1)*e` overestimates p by ~4x. The IC was
therefore **75% out of hydrostatic balance at t=0** and collapsed from the first step.

**How it was caught, and the lesson:** the offline `eoslib` inversion of the t=0 snapshot
returned T=4580 K where the code had *printed* base T=14712 K. The printed base state came
from the correct `get_wb_eos_arr` integration; the state actually WRITTEN to memory was
different. The decisive test was not comparing numbers but **checking
`(dp/dz + rho g)/(rho g)` on the t=0 snapshot**: ideal 3.5e-4, table 0.75. Do this on any
new stratified IC before trusting a run — a uniform residual is a smoking gun, and
`1 - 1/4.02 = 0.751` tied it straight back to the p/e ratio.

Fixed in `src/pgen/solar_convection.cpp` (UNCOMMITTED), all via `pgen_eos_utils.hpp`:
1. **IC** (~line 292/298): `p*igm1` -> `EintFromP(eos, igm1, den, p)`. Uses the PERTURBED
   density `den`, not `denwb`, so the 1% noise stays an entropy perturbation at fixed p.
2. **Well-balanced background**, 8 sites: the 6 `w0facewb_x{1,2,3}f` faces and
   `u0wb`/`w0wb` — `pwb*igm1` -> `EintFromP(eos, igm1, denwb, pwb)`.
3. **Bottom ghost hydrostatic extrapolation** (~line 904, a SEPARATE bug): the exponent
   `rho_is/eint_is*igm1` is `rho/p` only for an ideal gas. General branch now goes
   `PresFromEint` -> `TempKelvin` -> `p_g = p_is exp(-(rho/p) dphi)` -> `DensFromPT` ->
   `EintFromDensT` (constant specific e is constant T only for an ideal gas). On its own
   this made the ghosts ~4% under-dense and drained the box at ~2e-2 g/cm^2/s, which
   matched the observed drain — a real second bug, just dwarfed by the IC one.

Verified after the fix: table base T=14649 K, p=1.767e7 (code prints 14711.5 / 1.81481e7
at z=0; the snapshot is the cell centre at dz/2), **HSE residual 0.75 -> 4.3e-3**. Ideal
unchanged at 9.724e-4.

**`src/pgen/cooling_convection.cpp` has the SAME bug at 10 sites** — not fixed, not asked
for. Same one-line-per-site `EintFromP` treatment would do it. `hotbubble`, `convection`
and `hse_atm` are clean (0 sites).

## The tools now live IN THE REPO (2026-08-15, commit 8220c551)

`athenak/tools/solar_convection/` — `dumps.py eoslib.py compare.py slices.py
tau_surface.py lidtest.py acoustic.py eos_dump.cpp README.md`, committed on `general-eos`.
`soleos/tools/*` are now SYMLINKS to those files, so `cd soleos/tools && python3 ...`
still works exactly as before and there is one source of truth — edit either path.
`eos_grid.bin` (38 MB) and the `eos_dump` binary stay in `soleos/tools`, uncommitted.

Paths are no longer hard-coded: `ROOT = os.environ.get('SOLEOS_ROOT', <the old path>)`,
and `bin_convert` now comes from the repo's own `vis/python/` located via
`os.path.realpath(__file__)` — realpath, NOT abspath, because the symlinks would otherwise
resolve `../../vis/python` to `/orion/u/jinma/ATHENAK/vis/python`. Verified by rerunning
`compare.py -1 ideal table`, which reproduces the stored 64^2 summary exactly. Also added
a `__pycache__/` rule to `.gitignore` (the repo had none).

## Layout (all outside the repo and outside [[run-directory-untouchable]])

```
/orion/u/jinma/ATHENAK/soleos/
  athena              pinned binary (PROBLEM=solar_convection); REBUILT with the fix
  ideal/              64^2, job 189754, VALID, done (590 s wall)
  table/              64^2, job 189921, VALID, done (5150 s wall)
  ideal192/ table192/ 96x192x192, jobs 189934 / 189932 — the resolution upgrade
  table_brokenIC/     job 189755, the invalid run — keep for the writeup (164 MB)
  plots/              profiles/granulation/history .png, slices/<run>/, tau_surface/<run>/
  tools/  eos_dump.cpp  eoslib.py  compare.py  slices.py  tau_surface.py  eos_grid.bin
```
64^2 grid is 96 (vertical) x 64 x 64, box 2e8 cm, tlim=20000, 16 MPI x 4 OMP: ideal 590 s,
table 5150 s (the table EOS costs ~8x per zone-cycle). The broken run took 2.4 h and its
dt collapsed from 1.10 to ~0.2 — that was the collapse, not the EOS cost.

**192^2 = the historically best resolution** (job 179718, `run/sun_test/`): 96 x 192 x 192,
CUBIC 20.8 km cells, 3.54M zones. 64^2 is 3x anisotropic horizontally (62.5 km vs 20.8 km)
and under-resolves the intergranular lanes to 2-5 cells — adequate for the thermodynamic
comparison, not for granulation morphology. Meshblock 96x16x16 -> 12x12 = **144 blocks**
(nx1 MUST equal mesh nx1: the two-stream RT sweeps whole vertical columns), so MPI ranks
must divide 144.

**Cluster facts worth reusing** (orion, 2026-08-14): a node is **112 PHYSICAL cores**
(2 sockets x 56, ThreadsPerCore=2, so `sinfo` reports 224). With
`CR_CORE_MEMORY + CR_ONE_TASK_PER_CORE`, `--cpus-per-task` counts PHYSICAL cores and Slurm
doubles the allocation (`ReqTRES=cpu=64 -> AllocTRES=cpu=128`). p.exclusive queue depth is
misleading: of 390 pending jobs only ~2 waited on Resources, the rest on Dependency /
AssocMaxJobsLimit. QOS `normal` caps a USER at 30 nodes, no per-job limit. `sbatch
--test-only` gives real start estimates — 1-2 nodes started in 20 min, 4 nodes not for
10.5 h, so 2 nodes beat 4 by ~8.5 h wall.
The two inputs are IDENTICAL except `<hydro>` and `<units>`, so every difference is the
EOS. Both use `tfloor` at the SAME physical 200 K (2.76e10 ideal, 1.663e10 table).

## Base states at t=0 (AFTER the fix; the broken run's numbers were meaningless)

| | ideal (mu=0.602 fixed) | table |
|---|---|---|
| base T | 13718 K | 14712 K |
| base rho | 2.28e-6 | 1.70e-5 (7.5x denser) |
| base p | 4.32e6 | 1.81e7 |
| base grad_ad | 0.400 | 0.155 |
| initial dt | 0.35 s | 0.54 s |

## Offline EOS analysis tooling (the part worth reusing)

**AthenaK cannot output temperature for Newtonian hydro** — there is no `hydro_t` in
`var_choice` (`src/outputs/outputs.hpp`) — and the tabulated EOS lives in a `DvceArray`, so
it cannot be called from host or Python. Solution:

- `tools/eos_dump.cpp` — ~60-line standalone program that `#include`s the code's OWN
  `src/eos/eos_composition.hpp` and dumps a (log10 rho, log10 T) grid of
  `log10 e_spec, log10 p_spec, mu, xh2, xhii` to raw float64. No Kokkos needed:
  `g++ -O2 -I<athenak>/src -o eos_dump eos_dump.cpp`. Physics cannot drift from the code.
- `tools/eoslib.py` — `TableEOS` (interpolates that grid; `temperature()` inverts by
  60-step vectorised bisection on log10 T) and `IdealEOS` behind ONE interface.
  **API: every accessor takes `(rho, e_specific)`, not (rho,T) and not energy density.**
  Derivatives (chi_rho, chi_T, cv, Gamma_1, grad_ad) are `np.gradient` on the log-log grid.
  Also mirrors the pgen's `get_kapr` so tau matches.
- Validated against the code: at the table base state it returns p = 1.8148413e7
  (code printed 1.81481e+07) and grad_ad = 0.15529 (code 0.155273).
- `.bin` output is SINGLE precision, so agreement beyond ~1e-6 is not meaningful.
- `tools/tau_surface.py <run> [first] [last]` — granulation maps on the CORRUGATED
  tau=2/3 surface (T and v_z), into `soleos/plots/tau_surface/<run>/`. Adapted from
  `run/plot_tau_surface.py`, which is hardcoded to `sun_test/`, writes stray `.athdf`
  next to the data, and computes T with the ideal-gas formula. This copy reads the binary
  directly and takes T from `eoslib`, so it is valid for `general_eos=table` too.
  Colour limits default to the original 5200-5700 K / +-6 km/s; override with the
  `TMIN/TMAX/VMIN/VMAX` env vars.
- `tools/slices.py <run> [first] [last]` — vertical (x3-x1) 2D slices through mid-x2:
  vz, T, log10 dtau, Mach, into `soleos/plots/slices/<run>/<panel>/`. Adapted from
  `run/plotsun.py` the same way, and it also FIXES that script's stale opacity (plotsun
  used C=1e2 with no floor and no Kramers, which memory already flagged as untrustworthy;
  this uses `eoslib.get_kapr`). Colour limits kept at plotsun's values for comparability
  (vz +-2 km/s SATURATES — real velocities reach -10/+7 km/s); override with
  `VZLIM/TMAX/MACHMAX`.

## Notes / gotchas

- `eos_radiation = false` is the DEFAULT, so "without radiation" needs no extra flag.
- Table grid narrowed to `logd [-14,0]`, `logt [1.5,7]` (281 x 551, 14 MB).
- `<units>` is `length_cgs = mass_cgs = time_cgs = mu = 1`, i.e. code units ARE cgs.
  `T_code = 8.3145e7 * T_K`.
- The CO5BOLD bottom relaxation is **mass-neutral by construction** — step 4
  (`drho4 = (srho0-srho2)/sN`) restores the bottom plane's mean density, and step 5
  removes the mean vertical velocity. So a plane-averaged `rho*v_z` measured at the bottom
  CELL CENTRE is ~0 by design and tells you nothing about the leak; the leak is the
  Riemann flux through the `is-1/2` FACE, set by the ghosts.
- `Teff = 5778`, `T_top_fix = 3500`, opacity `C = 20` and the (disabled) bottom heating
  layer are unchanged and were tuned for the IDEAL run. Do not "fix" them before looking.
