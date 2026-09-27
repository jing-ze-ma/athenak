---
name: cs-mhd-lowbeta-fix
description: THE FIX for the cubed-sphere MHD low-beta instability -- cb1afd28 + 1d7c4e5a, a per-cell HLLD->HLLE fallback below beta 0.5, default ON for cs. The original blow-up now runs 9x past where it died with dt intact. Dissipation was the answer; accuracy fixes were not -- BUT IT DOES NOT SAVE THE PRODUCTION dhj RUN, which still dies at 0.21 rotations (see cs-dhj-production-retry). True for the test problems only
metadata:
  type: project
---

**SCOPE, added 2026-09-03: this is a fix for the TEST PROBLEMS, and it is INERT in production.**
The production dhj run's MINIMUM plasma beta is **26.3**, against this switch's 0.5
threshold -- zero faces out of 786,432 can trigger it, so it never fires there at all. The
earlier note that it "does not save the production run, an 8.5% delay" is retracted: the
arms were two identical schemes. See [[cs-dhj-production-retry]].

Re-run of the cubed-sphere deep-hot-Jupiter production case with the fallback in
(`bench/cs_prod_mhd_gate`, job 11353535) still collapses: dt -> 0 and simulation time
freezes at t = 6.44e4 s = 0.211 rotations, against 5.93e4 s = 0.194 without it. An 8.5%
delay, same cycle (7600), same pathology. The history row at t = 6.1e4 is finite rather
than NaN and dt runs ~1.7x higher through the approach, so the fallback is doing real
work -- it is simply not sufficient. [[cs-mhd-dhj-blowup]] is STILL OPEN.
See [[cs-dhj-production-retry]] for the numbers.

**2026-09-03. Committed as cb1afd28, then 1d7c4e5a.** Closes the thread characterized in
[[cs-mhd-instability-characterized]] and reproduced by [[cs-mhd-minimal-reproducer]].

## The fix

`<mhd>/cs_lowbeta_fallback` -- the plasma beta below which the flux falls back **per CELL**
from HLLD to **HLLE**. **Defaults to 0.5 on the cubed sphere, 0 (off) on every other grid.**
`src/mhd/rsolvers/cs_lowbeta_fallback.hpp`, applied in all three sweeps of `mhd_fluxes.cpp`
AFTER the Riemann solve and BEFORE `GnomonicEquiangleFlux*` rotates into covariant slots.
The L/R states it uses are the RECONSTRUCTED ones, so it costs dissipation, not order.

Switching per CELL rather than per run is the point: the deep atmosphere the cs grid exists
for is high beta, subsonic and long-lived, and that is where HLLD's smaller dissipation on
contacts and Alfven waves is worth keeping. `mhd_fluxes.cpp` already makes exactly this
trade at a POLAR boundary (`do_pole`), for a whole row; this is the same trade, finer.

**HLLE, not LLF (1d7c4e5a).** cb1afd28 used LLF only because `SingleStateLLF_MHD` was the
only per-CELL solver in the tree (written for FOFC); that is an accident, not a reason.
`hlle_mhd_singlestate.hpp` is the missing twin, ideal and general EOS. Either stabilises
it -- **at REST they are IDENTICAL**, since with v = 0 and symmetric states HLLE's fan is
+/- c_f and degenerates to Rusanov, which is also the check that the new code is right
(per-region errors 2.054e-05 / 3.994e-05 / 6.662e-05 against LLF's ...6.663e-05). They
part company once the flow MOVES, and these atmospheres have fast zonal winds: LLF spreads
both sides at |v| + c_f regardless of direction, HLLE keeps the signed fan and reduces to
full upwinding when supersonic. Same stabilisation, less smearing. Wave speeds are the L/R
estimate, not the Roe average -- which is what `hlle_mhd.hpp` itself does on its
general-EOS path, and errs toward MORE dissipation, the safe direction for a fallback.

**NOTHING TO DO WITH FOFC**, which is a different mechanism (it redoes a failed cell AFTER
the update from donor-cell states) and remains a startup FATAL on cs
(`src/mesh/mesh.cpp:380`, from 9f79a031). Only a leaf header was ever shared. This fallback
runs INSIDE the flux kernel on the RECONSTRUCTED face states, BEFORE
`GnomonicEquiangleFlux*` rotates into covariant slots, so the seam exchange and the CT EMFs
see it as they see any other flux -- confirmed: seam flux mismatch stays at round-off,
6.2e-19 against a flux scale of 3.4e-07, with the fallback firing on every face.

## Why dissipation, and not any of the accuracy fixes

Only the Riemann solver changed the BEHAVIOUR. Everything that reduced the truncation
error changed the seed and not the growth:

    hlle (or llf) alone      dies -> saturates.  THE stabiliser
    ppm4                     nothing (and it is the WORST arm when paired with hlld)
    well-balanced source     nothing alone; 10x once hlle has removed the instability
    GS07 corner EMF          30x over baseline; nothing once hlle is in

HLLD is the least dissipative solver on exactly the modes involved: its Alfven and
slow-mode structure degenerates as the NORMAL field weakens, which on a curvilinear grid
carrying a strong TANGENTIAL field is the common case. The radial monopole -- tangential
field exactly zero -- is stable at the same beta on cs AND on spherical polar, which is the
positive control for that reading.

## Measured (unstratified iprob=8, KE per unit volume, dt_end/dt0)

                        baseline                    with the fallback
    beta .05, nx=16   1.10e-05  dt 0.46           4.86e-06  dt 0.99
    beta .05, nx=32   2.00e-06  dt 0.11           3.78e-07  dt 0.99, never reaches 1e-6
    beta .02, nx=16   1.03e-04  DEAD at t=1.58    3.16e-05  runs to t=3
    beta .02, nx=32   2.99e-04  dt collapsed      2.73e-06  dt 0.98

**The stratified reproducer, which died at t = 0.2013, is past t = 1.02 with dt at 98.8% of
its initial value** (nx=16; nx=32 likewise past t = 0.39).

## The gates it passes

* field-free null unchanged, 1.0248e-17;
* per-step convergence at beta = 0.2, where the fallback IS active, unchanged --
  1.86/2.52 interior, 1.82/2.65 at a cube vertex;
* **at beta = 2 the history file is BITWISE IDENTICAL with it on and off**, so it is a
  strict no-op above the threshold and no existing high-beta result moves;
* cost +7.8% in the WORST case, where every cell in the domain triggers.

## What to do next

1. Re-run the dhj production reproducer on this grid -- [[cs-mhd-dhj-blowup]] is the run
   that motivated all of it and has not yet been retried with the fallback.
2. Decide whether to raise the default threshold. 0.5 is a guess with a factor ~10 of
   margin over the measured onset (beta ~ 0.05); the dhj atmosphere runs to beta ~ 1e-4 at
   the top, so nearly the whole magnetised region will use LLF either way.
3. `<mhd>/cs_wellbalanced_src` (c3992145) is worth turning ON alongside it -- once the
   instability is gone it buys 10x -- but it costs 1.68x wall clock until its face triads
   are cached per (m,k,j). That caching is the obvious follow-up.
4. The `<mhd>/cs_gs07_emf` corner EMF (2566a9f7) is the more correct CT scheme and worth
   30x on its own, but adds nothing once the fallback is on. Keep it default OFF until it
   has a full convergence gate.
