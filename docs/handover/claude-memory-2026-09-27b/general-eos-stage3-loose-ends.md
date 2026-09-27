---
name: general-eos-stage3-loose-ends
description: "How the last three TODO(stage3)s were closed: entropy floor, the EnergyFromPressure re-solve, and the well-balanced background's temperature hand-off"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2cb708dd-0a5f-46c0-b912-dd406f996581
  modified: 2026-08-13T15:14:46.060Z
---

Commit `c94504f7` on `general-eos` (2026-08-13). `grep -rn "TODO(stage3)" src/` is now
empty. None of the three changes a gamma-law answer. See [[general-eos-project]].

**1. Entropy floor.** The general ConsToPrim skipped it, so a run setting `<fluid>/sfloor`
silently lost the floor the ideal path applies. Now `EOS_Data::ApplyEntropyFloor()`,
reproducing the ideal arithmetic verbatim — INCLUDING upstream's quirk of adjusting `w.e`
without correcting `u.e`, which the pressure and temperature floors DO correct. Under
`general_eos = table` it is a no-op and `BuildGeneralEOS()` **refuses at startup** if
sfloor is above its `FLT_MIN` default: sfloor floors `p/d^gamma`, not an invariant of a
tabulated EOS. Same "refuse rather than fake" rule as the Roe solver.

**2. `EnergyFromPressure(d,p,t)`** — a three-argument form returning the T the inversion
found, so a floored cell does not re-solve. The ideal branch computes `e` first and then
`t` FROM that e (not the algebraically equal `p/d`), which is what keeps it bitwise.

**3. The well-balanced background.** `WBBackgroundStencil()` did five `Pressure(d,e)`
calls = five temperature solves. `WBAdvance()` now returns its temperature (every branch
has one free: isodensity from the p-inversion, isothermal IS t_ref, isentropic from
`WBEnergyFromEnthalpy`'s last iterate) → **one solve per stencil instead of five**.

  **The TODO's own suggestion — hand back the PRESSURE — is WRONG, and this is the reason
  to remember.** The ideal-gas `getWBq0()` (`hydro.hpp` ~line 897) derives the pressure
  channel as `q0 * gm1` from the SAME `e` it reconstructs. So the general path must also
  compute p from e or it stops reproducing the ideal path bit for bit. A TEMPERATURE is
  safe precisely because `Pressure(d,e,t)` **ignores t** under a gamma law — gamma mode
  literally cannot detect the change.

**Verification method worth reusing.** For a floor, comparing ideal vs general directly is
weak (the baseline gap is ~1e-14 anyway). Instead compare the *effect*: run
ideal-with-floor minus ideal-without, general-with-floor minus general-without, and
require the two DELTAS to agree. Got O(4) hydro / O(13) MHD effects agreeing to 7.0e-15 /
5.2e-15 relative. Needs `sfloor` chosen to actually bite — for Sod (d=1,p=1 | d=0.125,
p=0.1, gamma=5/3), `K = p/d^gamma` runs 1 to 3.2, so `sfloor = 1.5` works.

**Regression evidence:** all 16 cases of `eos_compare.py` still reproduce
`eos_compare_baseline.txt` exactly; `wbs_plain`/`wbs_perturb`/`wb_iso`/`wb_isen` (needs
`PROBLEM=hse_atm`) and `dhj` (needs `deep_hot_jupiter_rt`) unchanged to the last digit,
with `wb_isen` and `dhj` at exactly zero. Test suites `hydro` + `mhd` pass except the
pre-existing `mhd_linwave` threshold ([[branch-preexisting-breakage]]).

**Gotcha:** a floor parameter cannot be set on the command line unless the input file
already defines it — AthenaK exits fatally otherwise. Edit the input, don't override.
