---
name: general-eos-audit
description: 2026-09-04 audit of the general/tabulated EOS for bugs that could hurt sp or cs runs -- the (rho,e)->T inversion is BRACKETED, monotone and warm-start-independent (cold start = bracket midpoint), the radial BC inverts with the table, resistivity's T is consistent; no new live bug. The recorded issues are all FIXED. Two caveats: off-table states are silently log-linear-extrapolated, and ghost/active copies can differ by logtol = 1e-13
metadata:
  type: project
---

Asked: "is it possible that there is any bug involving general eos?" Checked with file:line.

**The inversion is branch-safe.** `EosTable::SolveLog` (`eos_table.hpp` ~300-360) is a
safeguarded Newton with a BRACKET on [ymin-3, ymax+3] and a bisection fallback whenever the
step leaves the bracket or fails to halve it -- the comment names the ionization-front
ping-pong it exists to catch. All three targets are monotone in their unknown, so (rho,e)->T
has ONE root. `SolveTemperature` maps `tguess <= 0` to a cold start at the bracket midpoint,
so a fresh, zero-filled ghost `wtemp` (Kokkos::realloc zero-fills) is handled, and
`TGuess()` returns -1 when the view is empty (ideal gas). **The warm start changes speed,
not the answer**, to `logtol = 1e-13`. The recorded 0.35-dex "wrong branch"
([[eos-inversion-nan-trap]]) was a NAIVE post-processing root find, not this solver.

**Consistent uses**: the outer-x1 user BC builds the ghost energy with the table
(`EintFromP -> EnergyFromPressure`, warm-started from the column's own solved T at ie,
`e0_ip < 0` guarded); resistivity's electron fraction takes T = tcode*temp_cgs from wtemp,
which RaiseVelMHD refreshes over ghosts too; the pgen's `EintFromP/PresFromEint` dispatch
to the table under IsGeneral().

**Already fixed, do not re-find**: stale p/Gamma_1/T cache before the gnomonic correction
([[cs-general-eos-stale-cache]]); tfloor_kelvin unrestartable ([[athenak-restart-tfloor-bug]]);
the C2P floor corrupting u.e on cs ([[cs-mhd-c2p-floor-corrupts-ue]]).

**Caveats, not bugs**
* Off-table states are SILENTLY continued "linear in the logs" beyond
  eos_logt_min = 1.5 (31.6 K) / eos_logt_max = 6.0 and the rho range. A cell that leaves the
  table gets an extrapolated p and T with no warning. The floors are meant to prevent it;
  worth an event counter if the top of the atmosphere ever cools below ~30 K.
* Ghost and active copies of the same cell can converge to roots differing by ~1e-13
  (different warm starts), so same-panel block-face conservation is 1e-13, not round-off,
  under the general EOS -- Fable's point. Real, tiny, and NOT a seed of the cs-vs-pert
  difference field (both arms share it). See [[cs-seam-ghost-eint]].
