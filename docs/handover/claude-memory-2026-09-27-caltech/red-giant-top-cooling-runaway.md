---
name: red-giant-top-cooling-runaway
description: "The red giant runs die from the OPTICALLY THIN TOP: the top cell cools 3364 -> ~1100 K within 2e4 s, the layers evacuate, and dt collapses in HYDRO at t = 3.7-4.4e5"
metadata:
  type: project
---

All three jobs left running on 2026-09-08 (193861 taufix100, 193862 taufix300, 193869
combo) died the SAME way: a **hydro** dt collapse, then a segfault a few cycles later.
`taufix100` t = 3.80e5, `combo` t = 3.66e5, `taufix300` t = 4.44e5 -- a deeper tau blend
DELAYS it and does not fix it. The dt-collapse line names the module but not the cell,
which is why this took a session; a hydro version of the conduction cell report is the
obvious next instrument.

## Where it comes from (measured, not inferred)

Mass and internal energy are conserved to 3e-5 over the whole run: **the deep envelope is
innocent** and the `face budget: inner = -7600 L` line is a gross flux, not a leak.
Everything happens in the top ~10 radial cells:

- t = 0 the top is uniform at 3311 K (offline EOS on the dump; the RT column says 3363.6).
- t = 5e4 the top cell is at ~1000 K and the 8 below it at 2200-2700 K.
- the top cell's tau_R goes 1.4e-5 -> 1.11e-1 (all of it in that ONE cell), density
  swings 27x, |v| grows 1.15 -> 25 km/s (Mach ~4), and the emergent flux falls
  1.074e10 -> 9.74e9 erg/cm2/s. **The emergent-flux decay (the old 0.877 -> 0.556) is
  this same top-layer collapse, not a deep-interior problem.**
- KE decays 6.4e5 -> 7.2e4 and then creeps back to 9.6e4: convection does NOT develop,
  and what KE there is at the end is the top runaway, not plumes.

## The mechanism, from `problem/rt_apply_debug` (commit 351a6f44)

At t = 0, top cell: emission Em = 5.325e-2 erg/cm3/s, net divF = -7.0e-4, i.e. it absorbs
**1.3 % less than it emits**. With a grey opacity that would settle 11 K lower and stop.
It does not: over the next calls, as T falls 3363.6 -> 3348.3 -> 3339.2, the cooling
divF *grows* (-7.02e-4 -> -7.26e-4 -> -7.27e-4). Cooling that increases as the cell cools
is a positive feedback and is what runs the top away.

Suspected cause: in the correlated-k IR sweep both the emission and the absorption use the
band opacities at the LOCAL temperature, but the incident field is at the (hotter)
radiation temperature. As the cell cools, its Planck weighting shifts red while the
incident flux stays blue, so absorption falls faster than emission and the imbalance
grows. NOT yet proved -- the clean test is the equilibrium curve divF(T) at fixed rho and
fixed incident flux, which needs one more instrument.

## What is ruled out
- the `rt_de_max = 0.5` limiter: no clip warning, and de/e is ~1e-3 per step
- the tau blend: w = 0 at every cell involved, the two-stream owns them outright
- mass or energy loss through either radial boundary
- the initial column: `column.txt` is smooth and correct, T -> 3364 K at the top

## Also seen
A finished run (`nlim` or `tlim` reached) **aborts in Kokkos finalize** -- signal 6, six
ranks, after the last output. Separate from the mid-run segfault; harmless but noisy.

Related: [[red-giant-dt-collapse-solved]] (the earlier, DIFFERENT collapse),
[[red-giant-session-2026-09-08]], [[correlated-k-design]].
