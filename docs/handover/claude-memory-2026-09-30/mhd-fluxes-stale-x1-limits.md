---
name: mhd-fluxes-stale-x1-limits
description: "3429f59f: four well-balanced call sites in MHD::CalculateFluxes x2/x3 passed il,iu (the x1 sweep's limits, never reset) instead of is-1,ie+1. Latent, not live -- measured to change no answer."
metadata:
  node_type: memory
  type: project
---

Found by auditing for more instances of [[dhj-run-to-run-nondeterminism]]. In
`MHD::CalculateFluxes`, `il`/`iu` are set ONCE near the top for the x1 sweep
(`is, ie+1`; `is-1, ie+2` with FOFC) and **never reset for x2/x3**, where every other
limit is `is-1, ie+1`. Four sites still used them:

* `WbLocalPiecewiseLinearX2/X3` -- the worst, because it is a *producer*: it sat one line
  above `PiecewiseLinearX2(..., is-1, ie+1, b0_, ...)`, so the same j-interface was
  reconstructed over two different ranges, leaving `wl_jp1/wr` unwritten at `i = is-1`
  while the Riemann solver reads it.
* `AddWbPrimFaceX2/X3` -- a consumer.

**Not live, and I overclaimed it before measuring.** `<mhd>/wellbalance_dynamic` and
`wellbalance_static_reconst` both default false and are set in NO input in the tree.
Measured with `wellbalance_static_reconst=true` on the general-EOS dhj setup (production
grid, 600 cycles, 5 replicates, `wb.athinput` in the GPU scratch): pre-fix and post-fix
**bitwise identical**, and both reproducible run to run. Two reasons: `w0facewb` is zero
at the ghost `i = is-1` (pgens fill active cells only), and the predicted lane-remap race
never tripped.

Hydro is clean -- `hydro_fluxes.cpp` DOES reset `il,iu` for x2/x3 and uses them
consistently everywhere in those kernels.

**Audit result, worth not repeating:** no second instance of the 6e600f12 bug exists.
Every other consumer of an X1 reconstruction respects the `ql(n,i+1)`/`qr(n,i)`
convention; `reconstruct/*` repeated inner loops are if/else branches; viscosity.cpp has
no barriers but needs none (index-i to index-i throughout); shearing box, orbital
advection, current_density, hydro/mhd_update, resistivity_update are all correctly
barriered; dyn_grmhd barriers immediately after reconstruction. The plain-`par_for`
class (read a neighbour you also write) is clean in the evolution path too.

**Method for both sweeps** (scripts were throwaway): parse each `par_for`/`par_for_outer`
body by brace matching, then (a) flag consecutive `par_for_inner` with no `team_barrier`
between and an offset scratch index, (b) flag arrays that are BOTH assigned and read at
an offset index in one launch. Distinguish assignment targets from reads or the face-field
idiom `x1f(i)=...; x1f(i+1)=...` buries the signal.
