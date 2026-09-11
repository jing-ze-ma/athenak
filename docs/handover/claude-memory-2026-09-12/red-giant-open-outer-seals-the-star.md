---
name: red-giant-open-outer-seals-the-star
description: "outer_bc = open stopped the drain but SEALED the star; the column above the domain now radiates (bd2b974d, problem/rt_top_re) and the seal is gone. Its closing claim -- that the upstream defect is the handover carrying 0.34 L -- is SUPERSEDED: see red-giant-flux-deficit-is-spinup"
metadata:
  type: project
---

Measured 2026-09-08 late on the 6x8x8x320 testbed in `/orion/ptmp/jinma/Athenak/red_giant`.
Tools kept there: `topdrain.py` (per-dump top-cell rho/T/v and the mass in the top 10
cells), `topprof.py` (radial profiles), and the offline EOS grid `eos_grid_rg.bin` built
with `tools/solar_convection/eos_dump.cpp` **patched to set `include_metal_ion = true`** --
the tool does not expose it and the red giant runs with it on. Worth adding as an argv.

## 1. The open outer boundary stops the drain and seals the star

`greyopen2` (open, b865cc83 + 6d999933) vs `greylong` (wall), t = 3.75e5, horizontal mean:
top-cell rho 3.3e-10 (22x rho0) against 4.8e-11; rho_min 1.3e-10 against 2.5e-12 (no
evacuated cells at all); T 6863 K against 2633 K; max |v| 1.9 km/s against 22 km/s; mass in
the top 10 cells 13x and climbing. Below i ~ 260 the two runs agree to 4 digits -- none of
this involves the deep envelope. By t = 5e5 the outer 30 cells are an ISOTHERMAL 6650 K
blanket and the emergent net LW flux is 4.0e7 erg/cm2/s against the 1.31e10 the star must
radiate at that radius: **0.3 %**. The drain was traded for a seal.

## 2. Why, and the fix (COMMITTED bd2b974d)

`two_stream_rt.hpp` fills the incident down-stream at the top face from the unresolved
column of mass p/g above the domain, radiating at the GHOST cell's temperature. The open
ghost IS the top active cell, so the column returns exactly what that cell emitted, the
net flux at the face goes to zero as the column thickens, and the loop closes.

`problem/rt_top_re` (grey path only, module default OFF, red_giant defaults it ON) makes
that column a slab in radiative equilibrium: it re-emits half of what it absorbs upward to
space, so its source is I_up/2 and the face always keeps half its flux. I_up at the top
face does not depend on I_down at all, so it is obtained EXACTLY by a scalar probe sweep
before the down-sweep -- no lag, no stored state, one extra pass.

Gates: `rt_top_re = false` reproduces the old path to the last digit (4.5238466416e+09).
The column above is **not thin even at t = 0** -- tau ~ 0.23 backed out of the flux change,
so it was already returning a third of B(T_ghost) on the first call -- and turning it on
raises the t = 0 emergent flux 4.524e9 -> 5.891e9.

> **SUPERSEDED, 2026-09-08 night.** Section 3 below concludes that the upstream defect
> is the handover carrying only 0.34 L. That is WRONG. The 0.34 L is a convective
> spin-up transient: see [[red-giant-flux-deficit-is-spinup]], which measures it against
> the layer-by-layer turnover time. Everything in sections 1 and 2 still stands, and
> `rt_top_re` is on in production.

## 3. What it did NOT fix -- read before calling this closed

`greyre` (job 194206, open + rt_top_re, to t = 6e5): the 6650 K blanket and the 22x
pile-up are GONE and dt is back to 58.9 s, the same as the wall run. But the atmosphere has
reverted to the COOLING failure: top-cell T oscillating 2100-2800 K, the top 10 cells'
mass swinging 0.3-2.4x, rho_min down to 4e-13, |v| 6-14 km/s, and the emergent flux at
t = 3e5 is **2.97e9, lower than the wall run's 4.44e9**. Two failure modes have been
traded, and the fix is right but not sufficient.

The upstream defect is the HANDOVER, and it is the thing to work on next: even the wall run
only ever carried **0.34 L** through the thin atmosphere (F flat at 4.44e9 vs 1.31e10
required). `icut` sits at r = 3.38e12, INSIDE the convection zone (RCB 3.426e12), where the
radiative flux genuinely is a fraction of L -- and the resolved convection that should carry
the rest never develops at 8x8 angular cells. With the RE blanket the skin temperature is
now self-consistent with whatever flux the solver actually delivers, which is why a
one-third flux gives a ~2500 K top instead of the 3364 K the initial column assumes.

## State at handover
- Job 194206 `greyre` (open + rt_top_re) RAN TO tlim = 6e5 with dt 58.9 and no collapse,
  then hit the known harmless abort in Kokkos finalize. 194207 `greyre_off` (the gate,
  the old boundary on the same binary) was still going, slower, as the old path does.
- `grey_prod` job 194160 (wall, full 6x32x32x320, tlim 3e6) was past t = 1.85e6 on the OLD
  binary, dt flat 30.64 -- unaffected by any of this and still the best long run.
- Build for the new binary: `build_rg_topre` (`module load gcc/13 openmpi/4.1`, and cmake
  needs `-D CMAKE_CXX_COMPILER=$(which mpicxx)` or it picks up the system gcc 7.5).

Related: [[red-giant-grey-two-stream]], [[red-giant-top-cooling-runaway]].
