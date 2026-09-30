---
name: nx1-ceiling-lds
description: "Radial resolution is capped at nx1 = 264 (general EOS) by the MHD flux kernels' 64 kB LDS budget, and x1 CANNOT be split across meshblocks because the RT is a whole-column solve. Production runs nx1=256, i.e. 8 cells of headroom. Escape hatch: scratch level 1 costs 1.1% and removes the cap entirely."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-22T00:15:00.000Z
---

## The cap

`mhd_fluxes.cpp` gives the x2 and x3 flux kernels NINE LDS scratch arrays of length
`ncells1` -- three staggered copies of `nvars` primitives, 3 magnetic components and `nder`
derived quantities -- against a hard **64 kB per workgroup**. That is a fixed number of
bytes of LDS per RADIAL CELL of the meshblock:

| EOS | rows | B per cell | measured ceiling |
|---|---|---|---|
| ideal (`nder = 0`) | (5+3)x3 = 24 | 192 | nx1 <= 326 (324 OK, 388 FAIL) |
| **general/table (`nder = NDERIVED = 2`)** | (5+3+2)x3 = 30 | 240 | **nx1 <= 264** (264 OK, 272 FAIL) |

Fitted: LDS = 192*ncells1 + 2144 for the ideal case. nx1 = 512 would want ~124 kB, hence
Kokkos's "could not find a valid team size" rather than a clean message.

**`inputs/mhd/deep_hot_jupiter_rt_eos.athinput` uses `eos = general`, and production runs
nx1 = 256. That is 8 cells, ~3 %, from a hard wall.** Note `correctness.athinput` in
bench/polar_ab uses `eos = ideal` (table only for x_e), so it has a HIGHER ceiling and does
not show the problem -- do not calibrate on it.

Three things push it over with no warning: raising nx1; raising `nghost` (ncells1 =
nx1 + 2*nghost, so 2 -> 4 costs 4 cells); adding passive scalars (`nvars = nmhd +
nscalars`, each scalar is 3 rows = 24 B/cell = ~7 cells of ceiling).

## Why the obvious workaround does NOT apply

`ncells1` is set by `meshblock/nx1`, not `mesh/nx1`, so normally one would split x1 across
meshblocks. **That is not available here: the RT is a whole-column solve.** It sweeps i
from is to ie+1 inside ONE meshblock, and the top boundary, the cumulative stellar optical
depth and the Planck bottom boundary at the cut all assume the block spans the atmosphere.
Split x1 and every block solves a truncated atmosphere with a fabricated top. So the LDS
cap is a cap on RADIAL RESOLUTION outright.

## The escape hatch: scratch level 1

`par_for_outer` takes a scratch level; level 1 is global memory with no 64 kB limit.
`mhd_fluxes.cpp` hardcodes `int scr_level = 0`, but **`dyn_grmhd` already does exactly this
properly**: `scratch_level = pin->GetOrAddInteger("mhd", "dyn_scratch", 0)` in
`dyn_grmhd.cpp:145`, used at `dyn_grmhd_fluxes.cpp:72`. The upstream authors hit this wall
and made it configurable, just not for non-relativistic MHD.

Measured with `scr_level = 1` (general EOS):

| | nx1=256 | nx1=272 | nx1=512 |
|---|---|---|---|
| level 0 (current) | OK | FAIL | FAIL |
| level 1 | OK | **OK** | **OK** |

Cost at nx1 = 256: CalculateFluxes 97.9 -> 103.4 ms (+5.6 %), **all GPU 614.9 -> 621.9,
just 1.1 %**. The cap disappears for ~1 % of runtime.

**NOT APPLIED, and the user decided 2026-08-22 not to go that far.** Keep it as a known
option, do not re-propose it unsolicited. Revisit only if a run actually needs nx1 > 264,
which is when the cap stops being theoretical.

If it is ever done, the change mirrors dyn_grmhd exactly: read `mhd/scratch_level`
(default 0, so no behaviour change) and pass it through to `par_for_outer`. Four lines
across `mhd.hpp`, `mhd.cpp` and `mhd_fluxes.cpp`. Caveat: level 1 was verified to RUN, not
verified bitwise identical -- check that before trusting any science from it.

**UPDATE 2026-09-12: the escape hatch IS APPLIED for HYDRO (commit 955c39df):** `<hydro>/scratch_level`
(default 0 = LDS, 1 = global) is passed to the three hydro flux kernels. First user: the red-giant
FOFC run (nx1 = 480, nghost = 3, general EOS) -- level 0 aborted at cycle 0 on MI300A with
"could not find a valid team size", level 1 runs (smoke 11625322). MHD still hardcodes level 0;
mirror the same four lines if an MHD run ever needs nx1 > 264. Bitwise level-0-vs-1 NOT checked.
