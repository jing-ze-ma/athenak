---
name: general-eos-mhd-polar-bug
description: "RESOLVED: eos=general vs eos=ideal differ under MHD+polar because the polar boundary forces HLLD->HLLE at the pole"
metadata: 
  node_type: memory
  type: project
  originSessionId: 85cd606f-82a3-4a02-a9d2-5e410a10595a
  modified: 2026-08-09T10:26:17.011Z
---

Investigated 2026-08-09 while verifying `deep_hot_jupiter_rt`. Under `eos = general`
evaluating a gamma law, MHD + `use_polar_boundary = true` diverged from `eos = ideal`
(~5e-3 in `bcc2`). **NOT A BUG — diagnosed and explained. Do not re-chase.**

**Cause:** `mhd_fluxes.cpp` substitutes **HLLE for HLLD at the pole** — `do_pole` at lines
216 (x1), 408 (x2), 604 (x3), true for `j == js/js+1/je/je+1` on a polar boundary. HLLE is
the one solver already known to differ between the ideal and general paths: the general
path uses a Davis/Einfeldt wave-speed bound because the **Roe average is only defined for
an ideal gas**. So the documented HLLE truncation-level difference surfaces exactly, and
only, in the pole cells. `hydro_fluxes.cpp` has NO `do_pole` — hydro never swaps solvers —
which is why hydro + polar is bitwise identical.

Evidence (1-cycle runs, difference localized per block by theta index):
- `rsolver=hlld` -> difference ONLY at the pole (blocks at logical j=0 differ at j=0,1,2;
  logical j=1 at j=5,6,7); the internal theta block boundary is exactly zero
- `rsolver=hlle` everywhere -> difference spread across ALL theta cells
- `rsolver=llf` -> NO pole difference at all
- hydro + polar -> bitwise; MHD + reflect BCs, same physics -> round-off

**Consequence to accept:** under a general EOS the pole cells use a different wave-speed
estimate than an ideal-gas run does. Both are valid discretizations; they cannot be made
bitwise equal, because the ideal HLLE's Roe average has no general-EOS analogue. This is
inherent, not fixable.

**Second finding, also RESOLVED and also benign — round-off amplified at the axis.** With
`rsolver=llf` (no pole solver swap) the general and ideal paths still differ slightly, and
the amplification comes from the polar axis, not from stretching or resistivity:

| grid (llf, 40 cycles) | B difference, relative to max\|B\| |
|---|---|
| reflect BC, no stretch, theta 0.35-2.79 | 0 (below float32 output) to 2.6e-14 |
| reflect BC + theta stretch, same range  | 0 to 1.2e-13 |
| POLAR BC + stretch, theta 0-pi          | 1.4e-7 to 5.2e-7 |

So the axis amplifies the seed by ~6 orders of magnitude, but the end result is still only
~5e-7 relative — at the resolution of the output itself. Resistivity is NOT involved
(removing it changes nothing) and stretching alone is NOT involved.

The seed is last-bit arithmetic: the general path evaluates pressure and the fast speed by
a different but mathematically equivalent route (`SoundSpeedFromP`/`FastSpeedFromP` from
reconstructed `wder`) than the ideal path does from `w0(IEN)`. A temporary diagnostic that
forced the general path to use the ideal's interface pressure
(`dl(IDPR)=(gamma-1)*wl(IEN)`) changed the answer NOT AT ALL, which rules out the
reconstruction non-commutation and points at the wave-speed evaluation.

**METHODOLOGY CAVEAT that applies to ALL the general-EOS verification in this project:**
AthenaK's `bin` output is written in SINGLE precision (`float *single_data`,
`src/outputs/binary.cpp:145`). "Bitwise identical" from comparing .bin files therefore
means agreement to within float32 rounding, ~6e-8 relative — strong, but not literally
bit-identical doubles. Differences below that floor cannot be seen this way.

See [[general-eos-project]]. Reproduce with
`/orion/u/jinma/ATHENAK/eostest/dhjmhd_{ideal,gen}.athinput` (needs `OMP_NUM_THREADS=1`;
the login node oversubscribes otherwise). Note `mhd/reconstruct=...` on the command line is
SILENTLY IGNORED on a stretched/spherical-polar grid — the code forces
`GridPiecewiseLinearX*` regardless of `recon_method`.
