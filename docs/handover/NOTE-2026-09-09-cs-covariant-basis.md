# Cubed sphere: covariant / contravariant basis bugs — audit, fixes, what remains

Written 2026-09-09 on orion. Companion to `NOTE-2026-09-09-wb-restart-bug.md`.

## The convention (verified, file:line as of polar-average-perf 32273ebb)

On the cubed sphere the two angular directions are NOT orthogonal; the angle between the
unit tangents has `cos = pcoord->cos_cell(m,k,j)` (up to |c| = 0.4755 at a cube vertex,
~0 at a panel centre, +0.03 mid-seam). The code uses two bases:

- `w0(IVY,IVZ)` are CONTRAVARIANT velocities v^xi, v^eta (coordinates.cpp 774-833).
- `u0(IM2,IM3)` are COVARIANT momenta m_i = rho g_ij v^j (coordinates.cpp 1026-1040,
  `GnomonicEquiangleLowerMom`): m2 = d (v2 + c v3), m3 = d (v3 + c v2).
- `bcc0` and the Riemann-solver states are ORTHONORMAL (gnomonic_kernels.hpp 47-100).

Correct kinetic energy:
- from momenta: `0.5 [m1^2 + (m2^2 + m3^2 - 2 c m2 m3)/(1 - c^2)] / d`
- from velocities: `0.5 d (v1^2 + v2^2 + v3^2 + 2 c v2 v3)`

Off the cubed sphere c = 0 and both reduce to the plain sum, so every fix below is
bit-identical on spherical-polar and Cartesian grids.

The failure mode: any code that forms `d*v2`, `d*v3` or `0.5 d sum v^2` from `w0`, or
`0.5 sum m^2 / d` from `u0`, is wrong by O(c) at the panel corners and by nothing
elsewhere. That is why these bugs show up as "the 8 cube-vertex columns and nothing
else" (red giant, 2026-09-09: 24 vertex cells 5-8 % cold in the sponge layer of prod11;
lethal in the open-top runs once transonic).

## Core code: CLEAN (audited 2026-09-09, read-only, every site checked)

C2P (`GnomonicEquiangleRaiseVel` over ghosts, `defer_cons_floors` on cs), all Riemann
solvers, `SrcTermsGnomonicEquiangleImpl`, history.cpp, newdt, conduction.cpp,
prolong_prims (refused on cs), fofc/viscosity/most srcterms (refused on cs), all bvals
seam transforms (`TransformMomentum`), two_stream_rt/atm_column (no velocity use),
deep_hot_jupiter_rt.cpp (Coriolis, ghosts; its drag/sponges scale the momentum VECTOR,
which is basis-independent — they thermalise the removed KE into e_int, a physics choice
worth a look but not a basis bug).

**The deep hot Jupiter runs are not affected by anything in this note.**

## Fixed on orion (working tree; commit + push pending — see "status" below)

`src/pgen/red_giant.cpp`:
1. `rg_sponge`: KE from the orthogonal momentum sum -> `CsKinetic(d,m1,m2,m3,c)`.
   The wrong KE, subtracted and re-added times fac^2, drained the TRUE internal energy
   by (KE_true - KE_orth)(1 - fac^2) every step. Binary 887c241e.
2. `rg_ic`, `fill`, `fill_open`, `fill_open_out` (radial ghost fills): `u0(IM2)=d*v2`
   etc. -> the LowerMom formulas with `c = cs ? cos_cell(m,k,j) : 0`. Note cos_cell is
   allocated over the full padded angular range, so it is valid in radial ghosts.
3. `rg_co50`, `rg_wallflux` inner-shell energy accounting -> `CsKinetic`.
   Tested (T14, restart of the failing run): active, changes the vertex columns in the
   3rd-4th digit, does NOT cure the vertex collapse (that is a separate mechanism still
   being bisected). Keep the fix regardless: it is a correctness bug.

## NOT fixed — for viper (or whoever gets there first)

All latent: none is in a configuration currently being run.

| site | what is wrong | fix |
| --- | --- | --- |
| `src/srcterms/ismcooling.cpp` ~225-232 (drag) | mixes contravariant `w0` velocity into the covariant `u0` momentum | use `u0(IM2,IM3)` directly (drag scales the momentum vector); the guard comment at mesh.cpp ~444-447 is wrong too |
| `src/pgen/solid_body_rot.cpp` PrimToCons on cs | orthonormal `PrimToCons` | `pcoord->GnomonicEquiangleLowerMom(...)` as cs_test.cpp ~996 does |
| `src/outputs/derived_variables.cpp` v_moments / sgs | orthogonal velocity sums | add the `2 c v2 v3` cross term (output only; does not touch the evolution) |
| `src/pgen/cubed_sphere.cpp` | legacy, dead | delete or ignore |
| `src/pgen/red_giant.cpp` `rg_co5B` ~2235-2286, `rg_co5D` ~2312-2340 | `u0(IM2)/rho` treated as contravariant v2, orthonormal KE, `r*v2` written back | needs the full raise `v2 = (m2 - c m3)/(d(1-c^2))` then LowerMom; partly self-cancelling, small residual |
| `src/pgen/red_giant.cpp` `rg_mlt_shell` ~1801 | orthogonal sum in a shell-mean diagnostic | add the cross term (mlt is off in production) |

### How complete is this list? (checked 2026-09-09 21:00)

A grep over all of `src/` for the tell-tale forms (`d*v2` into `u0(IM2)`, `SQR(IM2)+SQR(IM3)`
from `u0`, `v2*v2+v3*v3` / `SQR(v2)+SQR(v3)` from `w0`) hits, besides the sites above:
`pgen/{solar_convection,cooling_convection,turb,wb_column,field_loop,slotted_cyl,
disk-magnetosphere,sp_test}.cpp`, `pgen/tests/{linear_wave,diffusion,orszag_tang}.cpp`,
`srcterms/turb_driver.cpp`, `rsolvers/roe_hyd.hpp`, `coordinates.cpp:2137,2218`,
`cs_test.cpp:2112,5915`. All ruled out: the pgens and turb_driver never run on the cubed
sphere (no `use_cubed_sphere` path), roe_hyd gets the orthonormal state, coordinates.cpp
2137/2218 are the SPHERICAL-POLAR source (`spsrc`, orthogonal grid), and the two cs_test
lines are error norms of a static (v = 0) test. `deep_hot_jupiter_rt.cpp` has no hit.
So the table above is the complete list of cubed-sphere-relevant sites.

Proposed shared helpers (not written yet), to stop this recurring: in
`src/coordinates/gnomonic_kernels.hpp` device functions
`GnomonicKineticFromMom(d,m1,m2,m3,c)`, `GnomonicKineticFromVel(d,v1,v2,v3,c)`,
`GnomonicLowerXi/Eta(d,v2,v3,c)`; then replace red_giant's private `CsKinetic`.

Rule for any new cs pgen: never write `u0(IM2,IM3)` from a velocity or form |v|^2 from
either array with the plain sum. Go through `GnomonicEquiangleLowerMom` / the helpers.

## Status of the orion working tree (2026-09-09 evening)

Pushed: 5c0b98e4 (well-balanced restart cache; affects EVERY WB restart incl. dhj —
viper must pull), 32273ebb (its note). Uncommitted on orion: the red_giant fixes above,
the shared two_stream_rt.hpp changes (direct source, Newton solve — ON by default, NOT
yet A/B-tested on the hot Jupiter, so do not assume they are neutral there), conduction
and hydro options (default = old behaviour), the every-cycle NaN checker in driver.cpp,
two default-off debug flags (`mesh/cs_corner_poison`). They will be committed in pieces.
