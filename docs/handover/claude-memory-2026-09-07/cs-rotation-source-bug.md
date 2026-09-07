---
name: cs-rotation-source-bug
description: The cubed sphere fell through into the CARTESIAN BOX's equatorial beta-plane rotation source while spherical polar got full Coriolis+centrifugal -- so "same physics, cs dies and sp does not" was FALSE. Fixed f75ad783, <problem>/cs_full_rotation
metadata:
  type: project
---

**Found 2026-09-04, fixed in f75ad783.** `SourceFunc` in `deep_hot_jupiter_rt.cpp`: the
COORDINATE block branches three ways (spherical polar / cubed sphere / Cartesian box); the
FORCE block ten lines below branched only TWO. The cubed sphere therefore ran the Cartesian
box's equatorial beta-plane:

    Real omega1 = omega*lam;                       // f = 2 omega * lam, lam in RADIANS
    u0(IM2) += -2 rho omega1 (-v3) bdt;
    u0(IM3) += -2 rho omega1 v2 bdt;

while spherical polar got full Coriolis + centrifugal + the centrifugal work term.

**The beta-plane is not itself wrong** -- `f = 2 omega lam` IS the standard equatorial form
(f = beta*y, beta = 2 omega/R, y = R lam; Fromang+2016). It is wrong HERE three times over:
it is an EQUATORIAL approximation on a GLOBAL grid (11 % too strong at 45 deg, 21 % at 60,
57 % at the pole against 2 omega sin lam); it drops the centrifugal term entirely; and its
two terms assume an ORTHONORMAL (x,y) pair while IM2/IM3 on this grid are COVARIANT
components on the NON-ORTHOGONAL gnomonic tangent basis.

## Why this matters more than the physics error

**The premise of the entire cubed-sphere MHD investigation was false.** "The same problem
runs clean on spherical polar and dies on the cubed sphere" -- the two grids were never
solving the same problem. They differed in the ROTATION SOURCE, before any numerics. Every
cs-vs-sp production comparison in [[cs-mhd-dhj-blowup]], [[cs-dhj-production-retry]],
[[sp-hydro-vs-mhd-comparison]] and [[cs-mhd-instability-characterized]] is confounded by it.
A stronger high-latitude Coriolis and no centrifugal force change the jet, the shear, and
therefore WHERE AND WHEN beta crosses 1 -- which is exactly the shell where
[[cs-ulp-amplification]] measured the round-off growth.

Note this does NOT touch the idealised reproducers (cs_test iprob=8/13 run with omega = 0),
so [[cs-mhd-minimal-reproducer]] and [[cs-mhd-low-beta-divergent]] stand as measured.

## The fix, and how it was gated BEFORE running anything

Done in CARTESIAN and projected, rather than writing the curvilinear Coriolis by hand in a
non-orthonormal basis: `a = -2 Omega x v + Omega^2 (x,y,0)` with `Omega = omega zhat`, V
rebuilt as `v1 rhat + v2 e_xi + v3 e_eta` (w0's angular velocities are CONTRAVARIANT on that
pair -- that is what `GnomonicEquiangleRaiseVelMHD` builds, `v2 = (m2 - c m3)/(d det)`), and
the momentum source taken as `rho (a . e_i)` because `m_i = V.e_i` is COVARIANT.

Algebraic gate, `bench/../scratchpad/rotgate.py`, over all six panels and the full (xi,eta)
range INCLUDING the panel corners:

    cs branch vs the sp formula at the same physical point   3.0e-15  (vector AND
                                                                       each covariant slot)
    covariant -> contravariant -> Cartesian velocity round trip  4.9e-16

`<problem>/cs_full_rotation`, default TRUE. Kept switchable because this changes every
cubed-sphere answer and the old ones must stay reproducible. The Cartesian box is untouched.

## The bug CLASS, which has bitten before and will again

A three-way coordinate branch paired with a two-way physics branch. Audit of every
`use_spherical_polar` site in the pgen found:

* **live**: this one;
* **latent (dead code)**: `double_gray_two_stream_RT` and `double_gray_two_stream_RT_source`
  declare `Real z, r;` and assign `r` ONLY in the sp branch, then use it in `mucr`, `delta`,
  `drcor` -- an uninitialised read, i.e. UB, on any non-sp grid. Neither function is called
  today (checked: no call sites), so it changes nothing, but it would bite immediately if
  that path were re-enabled;
* **precedent, already fixed, comment still in the code**: `Real x1v = ...` inside the else
  branch SHADOWED the outer declarations, "leaving the outer ones uninitialised for every
  non-spherical-polar grid... the whole domain going NaN on the first source-term
  application";
* **checked and benign**: the restart `bcc0` fill takes the plain-average else branch on cs,
  but ConsToPrim -> GnomonicEquiangleRaiseVelMHD rebuilds bcc0 over the whole array
  including ghosts before first use.

**The telling contrast**: the stellar tide IS explicitly refused on cs with a fatal error,
commented "the gnomonic tangent pair is neither that basis nor orthonormal, so the
expressions would be silently wrong rather than merely untested". Someone reasoned about
exactly this hazard for the tide -- and the Coriolis term ten lines below had the same basis
problem, no guard, no comment. The codebase holds both patterns side by side; what separated
them was only whether anyone looked. The audit covered the pgen only, NOT the rest of `src/`.

Found by asking a Fable 5.1 agent for independent ideas on the instability; it read the code
and flagged the branch. Two of its other leads are unverified and worth pursuing: the seam
ghost cells carrying a total energy inconsistent with the separately-resampled ghost B (a
floor trigger IN THE HALO at low beta), and the tabulated-EOS warm start (`logtol = 1e-13`,
`eos_table.hpp:120`) making every block face a 1e-13 noise source, which may inflate the
block-face excess in [[cs-ulp-amplification]] -- though ctl and pert share that noise, so it
is not a direct seed of the difference field.


## FULL src/ SWEEP for the same bug class (2026-09-04)

Two passes: (1) every `use_spherical_polar` site in `src/`, checked for a cubed-sphere
companion; (2) a scripted search for the OTHER half of the pattern -- `Real a, b;` with no
initialiser, assigned only inside a grid branch, then read outside it.

**Core solver and infrastructure: CLEAN. No second live instance found.** Every
grid-dependent site in `bvals/` (flux_correct_cc, flux_correct_fc, prolongation),
`mesh/mesh_refinement`, `hydro/` (update, tasks, fluxes, newdt), `mhd/` (ct, update, tasks,
fluxes, newdt), `diffusion/` (resistivity, resistivity_ct, resistivity_update,
current_density, plus the dedicated resistivity_gnomonic.cpp) and `outputs/` (history,
derived_variables) either tests `use_cubed_sphere || use_spherical_polar` together or has an
explicit cs branch. That is a real negative result: **the rotation source was the only live
instance in compiled, reachable code.**

**Two EOS sites, second order, one half already known.** `ideal_mhd.cpp:82` and
`general_mhd.cpp:112` build bcc in ConsToPrim with sp's POSITION-WEIGHTED face
interpolation, and cs falls into the plain 0.5 average. The magnetic-energy half of this is
[[cs-mhd-c2p-floor-corrupts-ue]]. The half NOT recorded before: on a STRETCHED radial grid
the cell centre is not the midpoint of its two faces (x1v = S(midpoint), not
(S(left)+S(right))/2), so the two forms genuinely differ. **FIXED in ff99c522.**

MEASURED FIRST: on the production grid (128 radial cells, 4.07x spread in cell width) the
weight departs from 0.5 by at most 4.4e-3, mean 1.9e-3. It is a SECOND-ORDER ACCURACY fix
and will NOT move a dynamical result -- do not expect it to touch the instability. The
reason to make it is CONSISTENCY: ConsToPrim (ideal_mhd, general_mhd) and
`GnomonicEquiangleRaiseVelMHD` both rebuild bcc, and when they disagree the difference is an
energy mismatch the C2P floors bake into u.e -- which is exactly
[[cs-mhd-c2p-floor-corrupts-ue]]. All three sites now call one shared
`CellCenteredRadialFld` in `coordinates/cell_locations.hpp` so they cannot drift apart.

Applied ONLY when cubed sphere AND a stretched radial grid, because the weighted form equals
the plain average only up to round-off: applying it everywhere would perturb every Cartesian
and every unstretched cs answer in the last bits for nothing. x2/x3 keep the plain average
(the panel coordinates are uniform; `use_grid_stretch_theta` is refused).

Gated both ways: unstretched cs at 20 cycles is BIT-IDENTICAL to the pre-change binary; the
same problem with the production stretch on DOES change (2-mom by 3.8e-5, the radial
channel). A fix that changed nothing on the stretched arm would have been the failure.

**Dead code carrying the bug**: `mhd_corner_e_uct.cpp` is NOT in `src/CMakeLists.txt`, so it
is never compiled at all; `double_gray_two_stream_RT` / `_RT_source` are compiled but have
no call sites.

**Latent, saved by a guard 1500 lines away**: `SourceFunc`'s `Real lam, phi, z, theta, r;`
leaves `theta` and `r` unset on a CARTESIAN box. `theta` is only read inside the sp branch,
so it is fine. `r` is passed to `GravAccAt` unconditionally -- strictly an indeterminate
value passed by value -- but that function ignores r unless `grav_point_mass`, which is
REFUSED on any grid that is not sp or cs. Safe today only because of that guard.

**Sibling pgens carry the same idiom**: `deep_hot_jupiter_rt_old.cpp`,
`solar_convection.cpp`, `cooling_convection.cpp`, `solar_convection_old.cpp` all declare
`Real lam, phi, z, theta, r;` the same way; in solar_convection and cooling_convection the
beta-plane block is COMMENTED OUT, so there is no rotation source to get wrong.
`deep_hot_jupiter.cpp` (the older non-RT hot Jupiter) has the LIVE beta-plane and mentions
`use_cubed_sphere` ZERO times -- it predates cs support, so running it on a cubed-sphere
mesh would silently take the Cartesian path everywhere, not just for rotation. Only one pgen
is compiled at a time, so this is a trap rather than a bug.

## THE A/B: the beta-plane control DIES, the fixed arm LIVES

Run 2026-09-04, `bench/cs_rotfix` vs `bench/cs_rotctl`. Same binary (md5-identical), both
FROM SCRATCH at t = 0 (not restarts -- `start.rst` was evolved under the beta-plane, so
restarting it would confound the force change with the transient of switching forces),
nx1=128, nx2=nx3=64, everything else identical. The ONLY difference is
`<problem>/cs_full_rotation`:

    cs_rotctl  (beta-plane, false)   DIED at rot 0.2814   dt -> 0, hst all NaN
    cs_rotfix  (full Coriolis, true) ALIVE past rot 0.4577  dt 2.38, zero NaN,
                                     mass 3.45e26, every energy finite

**0.2814 lands exactly on the recorded band** (0.281 / 0.309 / 0.324 from
[[cs-dhj-production-retry]] and [[cs-ulp-amplification]]), so the control reproduces the
historical death with the historical value while the corrected arm walks past it. Every
earlier comparison differed from those deaths in three ways at once -- restart vs
from-scratch, older binary, ULP-perturbed vs clean -- and this A/B removes all three.

Still ONE trajectory per arm ([[validate-the-instrument]]), so the honest claim is: on a
paired from-scratch A/B the rotation source is the discriminator. Not yet: "the instability
is explained".

## The fallback is LOAD-BEARING: cs cannot run without it

`bench/cs_rot32_nofb` -- nx2=nx3=32 (2.812 deg, the EXACT angular match to sp_dhj_ctl's
180/64 and 360/128), `cs_full_rotation=true`, `cs_lowbeta_fallback=0.0` -- was meant to
remove the last cs/sp asymmetry, since the fallback is a cubed-sphere-only switch. It
**went NaN at rot 0.0898**, far earlier than any arm with the fallback on.

So the matched-resolution cs-vs-sp comparison CANNOT be run with the fallback off. That is
itself a result: **at matched angular resolution and with correct rotation, the cubed sphere
still needs dissipation that spherical polar does not.** The rotation bug was a real confound
and its removal does not remove the grid difference.

**THE NaN-GRIND SIGNATURE, now seen twice and worth memorising**: dt freezes at EXACTLY
1.04466e+01 and stays there for 20,000 cycles while `time` keeps advancing, and
`eos_efloor` overflows to EXACTLY -1135976192. Both values are identical to the ones in
[[cs-dhj-long-run]]. A naive "min dt over the run" check calls this a collapse-and-recovery;
a naive "last dt" check calls it healthy. Neither is true -- it is a dead run grinding. ALWAYS
grep the .hst for nan.

## AUDIT of deep_hot_jupiter_rt.cpp for cs-vs-sp deviations (2026-09-04, after the fix)

Method: enumerate EVERY site that reads or writes a tangential component (IVY/IVZ/IM2/IM3)
and every use of x2v/x3v as an angle, and classify each for "spherical-polar component
formula applied unbranched on the gnomonic basis" -- the rotation bug's class.

* **ICs** (UserProblem, HydrostaticEquilibrium): set IM2 = IM3 = 0. Zero in any basis. OK.
* **Sponge layers** (SourceFunc): isotropic scalar drag `m_i -> (1-f) m_i`, basis-free. The
  bottom sponge damps only IM2/IM3, which are the two TANGENTIAL components on BOTH grids. OK.
* **Every angle use of x2v/x3v** (`x3v - M_PI`, `sin(x2v)`) sits inside an sp-only branch,
  verified line by line. OK.
* **RT beam cosine** mu0 = sin(theta)cos(phi): every site takes theta/phi from the three-way
  grid branch (CSCellAngles on cs), including the precomputed cf(m,k,j,3). OK.
* **Gravity / potential**: GravAccAt(r), GravPotAt(r, z) with r = x1v, z = x1v - ap on both
  grids. OK.
* **x2/x3 user BCs**: all commented out. Dead.
* **OUTER-x1 USER BC: THE ONE LIVE INSTANCE, fixed in 8a7eae1d.** Ghost total energy was
  `e0 + 0.5(m1^2+m2^2+m3^2)/rho` -- exact on sp, but on cs it drops the `-2c m2 m3` cross
  term of `0.5 m_i G^{ij} m_j /rho`, up to 50 % of tangential KE at a panel corner
  (|c| = 1/2). Now metric-correct on cs, term-for-term identical elsewhere. Gated against
  RaiseVelMHD's contraction to 6e-16. Confined to the outer radial ghosts above the top
  sponge, so EXPECTED SMALL in production -- a correctness fix, not a candidate cure.
* The inner-x1 BC touches only the magnetic term, using bcc0, which IS orthonormal on cs. OK.

Not audited here: `src/` outside the pgen was swept earlier (clean); the RT column solve's
geometry (dx1, area, correct_spherical) uses pcoord quantities that are grid-aware.

## ENDPOINTS of the corrected-rotation cubed-sphere arms (chains exhausted 2026-09-04)

    cs_rot32   nx2=nx3=32, matched to sp_dhj_ctl at 2.812 deg, fallback 0.5, from scratch:
               rot 2.682, dt 11.64, 115993 cycles, ZERO NaN.  Through the dt trough (its
               minimum 3.35 at rot 0.64, vs sp's 1.20 at 1.10) and recovered to sp's
               plateau (~12).  Radial ME flat ~5.4-5.8e31 through the window where sp
               jumped 56x.
    cs_rotfix  nx2=nx3=64, from scratch: rot 0.7984, dt 1.71, 87480 cycles, ZERO NaN --
               past all three recorded deaths (0.281/0.309/0.324) and past sp's jump
               window with radial ME flat at 1.10e32.  Not extended further: the nx=32 arm
               already answers "usable through the trough", and the queue went to the sp
               polar test instead.

So with the rotation source corrected the cubed sphere runs the production problem for
2.7 rotations at matched resolution with no sign of the instability, and the one thing it
does NOT reproduce is spherical polar's own polar-axis field blow-up
([[sp-polar-field-blowup]]). One trajectory per arm; the beta-plane control died at 0.2814.

**SWITCH REMOVED (ae831665, 2026-09-05, user decision):** `<problem>/cs_full_rotation` no
longer exists; the cubed sphere ALWAYS applies the full Coriolis + centrifugal source. The
beta-plane answers (cs_prod_hyd 158 rot, cs_dhj_hyd 283 rot, cs_both_chain, the whole
low-beta arm matrix) are NOT reproducible on binaries after ae831665, by design.
