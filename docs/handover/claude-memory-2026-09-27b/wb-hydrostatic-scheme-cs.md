---
name: wb-hydrostatic-scheme-cs
description: problem/rot_potential (commit after 8c384e46): barotropic total-potential IC+ghosts+arrays, rotating column at rest deep v_r 4e-2 -> 3.4e-3, v_theta 2e-2 -> 4e-3 (sp), USE WITH wellbalance_dynamic. FULL-PHYSICS + GPU VERIFIED 09-07 (RT on: helps; rotation on: 'worse' deep v_r is the CONVERGED answer, plain scheme under-resolves the centrifugal adjustment at n128; HIP build OK). VERIFIED hydro/MHD/resistive on cs, sp AND Cartesian x1/x3 (pgen wb_column, ea0cd573); polytropic = T LINEAR in Phi (exact ideal polytrope). COST SOLVED (cbf38d20 + next commit): polytropic at 1.03x of a plain run with wb_cache_every=10 (1.34x every stage); ThermoAt one-eval accessor, 5-eval polytropic fast path, isentropic_dt/adaptive_fast, BuildWBCache cache. NEW closure wb_option=polytropic COMMITTED 140dbf9e (local d ln T/d Phi, no stored background). The Kappeli-Mishra deviation well-balanced scheme (<hydro|mhd>/wellbalance_dynamic + wb_x1 + wb_option) DOES support the cubed sphere and non-constant (point-mass) gravity -- measured 2026-09-06: a dhj column at rest gives radial KE identical to 4 digits on cs and sp, off and on; the scheme cuts the hydrostatic residual 33-58x but NOT to round-off; it is ~100x slower with the general EOS; never used in any production input
metadata:
  type: project
---

**Wiring (read 2026-09-06).** Potential-based: the pgen fills cell-centre and face
potentials (phicc0/phi0) from GravPotAt, which has the point-mass form
-g a (1 - a/r) under problem/grav_point_mass, at the STRETCHED faces (x1f_) and the
centroids (x1v_). The background walk (utils/wb_background.hpp, general-EOS capable:
isodensity/isothermal/isentropic/adaptive via wb_option) uses only potential
differences. The gravity momentum source under wellbalance_dynamic is the background's
own pressure difference, src = (A_r (p_r - p) + A_l (p - p_l))/V with area1/volume, and
the radial geometric source on BOTH grids uses z_ov_rE = (A_r - A_l)/V (coordinates.cpp
765 sp, 1628 cs) -- the same discrete curvature -- so the cancellation is consistent on cs
(default Impl source) and on the face-sum cs_wellbalanced_src path too (measured, below).
The x1 sweep on sp/cs (GridPiecewiseLinearX1) takes the WB flag. TWO flags are needed:
wellbalance_dynamic AND wb_x1 (else the reconstruction silently ignores it); wb_option
is required (GetString). Flags default off; NO input in the tree sets them; the old
static variant (wellbalance_static_reconst) was a no-op (ghost background zero).

**Test ($S/wbpair, $S/wbres; scripts wbgate.py, wbcmp.py; serial build_dhj_cpu).** dhj
column at rest: production hydro inputs, rt_ck=false, omega=0 (Rayleigh drag then = 0),
16x16 (cs) / 16x32 (sp) angular, nx1 128. hst 1-KE (radial KE) at t = 58/116/377 s:
    off: 7.018e30 / 2.658e31 / 2.187e32     on: 2.147e29 / 6.762e29 / 3.797e30
IDENTICAL on cs and sp to 4 digits (the column is 1D; the cs radial path == sp), and the
csmhd pair (cs_mhd_prod input, bbot 1e-6, cs_wellbalanced_src face-sum ON) matches to
5e-4. So cs is supported exactly as sp. The scheme reduces the residual 33x -> 58x but
the ON run still grows as t^2: NOT round-off. Deep band (i<16) rms v_r at t=1000 s (sp):
off 2.5e-2 km/s, on 1.7e-3 (15x). Bands 16-112 reach 0.1-0.7 km/s in BOTH arms by
2000 s and grow: the RT-free column CONVECTS (H2-dissociation region, cf.
[[upper-atm-mottling]]); use the deep band / early times as the gate. The ON runs are
~100x slower (general-EOS root finds in the background walk): sp_on ~5-10 s/cycle on one
core at 65k cells. Mass drifts ~1e-3 in 2000 s on both arms = the open top, not the scheme.
**Resolution scaling (sp, nx1 64 vs 128, 20 cycles, $S/wbres), 1-KE at t~405 s:**
off 3.72e33 -> 2.46e32 (15.1x, KE ~ dx^3.9, i.e. v ~ dx^2); on 5.0e31 -> 4.25e30 (11.8x,
KE ~ dx^3.6, v ~ dx^1.8); ratio off/on 74 (n64), 58 (n128). So the ON remainder is a
converging O(dx^2) closure error (the isothermal/isentropic family does not contain the
radiative profile; the adaptive selector in hydro.hpp:1023 uses e/rho^gamma with the
IDEAL gamma as its entropy proxy even under the general EOS), NOT a missing cs term. The
scheme is a constant-factor ~60x gain, not exactness, on this profile.

## 2026-09-06 late: the polytropic closure (140dbf9e) and what the column test CAN'T show

**wb_option = polytropic** (hydro + mhd; ideal closed form in getWBerho, general-EOS
branch wb_opt==3 in WBAdvance with dlntdphi from WBBackgroundStencil): ln T linear in
Phi with the gradient from the two neighbours; contains isothermal (a=0) and isentropic.
Ideal-gas T is (gamma-1) e/rho in code units (my first cut used e/rho: 1 km/s garbage).
Results, sp column at rest, 20 cycles, bottom-band rms v_r [cm/s] at ~450-900 s:
    ideal  n64: off 3.9e3  adaptive 5.8e2  poly 5.1e2 | n128: off 6.1e2  adaptive 1.1e2  poly 1.2e2
    table n128: off 1.6e3  adaptive 1.3e2  poly 1.6e2
Cost (20 cycles, n128 serial): table off 5.3 s, adaptive 49 s, poly 41 s; ideal 1.75/2.27/2.34 s.
So: equal to adaptive, cheaper on the table, but NO order gain here, because the residual
with WB on (~1e-4..3e-4 g, growing linearly in t over the 500 s window = the column's slow
global readjustment) is set by (i) the IC's own inconsistency with the discrete
equilibrium (offline walk between neighbouring IC cells: 1.6e-5 face mismatch even in
the exactly isothermal top, ~2-6e-5 in the radiative zone, the SAME for all three
closures) and (ii) a kink at the IC's radiative-convective boundary (i~24-32, nabla
0.29 -> 0.08) that no stencil closure fits. The cell-average anchor correction
((dx/H)^2/24 offset) was implemented, made things 2x worse in the interior and blew up
the wall cells, and was REVERTED: the offset is nearly equal for two neighbours and
cancels at their shared face at leading order -- its numerical match to the residual
was a coincidence. GravPotAt stores +g a(1 - a/r) because grav_acc is NEGATIVE in the
pgen (Phi increases outward): an offline reconstruction must use that sign.
**To actually show the closure's benefit:** coarser cells (dx/H > 0.3) or, better, the
production RT run with WB on vs off (deep-band spurious v_r, mass drift), where the
profile is smooth and evolving. Scripts: $S/wbgate.py, wbcmp.py, offline walks in the
session transcript (wbideal/, wbres/ dirs).

## 2026-09-06 23:10: the cost fix (commit after 140dbf9e)

Where the 9x went: FOUR full stencil walks per cell per stage (rho channel, e channel,
pressure channel in the Der reconstruction, and the pgen source) each with ~a dozen
table evaluations plus temperature root finds. Fix 1: WBT() hands the walk ConsToPrim's
cached T (wtemp; valid in ghosts too, C2P runs over the full array) -- only 25 % gain,
root finds were NOT the main cost. Fix 2 (the one that matters): Hydro/MHD::BuildWBCache
walks the stencil ONCE per cell (WBBackgroundStencil returns all three channels) into
wbq0(m,15,k,j,i) at the start of CalculateFluxes for the x1 sweep's rows; the four X1
reconstruction functions and the pgen source read WBReadCache(). First attempt called
getWBq0 three times per cell and over ghost rows and was SLOWER (57 s) -- one call, sweep
rows only, gives 15.7 s. Timings (20 cycles, n128, serial): table off 5.3 / adaptive
48.9 -> 15.7 / polytropic 40.6 -> 14.6 / isothermal 9.6; ideal adaptive 2.27 -> 2.32.
Remaining overhead is the isentropic Newton (adaptive picks it deep) and the table
lookups themselves. Dumps: density/energy bitwise, v_r 2.6e-9 (the cached T vs a
re-solved one). Memory: 15 doubles per cell (94 MB for production). GPU untested.
TRAPS this session: (a) `make` failing silently while the timing runs reuse the OLD
binary -- always check the binary timestamp; (b) MHD's x1 sweep is "mhd_flux1", hydro's
"hflux_x1"; (c) the parser refuses a cmdline flag not present in the input file.

## 2026-09-06 23:30: down to the cost of a plain run (commit after cbf38d20)

Levers, in order of what they bought (20 cycles, n128, table; plain 5.22 s): ThermoAt
(one EvalNoMu for e,p,chi_rho,chi_T,cv; walk returns its end pressure) poly 14.6 -> 8.7;
5-evaluation polytropic fast path (predictor / one eval at the predicted end / trapezoid
corrector, p and e moved first-order to the corrected d) 8.7 -> 7.0 (KE 5.4e30 -> 6.2e30;
the plain frozen-coefficient version was 6.7 s but KE 1.16e31); wb_cache_every=1 6.15,
=10 5.37 (KE 7.0e30 -- vs 2.46e32 off). isentropic_dt = (d,T) RK2 isentrope, no root
find; adaptive_fast uses it: 5.45 s at every=10. The floor is the deviation-PLM logic
itself (~0.3-0.5 s here, the ideal gas pays it too). OPTIONS: <hydro|mhd>/wb_option =
polytropic|adaptive_fast|isentropic_dt (+ the old four), wb_cache_every = N (0 = every
stage). Production recipe: wellbalance_dynamic=true, wb_x1=true, wb_option=polytropic,
wb_cache_every=10. UNTESTED: GPU build, MHD path, restart (the cache is rebuilt from w0,
so restarts need nothing), AMR (wbq0 is realloc'd with nmb? -- NO: allocated once in
the ctor at nmb; AMR/load balance would need a realloc hook).

## 2026-09-06 23:50: WORKS for hydro, MHD, resistive MHD on cs and sp (commit after adec51a5)

The ghost-FACE potentials were never filled by the pgen (active faces only): the
outermost ghost cell's background saw phi -> 0, the polytropic walk overflowed, every
MHD run with it NaN'd at cycle 2 (hydro survived by luck). Fixed in the wbgrav kernel
(all faces) + WBGuard (flatten a non-finite/non-positive background to the anchor).
Column at rest, 20 cycles, polytropic, wb_cache_every=10, EOS resistivity on:
cs hydro 23.2 -> 23.6 s (KE 2.46e32 -> 6.99e30, == sp to 4 digits), cs MHD bbot 1e-6
37.6 -> 38.3, cs MHD bbot 3 G 38.2 -> 38.9 (same KE), sp MHD 3 G 7.50 -> 7.70 (6.52e30).
Bisection trap: NaN in MHD but not hydro with the SAME walk -> look at what differs in
the INPUTS to the walk (ghost states), not the walk. Test dirs $S/wbmhd, $S/wbcs.

## 2026-09-07 00:10: Cartesian x1/x3 + the closure family fixed (commit ea0cd573)

New pgen src/pgen/wb_column.cpp + inputs/tests/wb_column.athinput (build with
-D PROBLEM=wb_column, build_wbcol): plane-parallel column, problem/axis = 1|3, g0, ap
(point-mass-like g if > 0), linear T(z) (tgrad), b0 (uniform B_x2), user walls that
continue the column into the ghosts (reflect walls alone generate a wall disturbance
that fills the box by 200 cycles and hides the balance: interior v was 1e-15 at cycle 1
with the scheme on, 1e-3 by cycle 200 from the walls). 64 vertical cells, 200 cycles,
interior rms v: point-mass g, hydro/MHD/resistive x1 and x3 ALL 4.6e-5 -> 3.7e-6;
constant g 3.1e-5 -> 1.4e-10; isothermal+isothermal closure KE 1.7e-11 -> 2e-28.
The polytropic closure is now T LINEAR in Phi (slope from the neighbours): for an ideal
gas the exact polytrope d ~ T^(-(1+b)/b), containing isothermal AND isentropic exactly
(the ln T form contained neither; it gave only 6x on this column). sp dhj column:
8.46e30 vs 6.99e30 (ln T) -- slightly worse there, kept for exactness on polytropes.
x3 uses the UNCACHED walk (1.2-1.4x cost); the cache is x1-only.
TRAP: a file-scope DvceArray in a pgen aborts at exit (destroyed after Kokkos::finalize)
-> release it in pgen_final_func.

## 2026-09-07 00:40: RT, rotation, GPU

RT on, rotation off (sp, 60 cycles): deep v_r 4.1e-2 -> 6.7e-3 km/s, WB helps as in the
column tests. Rotation on: WB on gives deep v_r 4-5e-2 vs off 1.9e-2 (n128), KE1 6x, and
mass +1.3e-3 vs -1.6e-3 -- the SAME with every closure and cache setting, fine at omega/10.
RESOLUTION SETTLES IT (rotation on, t=1700 s, deep mean v_r): off n64 -1.1e-1, n128
+3.3e-3, n256 +3.0e-2; on n64 +2.9e-2, n128 +3.7e-2, n256 +3.9e-2. The plain scheme
CONVERGES TOWARD the WB answer: the outward drift is the physical initial adjustment of
the (non-rotating-hydrostatic) IC to the centrifugal force (Omega^2 r ~ 6e-3 g), which
the plain scheme damps/under-resolves below n256. Long runs: both settle to +-1e-2
oscillations after ~3000 s. A Cartesian body-force reproducer (wb_column gsrc_fac=1.006,
ideal gas) shows WB-on 5x KE too, but its walls (fixed column = a reservoir; reflect =
3.8 % mass drift on this pgen) make it useless as a gate -- do NOT use it for that.
GPU: build_dhj_gpu rebuilt (00:16, HEAD), bench/wb_gpu_check{,_off,_1rank,_every0,_iso}
(cs 16x16, 20 cycles): GPU == fresh CPU to 5e-12 in rho/e and 1e-4 in the (tiny) residual
v_r for scheme off, polytropic, isothermal, cache every stage; 1 rank == 2 ranks BITWISE.
The first comparison showed 1.4e-5 / KE1 19 % -- that was a STALE CPU REFERENCE (run
before the T-linear closure + ghost-face fix); a CPU rerun matched the GPU. TRAP: re-run
the CPU twin on the SAME binary before blaming the device. Builds: an
ambiguous getWBq0 overload (defaulted temps on the member wrapper) broke HIP/clang while
gcc accepted it -- the wrapper now forwards -1s without defaults.

## 2026-09-07 01:00: rotation in the potential (problem/rot_potential)

TotPotAt = GravPotAt - Omega^2 (r sin theta)^2/2; ZEffFromPot inverts Phi_grav(z) so the
1D column is sampled at the equivalent height of each cell's TOTAL potential ->
barotropic equilibrium in r AND theta (IC and ghosts, one kernel). Radial centrifugal
force dropped from the explicit source when the WB scheme is on (rotpot_src); theta part
+ Coriolis explicit; energy: with etotgrav the flux carries the potential (radial work
term dropped), without it the explicit work stays. cs branch drops the radial projection
Omega^2 (x rh0 + y rh1). Rotating column at rest, 1700 s, deep band: WB on: old IC v_r
4.3e-2 / v_th 2.0e-2 -> rot_potential 3.4e-3 / 4.3e-3 (sp), 3.9e-3 / 2.7e-4 (cs), KE1 14x
lower. PLAIN scheme + rot_potential drifts inward 3.6e-2 (worse than plain + old IC): the
switch is for wellbalance_dynamic runs. Not restart-compatible with runs started without.
PRODUCTION RECIPE NOW: wellbalance_dynamic=true, wb_x1=true, wb_option=polytropic,
wb_cache_every=10, problem/rot_potential=true (from scratch only).
