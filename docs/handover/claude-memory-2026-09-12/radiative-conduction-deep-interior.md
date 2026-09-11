---
name: radiative-conduction-deep-interior
description: What the dhj RT does below the 10 bar cut (NOTHING radiative; sigma T_int^4 injected at the cut from an infinite reservoir) and the fix built 2026-09-07 -- <mhd>/isotropic_conduction = radiative (Freedman Rosseland diffusion, rad_pcut_bar, rad_flux_inner at the wall, ck_int_at_cut=false) -- with the column test that verified it; the tau_R BLEND (fad2f5db, rad_tau_lo/hi) tested on the sp dhj column: dt+cost equal, overlap 0.4-2 bar, but ck (no H-) vs Freedman disagree 6x at 3700 K so the handover flux jumps 1.7x. NOT yet run on a full dhj
metadata:
  type: project
---

**Status quo (deep_hot_jupiter_rt).** Correlated-k two-stream only where p < ck_pcut_bar
(10 bar). At the cut (icut per column) the upward intensity = B(T_icut) + band-weighted
sigma T_int^4/pi (T_int from the Thorngren fit, ~540 K). Below the cut: src = 0, no
diffusion, no energy extracted -- F_int comes from an infinite reservoir, the deep T
profile keeps the IC, and OHMIC HEATING AT DEPTH HAS NO RADIATIVE OUTLET.

**Fix (commit after 68eb2250):** conduction module branch `radiative`: kappa_rad =
16 sigma T^3/(3 kappa_R rho), kappa_R = RosselandFreedman2014 (utils/rosseland.hpp, the
pgen's get_kapr now calls it), rad_pcut_bar (no flux on faces with p below it),
rad_flux_inner (erg/cm^2/s through the inner x1 wall, replaces the gradient there),
rad_flux_limit (F/sqrt(1+(F/sigma T^4)^2); the first form F/(1+F/sigma T^4) cut 43 % at
F = 0.76 sigma T^4 and 7 % at 0.08 -- REJECTED), rad_met, rad_kappa_fac. Needs <units>.
Pgen: problem/ck_int_at_cut = false -> no sigma T_int^4 at the cut (the conduction carries
it; the two-stream bottom is then just B(T_icut)). Dt: the module's limit, dx from
pcoord on sp/cs; cells above the cut impose none.
**Test** ($S/radcol3, inputs/tests/wb_column_radcond.athinput, wb_column tprof=1 =
radiative-equilibrium column via the same fit, cgs + <units> mu 2.3, 64 cells, 3000 K,
1084 bar bottom): F scaled to 3.5e8 with kappa_R x 1e-5 (nabla_rad 0.17 -> 0.02): discrete
face flux == F_int to 5e-4 (python replica); equilibrium column with conduction ==
without (steady, Etot drift 0); linear-T column responds 190 cm/s; 2x wall flux -> +2 K
bottom, Etot +1.9e-5 = exact energy budget. TRAPS: (a) physical F_int (3.5e6) is
invisible over 3e4 s (divergence 1e-7 of e) -- scale F and kappa_R together; (b) the
unscaled Freedman kappa_R with rho0 1e-2 gives nabla_rad 17 (T collapses in 3e6 cm) --
this test column is NOT a hot Jupiter interior, it is an operator test; (c) a fixed
ghost column already supplies F_int through the wall gradient, so 'no inner flux' is
NOT a gate -- use 2x instead.
**Consistency with correlated-k:** the ck table lacks H- bf/ff, so its Rosseland mean is
BELOW Freedman above ~3000 K (pgen note near ck_continuum); the two agree below. At the
10 bar cut on the dayside (3400-4000 K) the diffusion opacity is therefore higher than
the table's -> a gradient kink at the cut; the deep interior (6000-12000 K) is outside
BOTH (Freedman clamped at 4000 K). Next: H- in the ck continuum or a higher-T Rosseland
table (Freedman/Lupu), and a dhj run with rad_pcut_bar = ck_pcut_bar, rad_flux_inner =
sigma T_int^4, ck_int_at_cut = false (production input NOT yet switched).


## The tau_R blend (fad2f5db, 2026-09-07 02:00-02:50)

`<mhd>/rad_tau_lo = 30`, `rad_tau_hi = 300` (hi > 0 turns it on, replaces rad_pcut_bar).
Conduction::BuildRadWeights: one fixed-length top-down sweep per column (every cell,
every column, no divergence), tau_R with Freedman kappa_R and the PHYSICAL dr; face
weight w = raised cosine in log tau; diffusion flux x w; the dhj RT reads pcond->rad_w,
multiplies its LW face fluxes by 1-w, takes icut = first face with w < 1, and forces
ck_int_at_cut = false. The ck column dump carries tau_R(face) and w(face) as cols 9-10.

**Measured, sp dhj column at rest, 128x16x64 (2 blocks), 40 cycles, scratch blend/:**
- cost 83 s vs 82 s old; dt 28.9 s both (after the fix below).
- dayside IC: overlap at 0.37-1.3 bar, icut 30 vs 19 (10 bar). Evolved cs production dump
  rot 28: tau_R=300 face at 1.6-2.1 bar in EVERY column (1 cell spread -> 1 % warp loss);
  the RT chain is 103 cells vs 111 with the 10 bar cut, i.e. the blend makes RT CHEAPER.
- GATE FAILED on opacity consistency, not code: in the overlap F_lw(two-stream) = 0.37
  F_int but Freedman diffusion = 0.62 F_int (1.7x). RETRACTED my first reading ("ck lacks
  H-"): the ck CONTINUUM in the code DOES include John 1988 H- bf/ff from FastChem n(H-)
  (deep_hot_jupiter_rt.cpp ~1269; a stale comment above it said otherwise, now fixed).
  The k-table's LINE-ONLY Rosseland is 6x below Freedman at 3730 K/1 bar, but with H- and
  CIA the code's Rosseland is 1.0-1.8x ABOVE Freedman at 3500 K and 1.7-4.4x at 4800 K
  (docs/correlated_k_rt.md) -- the right direction and size for the 1.7x flux jump. So
  the inconsistency is the Freedman FIT (capped at 4000 K) vs the table+continuum. Fix:
  derive kappa_R FROM the ck table + continuum (Planck-derivative-weighted harmonic mean
  over 11 bands x 8 g) for the diffusion side -> handover consistent by construction.
  Cheap: per cell, the RT already looks up all 88 g-points. Not built. What we omit vs
  Exo-FMS: the H2- and He- free-free files in data/exo_fms_ck/cia (present, unread).

**Time-step bug fixed in the same commit:** Conduction::NewTimeStep used size.dx2/dx3
(ANGLES on sp) and the unstretched dx1 -> dt 9.4e-15 s on the dhj grid. Now pcoord's
dx1/dx2/dx3 in the radiative branch (constant/spitzer conduction still use size.dx: a
pre-existing curvilinear gap). In the blend the limit is on max(w_faces)*kappa: no
constraint where there is no flux; explicit dt_diff is 1.5e3 s at the top of the overlap
vs the 29 s hydro dt, so it does not bind at this resolution.

Cartesian column (build_wbcol, scratch radtau/): pcut mode == tau mode to 1e-8 (the top
face has tau = 0 hence w = 0 by construction, the only difference).


## kappa_R from the ck table (rad_kappa_src = table), 2026-09-07 03:00-03:20

Built as a START-UP TABLE, not a per-cell refill (the first version refilled every cell
every stage, +28 % cost; the user pointed out the opacity is a fixed function of (T,p)).
ck_build_rosseland_table() in the pgen: 38x34 grid of the ck table, Planck-derivative
band weights (centred +-1 % in T on the Planck-fraction table), harmonic over 8 g-points,
continuum incl. H- at rho = p mu m_H/kT (mu from the CE table); handed to
Conduction::rad_kr_tab/lT/lP; RosselandTable() log-bilinear lookup in flux, sweep, dt.
Cost = Freedman. Result: handover CONSISTENT (two-stream 0.36-0.38 F_int vs diffusion
0.40-0.43 in the overlap; the residual is the two-stream's 2/D vs 4/3 factor). Source
flat across the overlap. Commit after fad2f5db/03467b7f (see git log).

ck-table Rosseland / Freedman fit, printed at start-up: 1000 K 1.7 (1 bar) / 0.8
(21 bar); 2500 K 0.74 / 0.83; 3700 K 1.5 / 0.76; 4900 K 4.8 / 4.5. Above 4000 K the
Freedman FIT is clamped at 4000 K (RosselandFreedman2014 clamps T), so the 4.5x there is
mostly the fit's ceiling while H- keeps rising. Below 3000 K the +-30-70 % is line lists
+ the 8-point g-quadrature under-resolving the low-k tail (biases kappa_R HIGH); not
separated. For the blend only consistency with the two-stream matters, and that is now
by construction. Input: <mhd> rad_kappa_src = table (needs rt_ck).

## Overhead on the PRODUCTION inputs (bench/rad_ovh, 2026-09-07 04:00)

1000 cycles, 2 GPUs, diffusion + tau blend 30-300 + rad_kappa_src=table + rad_flux_inner
vs off: sp 27.7 -> 30.4 ms/cycle (+9.5 %), cs 40.5 -> 42.9 (+5.8 %); dt unchanged on
both (the w*kappa limit never binds). The cs guard in mesh.cpp now ALLOWS only the
radiative branch. Operator-on-grid audit: x1 uses x1v differences (stretch-aware) and the
update multiplies by pcoord areas / volume on both grids -> radial operator exact and
conservative; x2/x3 use pcoord dx2/dx3 (r dtheta, r sin theta dphi; cs arc lengths) but
the gradient is taken as if orthogonal -> on cs the g^{xi eta} cross term is dropped.
Horizontal radiative flux is ~1e-2 of radial at most, so that is a <1 % error on a
negligible term. Trap of the night: `pgrep -f "make -j 8"` matches the polling shell
itself -- the GPU build had NOT been running for 25 min. Poll a log marker instead.

## Exact cs operator + DEFAULT ON (commit after 34e9b1fd, 2026-09-07 05:00)

Face-normal derivative on cs: (dT/dl_xi - cos a dT/dl_eta)/sin a, cos_cell/sin_cell at
the face = mean of the two cells, transverse derivative centred over the two cells;
dt carries 1/sin^2. `rad_cs_exact = false` = the old orthogonal form (diagnostic).
TEST: cs_test iprob 15 (l=2 field, uniform rho at rest; inputs/tests/
cubed_sphere_raddiff.athinput; tst/test_suite/rad/test_rad_cs_raddiff_cpu.py): the
measured source vs -6 kappa (T-T0) r_fac: amplitude 0.991/0.998/1.0003 at 16/32/64 once
the expectation uses the code's ATOMIC MASS UNIT (not m_H: 4.8 %!) and the in-run decay
(1.5 % at 16). Old form: L1 24-30 %, Linf 1.7 at panel edges. sp_test iprob 15 is the
twin (hydro allowed for 15). Note the sp polar cells make the diffusion dt 1e-7 at
nx2=32 -> float32 dumps too noisy for sp at 32; cs is fine.
DEFAULTS: _eos (rt_ck now TRUE) and _blowup inputs carry the block; rad_flux_inner = -1
means sigma T_int^4 from Teq (pgen sets it, prints it). Grey inputs (hydro, ideal_xe)
untouched: the pgen REFUSES rad_tau_hi > 0 without rt_ck. Regression gate rewritten
(two-stream vs diffusion flux 0.90 in the overlap, 30 % gate), both rad tests PASS.
