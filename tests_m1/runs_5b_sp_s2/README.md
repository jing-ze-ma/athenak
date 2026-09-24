# runs_5b_sp_s2: implicit M1 with the chi(f) closures on a spherical-polar wedge (stage S2)

Stage S2 of docs/dev/rad_m1_curvilinear_design.md (branch m1-curv-design, 2498ba73:
sect. 1.4 radial face with the integrating factor, 1.5 sp transverse faces, sect. 6 S2 row,
tests T-S2, T-S4, T-S5, risk 6), on the S1 wedge (no poles; runs_5a_sp_s1).
Branch m1-sp2 from rt-integration 4e428a30 (S1 merged).

Binaries, run trees and logs: /viper/ptmp2/jinma/s2_0924 (`bin/`, `cpu/`, `gpu/`,
`gate_*.log`).  ref = `git archive 4e428a30`, new = snapshot of this branch (`scripts/`
here: `build.sh`, `snap_new.sh`, `run.sh`, the gate and run scripts).  CPU: gcc 14 +
openmpi 5, Release, MPI.  Analysis: `s2lib.py` (loaders, the radial grid as
CoordSphericalPolar builds it), `ts2.py`, `ts4.py` (with the 1-D spherical reference),
`ts5.py`, `tsym.py`, `tr6.py`, `spread.py`.

## What changed

| file | change |
| --- | --- |
| `src/rad_m1/rad_m1_implicit.hpp` | `M1SphDrr`: the radial face weight of a cell with the integrating factor, `(1-chi)/2 + (3chi-1)/2 n_r^2 (r_c/r_f)^2`.  `M1SphCurv`: the LAGGED rest of `(div P)_d` at a cell centre in physical (r, theta, phi) components: always the curvature of the diagonal components (`-q(1-n_r^2)/r` radially, `cot(P_tt-P_pp)/r` on theta faces), and with `implicit_offdiag = lagged` the off-diagonal terms with their spherical factors (`(1/r) d_th P_rt + cot P_rt/r + (1/(r sin)) d_ph P_rp`; `(1/r^3) d_r(r^3 P_rt) + (1/(r sin)) d_ph P_tp`; `(1/r^3) d_r(r^3 P_rp) + (1/r) d_th P_tp + 2 cot P_tp/r`), over the coordinate distances of the neighbours, one-sided at physical boundaries as M1OffDiv. |
| `src/rad_m1/rad_m1_implicit.cpp` | `sph_q` overwrites inside the S1 `sph` blocks only: x1 row (`m1_impl_asm`) and x1 face flux (`m1_impl_face`) take `M1SphDrr` for both cells of each face and the face mean of `M1SphCurv(d=0)`; the x2/x3 face fluxes (`m1_impl_f2face/f3face`) overwrite `off` with the face mean of `M1SphCurv(d=1,2)`.  Every Cartesian expression is textually that of 4e428a30; the only removed lines are the two S1 sp lines `bb += nu{p,m}*df*wi`, now `wiu`/`wil` (equal to `wi` unless `sph_q`). |
| `src/rad_m1/rad_m1_sph.cpp` | `SphericalS1Check`: closure = m1, minerbo, kershaw allowed (vet_sc, tau refused); `sph_q = !eddington`; for these closures `implicit_offdiag` = auto -> none, lagged allowed, operator refused (Eddington keeps S1: none only). |
| `src/rad_m1/rad_m1.hpp`, `rad_m1.cpp` | `sph_q`; the constructor comment and fatal text. |
| `src/pgen/tests/rad_m1_tests2.cpp` | `sph_shell`: `e_in` (Dirichlet source cell, not used by the final gates).  New `m1_test = sph_atm` (T-S4, risk 6): gas `rho = atm_rho0 (r/r_in)^-atm_rho_n`, `atm_hold` (re-impose rho and v every stage, keep e_int), `atm_init = eddington` (spherical Eddington E, `F = F_in (r_in/r)^2` on cells and x1 faces, gas at `(E/a)^(1/4)` when kappa_p > 0), `atm_seed` (cell-to-cell perturbation of E), and for `force_reference = wb_arad` the reference `atm_aref = face` (default: mean of `kappa_t F_f/c` over the two x1 faces with the exact `F_f`) or `cell` (`kappa_t F(r_c)/c`). |

F stays in physical orthonormal components; the force deposit (`m1_impl_wb`) is unchanged.

## Gates

### T-bit (ref 4e428a30 vs new, bitwise = every output file, rst from `<par_end>`)

CPU, the 11 S1 cases (`scripts/gate_tbit.sh`, logs `gate_tbit.log`, `gate_tbit2.log`, the
second on the final binaries): 10/11 BITWISE (slab 1 and 2 ranks, new defaults, vet_sc
full, vimp, box3d old/new defaults, Marshak x1 and 3-D, sp hydro blast).  expl2 (slab,
transport = explicit): all 8 data files identical, the 2 restart files differ; this is
run-to-run: the explicit-transport restart comes out in one of two layouts (435877 or
435965 bytes) even for the same binary (ref vs ref rerun DIFFERENT, `cpu/t_expl2_r3`,
`_r4`), and every ref-vs-new pair with the same layout is BITWISE (new vs r3, n3 vs ref,
nn vs rr).  Pre-existing, not touched.

S1 sp Eddington gates, ref vs new (`scripts/gate_s1sp.sh`, `gate_s1sp.log`): T-S1
uniform and stretched, T-sym (200 cycles, 4 ranks), T-sym gas (200 cycles), Marshak
shell: 5/5 BITWISE.

GPU (apudev job 11957764, 2 x gfx942, rocm 6.3, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1,
`scripts/gate_gpu.sh`, log `gpu/log.out.11957764`): He slab 1 rank, box3d_nd 2 ranks,
box3d_nd + vet_sc 2 ranks: 3/3 BITWISE (10/10 files each).  Built before the last,
pgen-only commit (`atm_aref`, a branch no Cartesian run enters); the CPU gate was rerun
on the final binaries.

### T-S2 free-streaming point source (kappa = 0, M1, `inp/sph_fs*.athinput`)

r = 1..3, 4 x 4 wedge, imposed F_in = 1 at r_in, vacuum Marshak outflow with
marshak_q = 1 (F = c E_ie), implicit_cfl 1e4, 60 steps.  `ts2.py`; runs `cpu/F_fs*`.

| grid | n | max abs(r^2 E / mean - 1) | F/(cE) - 1 | order |
| --- | --- | --- | --- | --- |
| uniform | 32 | 7.4e-13 | 2.101e-2 | |
| | 64 | 5.9e-12 | 1.046e-2 | 1.006 |
| | 128 | 7.9e-12 | 5.220e-3 | 1.003 |
| he4 stretch | 32 | 2.0e-12 | 1.026e-2 | |
| | 64 | 3.7e-12 | 5.022e-3 | 1.030 |
| | 128 | 4.3e-12 | 2.481e-3 | 1.018 |

r^2 E is constant to the Picard floor (the solves stop at 2e-11 per pass).  F/(cE) - 1 is
evaluated with the face value E_f = (r^2E)/r_f^2 and r_f^2 F_f = r_in^2 F_in; it equals
(r_out/r_ie)^2 - 1 to all printed digits, i.e. the whole error is the Marshak end cell
using the cell E half a cell inside r_out (first order).  |F_theta|, |F_phi| / |F_r| <=
7e-27.  A first version (Dirichlet source cell, reflecting theta) diverged with
implicit_offdiag = lagged (Findings 3) and converged with `none`, the sp default.

### T-S4 grey extended atmosphere, gas coupled, M1 (`inp/sph_atm*.athinput`)

r = 1..5, rho = (r/r_in)^-2, kappa_t = kappa_f + kappa_s = 100, kappa_p = kappa_e = 20
(tau = 80), c = 100, a = 1e8, F_in = 1, Marshak q = 0.5 at r_out, 4 x 4 wedge periodic in
theta and phi, gas held (`atm_hold`), Eddington start, closure frozen per step
(`implicit_closure_lag = step`), lin/Picard tol 1e-10, 800 steps (steady: the dumps at 600
and 800 agree to all printed digits), a 1e-6 transverse seed (`atm_seed`).
Reference: `ts4.py` integrates `d(chi E)/dr = -(3chi-1)E/r - rho kappa F/c` with the same
chi(f) inward from F = c q E.  Runs `cpu/sa*`, `cpu/sas*`, `cpu/sae32`.

| run | n | L1(T_rad/T_ref - 1) | Linf | max abs(L_f/L_in - 1) | max abs(T_gas/T_rad - 1) |
| --- | --- | --- | --- | --- | --- |
| M1, uniform r | 16 | 1.93e-2 | 1.57e-1 | 9.4e-5 | 5.0e-6 |
| | 32 | 9.12e-3 | 9.98e-2 | 5.6e-5 | 5.4e-6 |
| | 64 | 4.40e-3 | 5.99e-2 | 3.0e-5 | 5.6e-6 |
| M1, he4 stretch | 32 | 7.88e-3 | 5.90e-2 | 7.0e-5 | 5.6e-6 |
| | 64 | 3.71e-3 | 3.29e-2 | 3.6e-5 | 5.8e-6 |
| Eddington (S1) vs its own reference | 32 | 8.99e-3 | 7.76e-2 | 5.6e-5 | 4.2e-6 |

L1 order 1.08, 1.05 (uniform), 1.09 (stretched): first order, set by the Marshak surface
(largest error in the top cell).  Richardson extrapolation of E from n = 32, 64 matches
the M1 reference to L1 2.6e-3 (1.3 < r < 4.4) and the Eddington reference only to 1.7e-2;
the two references differ by L1 6.4e-3, Linf 2.6e-2 in T.  The luminosity is rebuilt
from the cell means (F_{i+1/2} = 2 F1_i - F_{i-1/2}); its 3e-5..9e-5 is the enthalpy flux
of the held gas, which the hydro stage moves before the hold resets it (same in the
Eddington run).  The transverse seed decays: max over (theta, phi) of the E spread 2e-6 at
t = 0, 0 (float32 dump) at cycles 400 and 800 in every run (`spread.py`).

### Risk 6: radiation force vs a well-balanced reference, gas momentum on

`inp/sph_atm_r6.athinput`, `scripts/r6v3.sh`: uniform rho = 1, pure scattering kappa_s =
25 (tau = 100), T_gas = 1 (c_s = 1.3, c = 100), 64 x 4 x 4.  Stage A: gas held, 1000
steps to radiative equilibrium (L_f const to 2.5e-10); A3: 10 held steps at hydro CFL
0.3 (a restart keeps the last dt); B3: gas free, 100 steps (t = 1.08), sp hydro with
sp_wellbalanced_src.  No gravity: the reference acceleration stands for the external
well-balanced source (the static background is exact for the hydro, control run).

* The deposited force at the equilibrium (face mean of kappa_t F_f/c) equals the face-form
  reference kappa_t F_in r_in^2 (r_l^-2 + r_r^-2)/(2c) to 2.5e-10; against the cell-centred
  kappa_t F(r_c)/c it differs by max 4.0e-3 (inner cell), rms 1.2e-3: O((dr/r)^2).  One
  step with wb_arad: v_r = 5e-13 against 0.94 with force_reference = none.
* 100 steps, rms and max |v| (`tr6.py`, runs `cpu/r6B3_*`):

| run | rms abs(v) | max abs(v) | max abs(T/T0 - 1) |
| --- | --- | --- | --- |
| force_reference = wb_arad (face form) | 2.46e-5 | 5.8e-5 | 2.0e-4 |
| force_reference = none | 5.11e-2 | 1.29e-1 | 5.4e-2 |
| control, no radiation force | 3.1e-16 | 1.2e-15 | 0 |

  The wb_arad flow is not the momentum residual: it is the coupling's energy
  bookkeeping (rad_m1_coupling.cpp header, "MOMENTUM residual, WORK full"), which adds
  W = vbar dm = dm^2/(2 rho) per step with the FULL force while the momentum gets the
  residual; with no external source returning rho a_ref that heat is real.  Estimate:
  dt = 0.011, kappa F/c = 0.25 at r_in, W = 3.6e-6 per step, T/T0 - 1 = 2.4e-4 after 100
  steps; measured 2.0e-4, and it drives the 2.5e-5 flow.  With the stage-A time step
  carried into the first free step (dt = 4) the same W heated the inner cell by 29 % in
  one step (`cpu/r6C_wb`); a first version with cold gas (T = 1e-4, dt = 1.2) reached rms
  |v| 5.4e-2 in 100 steps for the same reason (`cpu/r6B_wb`).

### T-S5 radiation-modified acoustic wave (T10) in a thin shell (`inp/rw_*.athinput`)

The runs_3b8 radwave (M1, c/c_s = 1e3, kappa rho = 1e5, one wavelength = 1), run along
theta in a shell r = R -+ 0.05 (4 x 64 x 4, theta range 1/R, phi range 0.1/R, reflecting
in r), vs the same wave along x2 of a Cartesian slab 0.1 x 1 x 0.1 whose x2 origin gives
the IC the same phase sampling (the IC phase is absolute; unaligned twins differ by 3e-4
from IC sampling alone).  Two periods; `ts5.py` (Fourier amplitude Z of the density
along the slice, phase speed and decay from fits).  Runs `cpu/F_rws*`, `cpu/F_rwc*`.

| R | abs(Z_sp/Z_sp(0) - Z_c/Z_c(0)) / abs(Z_c/Z_c(0)) at t = 2 periods | v_phase rel. diff | decay diff (of 0.1106) |
| --- | --- | --- | --- |
| 10 | 1.60e-3 | +3.6e-5 | +5.9e-4 |
| 100 | 4.97e-4 | -8.6e-6 | +1.5e-4 |
| 1000 | 4.97e-5 | +1.3e-6 | +6.0e-5 |
| 10000 | 3.98e-6 | +2.0e-7 | +1.9e-7 |

1/R from R = 100 on (x10 per decade), 3.2x from 10 to 100.  implicit_offdiag = lagged at
R = 100: 4.975e-4 (none: 4.974e-4).

### T-sym and restart with the M1 closure (the S1 inputs, closure = m1, precond = line)

`sph_sym` (stretched r, tau 2000, reflecting theta, 4 ranks, 500 steps; E changes by 0.80
of its peak): max spread(E)/E over six (j, k) slices 2.4e-15 with lin tol 1e-14, 8.5e-10
at the default tolerances, 8.3e-10 with closure = kershaw (`cpu/sym_m1d`, `sym_m1a`,
`sym_ker`, `tsym.py`, final binary).  `sph_sym_gas` with
c = 100 (absorbing, moving gas): E 2.5e-14, rho 5.2e-15, e 6.4e-15, v_r 3.8e-15, v_theta
and v_phi / v_r <= 5e-16 (`cpu/symg_m1`).  (With the input's c = 1 the gas reaches v ~ c,
f -> 1 and the O(v/c) system fails, NaN at cycle ~230; Eddington survives that only
because chi is fixed.)  Restart: `sph_sym_gas_rst` (c = 100), 4 ranks, restarted at cycle
250: 24 of 24 tab slices and both later restart payloads identical (`cpu/rstA_m1`,
`rstB_m1`).

### Refusals (all rc 1, `cpu/fatal_*`)

closure = tau, vet_sc (periodic theta), implicit_offdiag = operator, eddington + lagged,
transport = explicit, a theta range reaching 0.  cpplint (CI filter) clean on every
changed C++ file.

## Findings outside S2 (pre-existing, Cartesian too)

1. **Transverse instability of the multi-D implicit M1 (chi(f) closure, central flux) in
   optically thin cells.**  Cartesian reproducer: `sph_atm_bin` with
   use_spherical_polar = false, kappa_s = 25 (rho kappa = 25/r^2), atm_seed = 1e-6
   (`cpu/zc_base`, `scripts/zt12.sh`): the (y, z) spread of E grows from 2e-6 to 6e-3 in 10
   steps and to O(1) by 20.  Same with 10x smaller dt, with dbg_trans_memory = 0, with
   closure_lag = pass or relax 0.3, with Kershaw, with c = 1000, with 4x wider transverse
   cells; gone with kappa x4 (tau_cell >= 0.5 at the surface) and with Eddington.  The
   1-D radial problem is stable (a Python replica of the column scheme converges from a
   seeded start).  On the wedge it appears unseeded because the columns differ at
   round-off; on a Cartesian mesh the columns stay bitwise identical, which hides it.
   T-S4 is therefore run in the stable regime (kappa_t = 100).
2. **Reflecting transverse walls** make the thin M1 atmosphere unstable on both meshes
   (`cpu/z_scs_cr` Cartesian, F_2 -> c E), periodic walls do not.  The He wedge is
   periodic in theta and phi.
3. **Lagged off-diagonal terms** have no fixed point in thin cells at c dt >> dx:
   Cartesian with reflecting x2 walls diverges (`cpu/x_cartrefl`), as does the sp form
   (free-streaming source).  Hence `none` is the sp default; `lagged` is available.
4. **Linear tolerance**: implicit_lin_tol 1e-13 gives ~3 BiCGStab breakdowns and a
   line-Jacobi fallback per solve in the M1 atmospheres (both meshes); 1e-10 converges
   (Picard 3.6 passes mean).
5. **wb_arad energy bookkeeping** (risk 6 above): W = dm^2/(2 rho) per step with the full
   force; harmless at the explicit-hydro dt, not with an implicit-radiation-sized dt.
6. **Explicit-transport restart layout** is nondeterministic (two sizes), T-bit above.

## Still refused on sp

The poles (polar boundary, theta range reaching 0 or pi), the cubed sphere, the theta
stretch; closure = vet_sc and tau; implicit_offdiag = operator (M1OffDiv on mb_size, the
od cache and the 19-point stencil are not converted); transport = explicit;
time_scheme = hesdirk2; more than one MeshBlock along x1; SMR/AMR;
implicit_halo_mpi = true; implicit_flux != central (the x1 HLL/blend uses mb_size
`tauf`); implicit_recon != dc; implicit_trans_limit != none; implicit_vimp; dbg_tensor.

## HANDOVER

State: m1-sp2 = rt-integration 4e428a30 + one commit (code, inputs, scripts, this README),
NOT merged.  Final binaries (md5): new none 2edd0f1e, new box_convection 73c1e037, new
sp_test 44bbe0b5, new box_convection GPU 2e355982 (pre pgen-only change), ref none
6b9f9094, ref box_convection 944a30d1, ref GPU d12c9663.  Build trees deleted, binaries
kept in /viper/ptmp2/jinma/s2_0924/bin.

Next steps, in order:

1. Finding 1 decides whether M1 (chi(f)) can be used in the thin top of the He wedge at
   all.  Look first at the transverse face flux of a chi(f) cell (M1SphTransRow /
   m1_impl_tcell and the closure's n from the lagged cell F2, F3): the instability needs
   the multi-D closure and thin cells, is dt-independent, and is absent in 1-D.  The
   Cartesian cure is S2-independent; the wedge takes it unchanged.  Until then use
   Eddington (S1) or vet on the wedge.
2. **S5 `vet_col`** (design sect. 3, option D): per column a 1-D spherical short-
   characteristics formal solution (mu grid, radial axis projection) giving chi and n_r
   for the lagged closure; gate T-S6 (f_K(r) vs the 1-D reference, e.g. ts4.py's ODE with
   the exact sphericity factors) and T-S4 with vet_col.  vet_col makes the tensor a
   frozen VET per step (like Eddington in the solver), which is also the way around
   Finding 1 if it turns out to be the chi(f) feedback.
3. Cheap lifts, each with its own bitwise test: implicit_halo_mpi on the wedge,
   hesdirk2, implicit_flux = blend (x1 `tauf`/`wdc` over dxface), implicit_offdiag =
   operator (M1SphCurv as the od cache, stencil slots from the sp coefficients).
4. The poles (S1 HANDOVER item 2) and the tau closure on spherical columns.
