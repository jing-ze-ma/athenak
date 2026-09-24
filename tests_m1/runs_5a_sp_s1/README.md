# runs_5a_sp_s1: implicit M1 on a spherical-polar wedge (stage S1, branch m1-sp)

Stage S1 of docs/dev/rad_m1_curvilinear_design.md (branch m1-curv-design, 2498ba73),
with the scope change of 2026-09-24: a WEDGE clear of the poles.  The pole face, the polar
ghosts through the slot-role table and the polar-row Krylov counts are NOT in S1.
Base: rt-integration 90b01f2f (S0 = 22418492 is in it: rad_m1 geometry fatal and the
per-exchange vector-slot table `SetVectorPairs`).

Binaries, inputs, scripts and run trees: /viper/ptmp2/jinma/s1_0924 (`build.sh`,
`snap_new.sh`, `run.sh`, `gate_tbit.sh`, `cmp.py`, `inp/`, runs in `cpu/`).  ref =
`git archive 90b01f2f`, new = snapshot of this branch.  CPU: gcc 14 + openmpi 5, Release,
MPI.  Analysis scripts here: `s1lib.py` (tab/bin loader, the radial grid exactly as
`CoordSphericalPolar` builds it), `ts1.py`, `ts3.py`, `tsym.py`.

## Target grid (He presupernova model, bench/wt_he4, read only)

`inputs/hydro/he4_presn_sp.athinput` (the sp variant of he4_presn_cs, 2026-09-19):

* radial: nx1 = 96, r = 1.18585e11 .. 2.4057e11 cm (0.50 R .. the tau = 1e-2 level),
  **stretched**: `use_grid_stretch_r_poly`, c1..c4 = +0.134559, +4.942505, -8.950106,
  +4.403042; meshblock nx1 = 96 (one radial block per column);
* theta: nx2 = 256 over pi/2 +- 45 deg, **periodic** (a wedge, no poles);
* phi: nx3 = 256 over 90 deg, periodic; meshblocks 96 x 32 x 32 (64 blocks).

So the target is exactly what S1 now supports (wedge, stretched r, one block in r).  The
stretch made the per-cell spacing part of S1.

## What changed

| file | change |
| --- | --- |
| `src/rad_m1/rad_m1.cpp` | geometry guard: sp is no longer refused outright; still refused: the polar boundary and any theta range reaching 0 or pi ("the poles are not supported yet"), the cubed sphere, the theta stretch, the radial stretches on a Cartesian mesh, relativity.  Sets `sph_geom` and calls `SphericalS1Check` after `ImplicitInit`. |
| `src/rad_m1/rad_m1_sph.cpp` (new, in `src/CMakeLists.txt`) | `SphericalS1Check`: on sp everything but transport = implicit, closure = eddington, time_scheme = be, one block along x1, no SMR, implicit_halo_mpi = false, implicit_flux = central, implicit_recon = dc, implicit_trans_limit = none, no implicit_vimp, no dbg_tensor is a fatal.  implicit_offdiag is set to none (Eddington D_ab = 0 off the diagonal: the dropped terms are identically zero, and M1OffDiv's uniform dx stays off this mesh). |
| `src/rad_m1/rad_m1.hpp` | `sph_geom`, `SphericalS1Check`. |
| `src/rad_m1/rad_m1_implicit.cpp` | the geometry, in separate `if (sph)` blocks that OVERWRITE the Cartesian results; every Cartesian expression is textually that of 90b01f2f (the first version, 600c4fed, used per-face run-time selects in the Cartesian expressions and was NOT GPU bitwise, see below): x1 row (`m1_impl_asm`): per face `nu = dt A1_f/V_i` instead of `dt/dx1`, and the face-flux gradient over `dxface.x1f` (centroid to centroid) instead of `dx1`, including the Marshak / flux boundary faces (the imposed flux enters times its area); x1 face flux (`m1_impl_face`): gradient over `dxface.x1f`; x2/x3 face fluxes (`m1_impl_f2face/f3face`): gradient over the arc lengths `dxface.x2f = r dtheta`, `dxface.x3f = r sin(theta) dphi`; transverse cell terms (`m1_impl_tcell`, rebuilt by `M1SphTransRow`): `dt A2_f/V`, `dt A3_f/V` per face in TDIA / TRHS / CJM..CKP.  The Krylov operator, the stencil, the line preconditioner and PCR read those stored coefficients, so they need nothing.  `implicit_halo_mpi` defaults to false on sp. |
| `src/rad_m1/rad_m1_newdt.cpp` | on sp `implicit_cfl` uses the smallest physical width (Coordinates dx1, r dtheta, r sin(theta) dphi), not mb_size (angles, unstretched dr). |
| `src/pgen/tests/rad_m1_tests2.cpp` | `m1_test = sph_shell` (built-in `pgen_name = rad_m1_beam`): static uniform gas, `E = e_out + e_amp exp(-((r-r0)/w)^2)`, F = 0; also runs on a Cartesian mesh (r = x1) as the planar twin. |

F is stored in physical orthonormal components (r, theta, phi), the hydro momentum's
basis, so the radiation force deposit (`m1_impl_wb`: face flux times the face opacity,
half to each cell) needs no change on the wedge.  The cell F is the mean of its two face
fluxes as on the Cartesian mesh.

## Gates

### T-bit (Cartesian M1 and sp hydro, ref vs new, CPU)

All 11 cases bitwise; table at the end.

### T-S1 steady spherical diffusion (Eddington, rho kappa = 100 scattering, thick)

`inp/sph_diff.athinput`: r = 1..3 (tau = 200), wedge theta = pi/2 +- 0.1 (reflect),
phi 0..0.2 (periodic), 4 x 4 angular cells, imposed flux F_in = 1e-2 at r_in
(`implicit_bc_x1min = flux`), outer end cell held at e_out = 1 (`efix`), flat start,
implicit_cfl 1e6, 20 steps.  Analytic: `E = e_out + b (1/r - 1/r_ie)`, `b = 3 k F_in
r_in^2/c`.  Face luminosity rebuilt from the cell means with F_in at r_in.
Double-precision radial slice (the bin writer is float32).

| grid | n | L1(E) | Linf(E) | max abs(r_f^2 F_f/(r_in^2 F_in) - 1) | order L1 |
| --- | --- | --- | --- | --- | --- |
| uniform r | 32 | 2.805e-5 | 9.84e-5 | 6.1e-14 | |
| | 64 | 7.021e-6 | 2.53e-5 | 1.3e-13 | 1.998 |
| | 128 | 1.757e-6 | 6.43e-6 | 3.9e-13 | 1.999 |
| he4 poly stretch | 32 | 1.610e-4 | 9.48e-4 | 7.7e-11 | |
| | 64 | 4.061e-5 | 2.75e-4 | 9.2e-12 | 1.987 |
| | 128 | 1.019e-5 | 7.34e-5 | 4.9e-12 | 1.994 |

max|F_theta|, max|F_phi| / max|F_r| <= 1.7e-13 in every run.  Runs `cpu/sd{32,64,128}`
(line_jacobi) and `cpu/sdp{32,64,128}` (bicgstab defaults), `ts1.py`.

### T-S3 Marshak wave in a thin shell

`inp/marshak_cart.athinput` / `inp/marshak_sph.athinput`: the I4 Marshak input
(inputs/tests/rad_m1_marshak.athinput, 128 cells, t = 202) with transport = implicit,
closure = eddington, implicit_cfl 10, Marshak ends; Cartesian with 4 x 4 transverse cells,
sp shell r = R..R+1 with a 4 x 4 wedge (dxmin = dr, same dt).  L1(E) against the S_N
reference `runs_3a/t6_ref_sn.txt` (`ts3.py`):

| run | L1(E) vs S_N | shell vs Cartesian L1 |
| --- | --- | --- |
| Cartesian | 0.01425 | |
| shell R = 100 | 0.01151 | 2.70e-3 |
| shell R = 1000 | 0.01398 | 2.71e-4 |
| shell R = 1e4 | 0.01422 | 2.72e-5 |

The shell reproduces the Cartesian gate and the difference is the curvature, exactly
1/R.  (G3 of runs_3c, central dc with the M1 closure: 0.0143 at CFL 10.)

### T-sym (spherically symmetric state on the wedge)

`inp/sph_sym.athinput`: he4 poly-stretched r = 1..3, 64 x 8 x 8, theta = pi/2 +- 0.4
(reflect), phi 0..0.8 (periodic), 4 MeshBlocks on 4 ranks, Gaussian shell E (e_out 1e-2,
amplitude 1, r0 = 2, w = 0.2), rho kappa = 1000 (tau 2000), implicit_cfl 100, bicgstab
with the default fused path (rbgs_fwd), 500 steps; E changes by 0.80 of its peak.
Radial double-precision slices at (j, k) = (0|3|7, 0|7) (edge and interior theta rows,
both phi edges), `tsym.py`:

| run | max_i spread(E)/E | max abs(F_theta), abs(F_phi) / max abs(F_r) |
| --- | --- | --- |
| sym_d: implicit_lin_tol 1e-14, implicit_tol 1e-12 | **1.96e-14** | 2.1e-14, 5.2e-14 |
| sym_a: default tolerances (lin 1e-10) | 1.64e-11 | 2.9e-11 |
| csym_a: CARTESIAN twin (E(x1) only), default tolerances | 8.57e-11 | 1.0e-10 |
| symg: + absorption (kappa_p 100), gas_feedback, moving gas, lin 1e-14 | 7.5e-13 (gas rho, e, v_r: <= 2.4e-13; v_theta, v_phi / v_r <= 4.8e-14) | 8.6e-14 |

At round-off with a tight solve.  At the default tolerance the spread is the linear
solver's: the red-black preconditioner is not symmetric in (j, k), and the Cartesian twin
shows the same (larger) spread.  Inner BiCGStab iterations per solve: mean 1.76 (sym_a),
2.37 (sym_d), max 5.  (The polar-row counts of design sect. 7 item 2 are out of S1.)

### Restart bitwise on sp

`symg` configuration + a restart dump every 250 cycles (`cpu/rstA`, 4 ranks); restarted
from `sd.00001.rst` (cycle 250) to 500 (`cpu/rstB`): 12 of 12 final tab slices (m1 and
hydro_w) and both later restart files (payload after `<par_end>`) identical.

### Other checks

* implicit_halo_direct on vs off, 1 rank x 4 MeshBlocks, sph_sym 100 cycles: 6/6 slices
  identical.
* Refusals (all rc 1, `cpu/fatal_*`): polar boundary ("the poles are not supported
  yet"), a theta range reaching 0 without the polar flag (x2min = 0), cubed sphere, sp + transport = explicit, sp + closure = m1, sp +
  implicit_halo_mpi = true, sp + time_scheme = hesdirk2, sp with two blocks along x1,
  radial poly stretch on a Cartesian mesh.
* Style: cpplint (CI filter) clean on every changed C++ file; flake8 (90) on the scripts.

### Not a sp effect: BiCGStab breakdowns at c dt/dx = 1e6

T-S1 at implicit_cfl 1e6 with bicgstab: 34 breakdowns / 6 line-Jacobi fallbacks in 20
steps (answer unaffected, L1 as above).  The Cartesian twin (same input,
use_spherical_polar = false) gives 74 / 24.  At that step the operator's condition number
puts the attainable relative residual near implicit_lin_tol.  Pre-existing; not touched.

## T-bit table

ref = 90b01f2f, new = m1-sp, same inputs and ranks, bitwise = every output file identical
(bin/rst from `<par_end>` on), `scripts/gate_tbit.sh` + `gate_tbit2.sh`, logs
/viper/ptmp2/jinma/s1_0924/gate_tbit*.log.  Inputs of the M1 He cases:
/viper/ptmp2/jinma/defaults_0923/cpu (the runs_4b/4c set).

| # | test (CPU) | ranks | result |
|---|---|---|---|
| 1 | M1 He slab 2-D `slab2d_old` (implicit, Eddington, bicgstab), 200 s, 1241 cycles | 1 | bitwise, 12 files |
| 2 | same | 2 | bitwise, 12 files |
| 3 | He slab 2-D `slab2d_nd` (new defaults) | 2 | bitwise, 12 files |
| 4 | slab + closure = vet_sc, vet_tensor = full | 2 | bitwise, 12 files |
| 5 | slab + implicit_vimp (s0_0923/inp/slab2d_vimp), 50 s, 310 cycles | 2 | bitwise, 10 files |
| 6 | slab, transport = explicit, 40 cycles | 2 | bitwise, 10 files |
| 7 | M1 3-D box `box3d_old`, 60 cycles | 4 | bitwise, 10 files |
| 8 | M1 3-D box `box3d_nd` (halo_mpi on), 60 cycles | 4 | bitwise, 10 files |
| 9 | Marshak I4, implicit_x1, closure m1, CFL 10, 2586 cycles | 1 | bitwise, 4 files |
| 10 | Marshak, transport = implicit, Eddington, 3-D 128x4x4 (`inp/marshak_cart`) | 1 | bitwise, 4 files |
| 11 | sp HYDRO: sp_test blast across the pole (polar boundary, 8x16x32), 79 cycles to t = 0.2 | 2 | bitwise, 7 files (hst + 6 bin) |

## HANDOVER

State: m1-sp = rt-integration 90b01f2f + 600c4fed + ced90b9a (the overwrite
restructure) + this README commit, NOT merged.  CPU numbers above: 600c4fed; after
ced90b9a the sp gates were rerun (`cpu/v2*`): T-S1 L1 identical to 5 digits (order 1.998 /
1.999 uniform, 1.987 / 1.994 stretched), T-sym 1.86e-14, symg 7.5e-13, Marshak R = 1e3
0.01398, restart 12/12 identical.

GPU Cartesian bitwise gate (apudev, 2 x gfx942, hipcc/rocm 6.3, HSA_XNACK=1,
HSA_NO_SCRATCH_RECLAIM=1, one job, interleaved; `scripts/gate_gpu*.sh`):
* 600c4fed vs 90b01f2f (job 11956456): He slab 2-D 1 rank bitwise; 3-D box3d_nd 2 ranks
  and box3d_nd + vet_sc 2 ranks DIFFERENT from the first hst row (time 1.1298793033168635
  vs ...184462): the per-face selects changed hipcc's FMA contraction.
* ced90b9a vs 90b01f2f (job 11956732): slab 1 rank, box3d_nd 2 ranks, box3d_nd + vet_sc
  2 ranks all BITWISE (10/10 files each); 90b01f2f rerun vs itself also bitwise (no
  run-to-run noise).

CPU T-bit of the merge: ref = rt-integration e89954e2, new = `git merge-tree
--write-tree ced90b9a e89954e2` (tree e9629506, no conflicts), the same 11 cases as the
table below (`scripts/gate_tbit3.sh`, log s1_0924/gate_tbit3.log): all 11 BITWISE.

Next step, S2 (design sect. 6), plus what S1 left refused:

1. M1 / Kershaw closures on the wedge: the radial integrating factor (sect. 1.4:
   `(r_R^2 Q_R E_R - r_L^2 Q_L E_L)/(r_f^2 dr_f)` in TA/TC and the lagged `-q(1-n_r^2)/r`),
   the lagged curvature by the Cartesian face sum (sect. 1.6, the `wb_geom` cache);
   gates T-S2 (free streaming 1/r^2), T-S4, T-S5.  The off-diagonal terms (`M1OffDiv`,
   stencil edge slots, od cache) still divide by mb_size and must take dxface there.
2. The poles: zero-area pole face, polar ghosts through `SetVectorPairs` (pbval_th pairs
   N2/N3, A2/A3, V2/V3 are already declared), and the polar-row Krylov counts at thick
   optical depth (design risk 2).
3. Cheap lifts, each needing its own bitwise test: implicit_halo_mpi on the wedge (no
   pole: the plain copy), time_scheme = hesdirk2 (rad_m1_time2.cpp has no dx), and the
   options with a uniform dx left inside: implicit_flux != central (`tauf`, `wdc` in
   `m1_impl_aphll`), implicit_recon = plm_dc, implicit_trans_limit (ImplicitTransTheta),
   implicit_vimp (ImplicitVimpBuild), closure = tau (rad_m1_tau.cpp).
4. implicit_enthalpy = plm (the default for Eddington) reconstructs on index space; on the
   stretched grid its deferred correction is not the exact non-uniform PLM (only matters
   with gas velocity; untested beyond T-sym's moving gas).
5. Design risk 6: the radiation force vs the hydro's well-balanced sp source in a
   hydrostatic radiative envelope (T-S4 with gas momentum on) before the He model runs.
