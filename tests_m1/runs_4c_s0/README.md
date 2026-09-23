# runs_4c_s0: rad_m1 geometry fatal + per-exchange component-role table (branch m1-s0)

Stage 0 of docs/dev/rad_m1_curvilinear_design.md (branch m1-curv-design, 2498ba73), sect. 0.
Ref = rt-integration be02c647, new = m1-s0 (cd117ac3 + b7be8e57). Binaries, inputs,
scripts and run trees: /viper/ptmp2/jinma/s0_0923 (`build.sh`, `run.sh`, `gate_cpu.sh`,
`gate_cpu2.sh`, `gate_dhj.sh`, `gpu/gate_gpu.sh`, `cmp.py`; runs in `cpu/`, `gpu/`).
Every build is a `git archive` snapshot (`src_ref`, `src_new`).

## What changed

1. **Geometry fatal** (`src/rad_m1/rad_m1.cpp`, the GEOMETRY GUARD block after the transport
   parse, ~l.70-107).  `<rad_m1>` (explicit and implicit transport) now refuses
   `mesh/use_spherical_polar`, `mesh/use_cubed_sphere`, `mesh/use_polar_boundary`, the
   radial stretches `use_grid_stretch_r` / `use_grid_stretch_r_poly`, `use_grid_stretch_theta`,
   and SR / GR / dynamical-GR coordinates.  The M1 kernels divide by the uniform
   `mb_size.dx1..3` and carry a Cartesian flux vector.  What the hydro allows: Cartesian;
   sp and cs; the radial stretches on sp and cs; the theta stretch on sp only; SR/GR.
   The stretches act only inside the sp/cs coordinate kernels
   (`coordinates.cpp` ~l.500, ~l.1596).  On a plain Cartesian mesh the hydro ignores them,
   and so would M1.  So a "stretched Cartesian" M1 run was not wrong: it was silently
   uniform.  It is refused anyway, so no stretched input reaches the uniform-dx kernels.
2. **Component-role table** (`src/bvals/bvals.hpp` `SetVectorPairs`, `vrole_set_`, `vrole_`;
   `src/bvals/bvals_cc.cpp` `SetVectorPairs()` after the constructor, and in
   `PackAndSendCC` the role lookup (`vec_`, `va_`, `vb_`) used by the seam transform switch
   (`cs_xform = vec_`), the polar flip (`if (vec_) signvar = -1`) and the seam transform
   loads (`a(m,va_,...)`, `a(m,vb_,...)`)).  Per `MeshBoundaryValuesCC` object:
   `SetVectorPairs(nvar, {{a,b},...})` marks each (a = x2 member, b = x3 member) as a
   tangential pair and every other slot as a scalar.  An object that never calls it keeps
   the (IVY, IVZ) rule unchanged (hydro, MHD, dyn_grmhd, radiation, z4c, two-stream,
   conduction, shearing box, ...).  The cs flux-seam path (`flux_seam_cc.cpp`, IM2/IM3)
   is only used by hydro/MHD and was left alone.
3. **rad_m1 tables**: `pbval_u` {(F2,F3)} (= today's rule), `pbval_th` (M1HaloCompT)
   {(N2,N3),(A2,A3),(V2,V3)} with F1 and KT in slots 2,3 as scalars, `pbval_tq` all
   scalar, `pbval_kr` scalar, `pbval_vm` {(DV2,DV3)} with the nine Jacobian rows P
   treated as scalars, and the vet_tensor = full D exchange (`rad_m1_vet.cpp`) with all
   six components scalar.  **Not done (later stage)**: the tensor components (D, and the
   directional P rows) would need a tensor transform at seams and poles, not a
   scalar copy.  They are latent while (1) forbids sp/cs.
4. **Unit test**: `src/pgen/bvals_roles_test.cpp` + `inputs/tests/bvals_roles_{cs,sp}.athinput`
   (`-D PROBLEM=bvals_roles_test`, nlim = 0, prints one `BVALS_ROLES ... PASS|FAIL` line
   and exits 1 on FAIL).  A 6-slot array (4 scalars + 1 pair) goes through four
   exchanges.  A uses the table {(4,5)}.  B uses the default rule with the fields permuted
   so the pair sits in (2,3), i.e. exactly how the hydro velocity is exchanged.  C uses
   the default rule on A's layout.  D uses the table {(2,3)} on B's layout.  Checks:
   A[v] == B[perm v] bitwise in every cell; D == B bitwise; ghosts were filled; and C
   differs from A in slots 2,3 and in slots 4,5, so the table is not vacuous.

## Gates

CPU: gcc 14 + openmpi 5, Release, MPI.  Bitwise = every output file identical (hst, bin,
rst, rt_profile/column dumps; bin/rst compared from `<par_end>` on, `cmp.py`).

| # | test (CPU) | ranks | result |
|---|---|---|---|
| 1 | M1 He slab 2-D (runs_4b `slab2d_old`: implicit, Eddington, bicgstab, halo_mpi off, so the pbval_th/tq/kr objects are used), 200 s, 1241 cycles | 1 | bitwise, 12 files |
| 2 | same | 2 | bitwise, 12 files |
| 3 | same + closure = vet_sc, vet_tensor = full (the D exchange) | 2 | bitwise, 12 files |
| 4 | same + implicit_vimp = true (pbval_vm), 50 s, 310 cycles | 2 | bitwise, 10 files |
| 5 | same, transport = explicit (pbval_u), 40 cycles | 2 | bitwise, 10 files |
| 6 | M1 3-D box 84x32x32 (`box3d_old`), 60 cycles | 4 | bitwise, 10 files |
| 7 | M1 3-D box, `box3d_nd` (new defaults: plm, halo_mpi on), 60 cycles | 4 | bitwise, 10 files |
| 8 | dhj HYDRO on cs, from cs_hyd4_prod/rst/dhj.00129.rst (t = 1.96725e7 s, rotation 64.5), read in place, prod input, tlim = t_rst + 80 s (5 cycles, to cycle 1003314) | 8 | bitwise (hst, rst, bin, log output) |
| 9 | dhj MHD on cs, from cs_mhd_prod4/rst/dhj.00097.rst (t = 1.47925e7 s, rotation 48.5), read in place, prod4 input, 5 cycles (to cycle 951489) | 8 | bitwise (hst, rst, bin, log output) |
| 10 | dhj on sp (`inputs/tests/dhj_ck_spherical.athinput`: sp + polar boundary + poly r stretch, ck two-stream), from scratch, 20 cycles | 2 | bitwise, 7 files |
| 11 | unit test `bvals_roles_cs` | 1, 4 | PASS: ghosts_filled 25344, table vs hydro-rule mismatches 0, (2,3)-table vs default mismatches 0; default rule would move scalars 2,3 in 14400 ghost cells and miss the pair 4,5 in 14400 |
| 12 | unit test `bvals_roles_sp` | 1, 4 | PASS: ghosts_filled 2560, mismatches 0 / 0; default rule: 768 / 768 polar ghost cells |
| 13 | M1 on cs, explicit and implicit (`bvals_roles_cs` + `<rad_m1>`) | 1 | FATAL "... this input sets: mesh/use_cubed_sphere.", rc 1 |
| 14 | M1 on sp, explicit and implicit | 1 | FATAL "... mesh/use_spherical_polar mesh/use_polar_boundary.", rc 1 |
| 15 | M1 slab + mesh/use_grid_stretch_r (Cartesian) | 1 | FATAL "... mesh/use_grid_stretch_r.", rc 1 |

cpplint (tst/scripts/style/cpplint.py) on bvals.hpp, bvals_cc.cpp, bvals_roles_test.cpp,
rad_m1.cpp, rad_m1_implicit.cpp, rad_m1_vet.cpp, with the CI filter
(`--filter=-build/include_subdir`, check_athena_cpp_style.sh): 0 errors.

GPU (apudev, 2 x gfx942, hipcc 6.3.4, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, ref and new
interleaved in one job): the same cs hydro / cs MHD restarts (2 ranks), the sp dhj
(20 cycles), and the 3-D M1 He box `runs_3k_gpu3d/he_slab_m1_3d.athinput`, 20 cycles, with
halo_mpi on (default) and off.  Script `gpu/gate_gpu.sh`, job 11953801, summary
`gpu/analyze_gpu.sh`.
