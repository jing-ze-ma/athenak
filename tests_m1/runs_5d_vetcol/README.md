# runs_5d_vetcol: `closure = vet_col`, a per-column spherical VET (stage S5)

Stage S5 of docs/dev/rad_m1_curvilinear_design.md (sect. 3 option D, sect. 6 S5 row, tests
T-S6 and T-S4 with vet_col), on the S1/S2 wedge.  Branch m1-vetcol from m1-sp2 (b5f3afea).
The cubed sphere is not in this stage.

Run trees, binaries and logs: /viper/ptmp2/jinma/vetcol_0924 (`bin/`, `cpu/`, `gpu/`).
ref = the m1-sp2 binaries of runs_5b (/viper/ptmp2/jinma/s2_0924/bin/athena_new_*: none
2edd0f1e, box_convection 73c1e037, sp_test 44bbe0b5, GPU box_convection 2e355982, see the
S2 HANDOVER; symlinked as `bin/athena_ref_*`).  new = snapshot of this branch
(`scripts/snap_new.sh`, `build.sh`): none 11fb6d07, box_convection ef54daf3, sp_test
a70cc991, GPU none c88d205b, GPU box_convection 5ff26743.  The committed source adds one
`#include <utility>` (cpplint) after the snapshot; nothing else.  CPU: gcc 14 +
openmpi 5, Release, MPI.  GPU: apudev, 2 x gfx942, rocm 6.3, HSA_XNACK=1,
HSA_NO_SCRATCH_RECLAIM=1.

## Algorithm (src/rad_m1/rad_m1_vetcol.cpp, header comment has the details)

* Per radial column (i-line at fixed m, k, j), laterally homogeneous: the grey transfer
  equation with the column's own chi = rho (kappa_F + kappa_s) and the vet_sc source
  S = eps_th a T^4 + (1 - eps_th) E^n, eps_th = min(kappa_P/chi, 1).
* sp: **impact-parameter (p-ray) rays**, no mu interpolation.  Rays ordered by p:
  `vet_col_ncore` + 1 core rays (uniform in mu at the inner face), then per shell
  `vet_col_nsub` - 1 rays whose tangent point lies inside the interval below the shell and
  the ray tangent at the shell.  A shell sees the prefix of rays with p <= r (mu from 1 to
  0), so the beam of a free-streaming core is resolved at every radius.
* Sweeps: incoming top down from vacuum at the outer face; outgoing bottom up: core rays
  from the DIFFUSION intensity E + 3 F_r mu/c of the bottom cell at the inner face (as
  vet_sc); a ray tangent at the shell takes its own incoming intensity (mu = 0); a ray with
  its tangent inside the interval below goes from its incoming intensity down to the
  tangent point and back (the MIRRORED incoming intensity), chi and S interpolated in r
  there.  Segments: chi mean of the ends, S linear in tau (vet_sc's positive first-order
  weights), end half-cells with the face S extrapolated as vet_sc.
* Geometry in relative radii x = r/r_top (the column top, as the two-stream spherical code
  uses areas relative to the column top), segment lengths (x_b^2 - x_a^2)/(z_a + z_b)
  (no cancellation near a tangent point), times r_top once.
* Quadrature per shell: trapezoid in mu corrected by a factor (a + b mu + c mu^2) so
  that 1, mu, mu^2 are exact per hemisphere: f_K = 1/3 exactly for isotropic and
  diffusion intensities.
* Tensor: D = diag(f_K, (1-f_K)/2, (1-f_K)/2) about r_hat (`vet_col_axis = radial`,
  default) or about the M1 cell flux (`= flux`), f_K = K/J clamped to [1/3, 1].
* Cartesian: `vet_col_nmu` Gauss nodes per hemisphere, the plane-parallel 1-D version of
  vet_sc's sweep.
* When: in ImplicitSolve where vet_sc's formal solution runs (start-of-step E^n, T^n,
  before the predictor), every step or every `vet_col_every` steps.  The tensor goes to
  the solve through the tau closure's `tau_ten` (vet_col sets `tau_closure` as well), so
  rad_m1_implicit.cpp got ONE line (the call) and no Cartesian expression changed.
* GPU layout: one thread per column, rays looped inside, the running intensity per ray in
  `vcol_buf(ray, column)`; moments summed in a fixed order (deterministic, identical
  columns give identical tensors).  Column-local, no communication.
* Diagnostics: `vet_col_dump = <prefix>`, `vet_col_dump_every` (rank 0, block 0, column
  (ks, js): r, chi, S, J, H, K, f_K, chi used, E^n, F1).

Files: new `src/rad_m1/rad_m1_vetcol.cpp` (in src/CMakeLists.txt); `rad_m1.hpp` (members),
`rad_m1.cpp` (closure parse, options), `rad_m1_tau.cpp` (3 dispatch lines),
`rad_m1_sph.cpp` (vet_col allowed on the wedge), `rad_m1_implicit.cpp` (1 call).

## Gates

### T-bit (ref m1-sp2 vs new, bitwise = every output file)

* CPU, the 11 S1/S2 Cartesian + sp-hydro cases (`scripts/gate_tbit.sh`, `gate_tbit.log`):
  **11/11 BITWISE** (this time the explicit-transport restart layout also matched).
* S1/S2 sp gates (`scripts/gate_sp.sh`, `gate_sp.log`): T-S1 uniform and stretched, T-sym,
  T-sym gas, Marshak shell (Eddington), T-S2 free streaming, T-S4 atmosphere n = 32, T-S5
  radwave R = 100, T-sym M1: **9/9 BITWISE**.
* GPU, job 11958342 (`log.out.11958342`): He slab 1 rank, box3d_nd 2 ranks, box3d_nd +
  vet_sc 2 ranks: **3/3 BITWISE** (10/10 files each).

### T-S6 VET sphericity (`ts6.py` vs the independent reference `sphref.py`)

`sphref.py`: long characteristics per (r, mu), exact geometry (tau(s) = (C/p) atan(s/p)
for chi = C/r^2), Simpson on a grid uniform in tau, Gauss in mu split at the core edge;
checked against the analytic free-streaming J and K/J to 1e-15.  Code: one step of
`inp/sph_atm_vc` (r = 1..5, 4 x 4 wedge, analytic initial state, S = E = the spherical
Eddington profile, pure scattering), the first build dumped.  Runs `cpu/t6_*`, result
`ts6_result.txt`.  f_K error, L1 and Linf over the column:

| case | n_r (ncore 8, nsub 1) | L1 | Linf |
| --- | --- | --- | --- |
| (i) free streaming (kappa 1e-12; ref uses the code's bottom-cell E, F) | 32 / 64 / 128 / 256 | 6.9e-3 / 2.6e-3 / 7.9e-4 / 4.5e-4 | 5.4e-2 / 3.1e-2 / 1.9e-2 / 1.2e-2 (first shell) |
| (ii) extended atmosphere, chi = 10/r^2 (tau 8), ref f_K 0.337 (r=1) .. 0.389 (3) .. 0.545 (5) | 32 / 64 / 128 / 256 | 1.9e-3 / 7.7e-4 / 3.1e-4 / 1.2e-4 | 5.3e-3 / 2.6e-3 / 1.3e-3 / 6.9e-4 (top) |
| (iii) thick, chi = 1e4/r^2 | 32 / 64 / 128 / 256 | 3.4e-5 / 2.1e-5 / 1.3e-5 / 7.4e-6 | 2.9e-4 / 2.4e-4 / 2.6e-4 / 1.7e-4 (top cell) |

Top f_K at n = 256: (i) 0.98265 vs 0.98255, (ii) 0.53801 vs 0.53870.  (iii): the
interior |f_K - 1/3| ~ 1e-5; the max 2-4e-4 is the top cell.
Angles at n_r = 128: ncore 2/4/8/16/32 gives (i) L1 1.05e-2 / 2.0e-3 / 7.9e-4 / 1.3e-3 /
1.4e-3 (saturated by the first shell above the core, where the core edge is a
discontinuity), (ii) and (iii) insensitive (3.1e-4, 1.3e-5).  nsub 1/2/4: (i) Linf
1.9e-2 / 5.4e-3 / 2.0e-3, (ii) L1 3.1e-4 / 1.4e-4 / 1.4e-4.
Plane-parallel check (Cartesian): vet_col (nmu 4) vs vet_sc (nmu 4, nphi 8) on a
laterally uniform atmosphere, first tensor to 5e-11 (the precision of the vet_sc dump),
E after 20 steps to 1.3e-10 with implicit_predictor = none (`cpu/pp_*n`).

### T-S4 with vet_col (`ts4v.py`, `ts4rich.py`, result `ts4_result.txt`)

The runs_5b atmosphere (r = 1..5, rho kappa_t = 100/r^2, tau = 80, S2 settings, 800
steps, steady), closure = eddington / m1 (ref binary) / vet_col.  Reference: EXACT grey
spherical transfer by VET iteration (moment equations with f and g = H/J(r_out) from
`sphref.py`, S = J, converged to 5e-9 in 10 iterations): g_top = 0.597, f_top = 0.429.
The runs use `marshak_q = g_top` (with the code's q = 0.5 the outer BC error, identical
for all closures, hides the closure; `cpu/a4*` are those runs).  L1 of T_rad/T_ref - 1:

| n | Eddington | M1 | vet_col |
| --- | --- | --- | --- |
| 32 | 1.37e-2 | 5.8e-3 | 1.21e-2 |
| 64 | 8.8e-3 | 2.4e-3 | 6.7e-3 |
| 128 | 6.5e-3 | 2.7e-3 | 3.8e-3 |
| 256 | 5.3e-3 | 3.3e-3 | **2.1e-3** |
| Richardson 128/256 (end cells excluded) | 4.2e-3 | 4.3e-3 | **4.9e-4** |

Eddington and M1 converge to their own closure solutions (M1's error grows from n = 64 on),
vet_col converges to the exact transfer at first order (the Marshak end cell); with the
first-order term removed vet_col is 9x closer than M1 or Eddington.  Linf is the top cell
(same for all: set by the Marshak BC).  Luminosity L_f/L_in - 1 within 7.8e-6 (n = 256)
.. 5.6e-5 (n = 32), the same for all closures (the held gas's enthalpy flux, runs_5b).
nsub = 2 changes nothing (3.76e-3 at 128).

### Thin-cell stability (the runs_5b Finding 1 reproducer, `scripts/zt.sh`, `zt_result.txt`)

`sph_atm_bin`, kappa_s = 25 (rho kappa = 25/r^2), seed 1e-6, lin/Picard tol 1e-10.  Max
over i of the transverse spread of E (float32 bin dumps):

| run | cycle 10 | 20 | 50 | 400 |
| --- | --- | --- | --- | --- |
| Cartesian M1 | 6.1e-3 | 4.4 | 5.5 | |
| Cartesian Eddington | 0 | 0 | 0 | |
| Cartesian vet_col (radial / flux axis) | 0 | 0 | 0 | 0 / 0 |
| wedge M1 | 3.2e-5 | 1.1e-2 | 4.0 | |
| wedge vet_col (radial / flux axis) | 0 (50) | | 0 | 0 / 0 |

CONFIRMED: the fixed vet_col tensor is stable where M1 blows up (0 = below float32).

### Symmetry, restart, fast path (sp)

* T-sym (`sph_sym`, tau 2000, stretched r, 4 ranks, 500 steps, E changes by 0.80):
  spread(E)/E 5.2e-14 with lin tol 1e-14 and precond = line (`cpu/sym_vcl`); 3.8e-11
  with the default rbgs_fwd at the same tolerance (`sym_vcd`, the rbgs_fwd asymmetry of
  runs_5a), 1.6e-11 at default tolerances.  `sph_sym_gas` with c = 100: E 2.2e-14, rho
  5.6e-15, v_theta, v_phi / v_r <= 5e-16 (`cpu/symg_vc`).
* Restart (`sph_sym_gas_rst`, c = 100, 4 ranks, restart at cycle 250): all 24 later tab
  slices and both later restart files identical (`cpu/rstA_vc`, `rstB_vc`).  Holds for
  vet_col_every = 1; with every = k it holds only for a restart at a multiple of k.
* vet_col turns the fast default set on (a fixed tensor): fast (od_cache, krylov_fuse 3,
  op_stencil, rbgs_fwd) vs all off, 200 steps n = 64: E to 2.2e-9, F1 8.9e-11 (tolerance
  1e-10; `cpu/fast64`, `slow64`).  BiCGStab 3.7 vs 5.0 inner iterations per solve.

### Cost

* GPU (job 11958342, He wedge grid at half angular resolution: 96 x 128 x 128, stretched
  r, 16 blocks of 96 x 32 x 32 on 2 GPUs, 40 steps, same binary, arms interleaved, 2
  repeats; `inp/hewedge.athinput`): wall 1.844 / 1.848 s Eddington, **2.273 / 2.243 s
  vet_col (+22 %)**, 5.85 / 5.77 s M1.  vet_col build 11.9 / 11.8 ms per step = 0.47 s,
  21 % of the vet_col run; its solve is slightly cheaper than Eddington's (27.3 vs 29.5
  inner iterations per solve).  One thread per column (8192 per GPU) leaves the GPU
  mostly idle: the obvious speed-up is a team per column with the rays over the lanes.
* CPU share (correctness machine, 96 x 32 x 32, 1 rank, 20 steps, `cost_cpu.sh`, loaded
  node): build 103 ms of 11.9 s total = 17 %.

### Refusals

vet_col: more than one MeshBlock along x1 (also the implicit solve's own fatal), SMR/AMR,
periodic x1, cubed sphere, time_scheme = hesdirk2 (tau closure rule), transport !=
implicit.  Everything runs_5b still refuses on sp stays refused.

## HANDOVER

State: m1-vetcol = m1-sp2 (b5f3afea) + one commit (code, scripts, inputs, this README);
not merged (m1-sp2 is not merged either).  Build trees deleted; binaries kept in
/viper/ptmp2/jinma/vetcol_0924/bin.

Next steps:
1. GPU: team per column (rays over lanes, a fixed-order team reduction per shell) to take
   the build from 21 % to a few %; vet_col_every = k is already there.
2. Outer boundary: vet_col's f_K at the top is ~0.43-0.55 while the Marshak face uses q =
   0.5; feeding the formal solution's g = H/J at the top into the Marshak q would remove
   the dominant first-order top error seen in T-S4 (the exact g was 0.597).
3. S3 (cubed sphere) will reuse this file unchanged in the radial direction.
4. Default vet_col_nsub = 1: nsub = 2 only helps the artificial core-edge case (i).
