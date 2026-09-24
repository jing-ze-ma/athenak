# runs_5e_vetcol2: vet_col surface Marshak q from the formal solution + team-parallel build

Branch m1-vetcol2 from rt-integration 9c1a12d7 (S5 vet_col merged).  Work dir, binaries,
run trees: /viper/ptmp2/jinma/vetcol2_0924 (`bin/`, `cpu/`, `gpu/`).  ref = snapshot of
9c1a12d7 (`git archive`), new = snapshot of this branch (`scripts/snap_new.sh`,
`build.sh`).  Final new binaries: none cpu 02c6d5fb, none gpu c8a0f20b, box_convection gpu
324ce55f; box_convection cpu ec29f9da and sp_test cpu 21aab742 are from the snapshot
before the last team-only edit (chunk default and the knobs below; the Cartesian and box
paths did not change).  CPU: gcc 14 + openmpi 5, Release, MPI.  GPU: apudev, 2 x gfx942,
rocm 6.3, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1.

## A. `<rad_m1>/vet_col_surface_q` (default false)

Under closure = vet_col, the outer-x1 Marshak face flux F = c q E(top cell) takes q per
column from the formal solution:

    q = H(top face) / J(top cell),   clamped to [vet_col_surface_qmin, _qmax] = [1e-3, 1]

The outgoing intensities of the top shell are carried over the last half segment to the
face (the incoming sweep's first segment backwards).  H(face) = 1/2 sum wf mu_f I_f: sp:
every ray at mu_f = sqrt(1 - p^2) plus a node at mu = 0 with zero intensity, the same
moment-corrected trapezoid as the shells; Cartesian: the Gauss nodes.  The ratio uses J
at the top CELL, not at the face.  With E(top cell) = J, the discrete face flux is then
the formal solution's own H, so the BC carries no O(dx) mismatch.  q is built in the same
call as the tensor (start-of-step state, lagged the same way), stored in `vcol_q(m,k,j)`,
and read only at the OUTER x1 face.  In rad_m1_implicit.cpp, `mq` is shadowed in the 4
outer-Marshak branches (m1_impl_asm x2, m1_impl_face, m1_vimp_p) by
`const Real mq = vqs ? vq_(m,k,j) : mqo;`.  The Cartesian expressions are textually
unchanged.

### T-S4 (the S2/S5 grey spherical atmosphere, tau 80, `scripts/ts4sq.sh`, `ts4_result.txt`)

L1 of T_rad/T_exact - 1 against the exact-transfer reference of runs_5d (`ts4v.py`,
`ts4rich.py`):

| n | S5 vet_col, q = exact g_top 0.597 | **vet_col + surface_q** | Linf (surface_q) |
| --- | --- | --- | --- |
| 32 | 1.21e-2 | **4.8e-4** | 2.4e-3 |
| 64 | 6.7e-3 | **2.5e-4** | 1.8e-3 |
| 128 | 3.8e-3 | **1.18e-4** | 1.0e-3 |
| 256 | 2.1e-3 | **5.2e-5** | 4.9e-4 |
| Richardson 128/256 | 4.9e-4 | **2.3e-5** | 2.4e-4 |

First order, with an error 25-40x smaller than S5 with the exact g_top (which is the
wrong number for a cell-centred E: q_eff = F/(c E_top) settles at 0.38 / 0.45 / 0.51 / 0.54
for n = 32..256).  With q = 0.5 the top error dominated every closure (runs_5d).
L_out/L_in - 1 is unchanged (7.8e-6 at n = 256).  Richardson 32/64, 64/128: 1.2e-4,
5.0e-5.

### Plane-parallel Milne on Cartesian (`inp/milne_vc.athinput`, `milne_v.py`, `milne_result.txt`)

The same atmosphere as a Cartesian slab: rho kappa_t = 100, x = 0..0.8 (tau 80), flux
bottom, Marshak top, vet_col nmu 4 (and 8).  The reference is the exact Milne solution
J = 3H(tau + q(tau)), with the Hopf function from a VET iteration with the exact E_n
kernels in `milne_v.py`: q(0) = 0.577350, q(0.01) 0.588235, q(0.1) 0.627919, q(inf)
0.710446, which reproduces the gate3_milne.txt table of runs_3j_vet.  L1 (Linf) of T:

| n | q = 0.5 (nmu 4) | surface_q (nmu 4) | surface_q (nmu 8) |
| --- | --- | --- | --- |
| 32 | 2.0e-2 (0.24) | 5.4e-4 (5.1e-3) | 5.4e-4 (5.1e-3) |
| 64 | 9.7e-3 (0.16) | 2.3e-4 (3.2e-3) | 2.4e-4 (3.1e-3) |
| 128 | 4.2e-3 (0.091) | 7.5e-5 (1.0e-3) | 7.4e-5 (1.1e-3) |
| 256 | 1.8e-3 (0.044) | 1.1e-4 (5.3e-4) | 1.1e-4 (2.8e-4) |

The top-cell E error drops from -16 % (q = 0.5, n = 256) to +0.2 % (nmu 4) or -0.03 %
(nmu 8).  L1 floors at ~1e-4 from n = 128 on, with nmu 4 and 8 alike.  That floor is not
the surface BC (Linf keeps falling); it is the deep interior, likely the reference's
spline or the end cell.  Improved 17-40x in L1, never worse.

### Gates (A)

* Off: bitwise.  CPU Cartesian 11/11 (`gate_tbit.log`); sp S1/S2/S5 9 + 6 vet_col cases
  15/15 (`gate_sp.log`); GPU 3/3 (`log.out.11959508`).
* Restart with surface_q on (`sph_sym_gas_rst`, c = 100, 4 ranks, restart at 250, rebuild
  every step): 26/26 files bitwise (`rst_sq.log`).
* tst/test_suite/rad_m1: slab_cpu, opcheck_mpicpu, restart_mpicpu all pass
  (ATHENAK_M1_DATA = /viper/ptmp2/jinma/faces_0924/m1data).
* Cost: +0.17 ms per build (3.48 vs 3.31 ms, GPU).

**Recommendation: make vet_col_surface_q = true the vet_col default.**  It removes the
dominant top error in both geometries, at no measurable cost.  This commit keeps it
default off (bitwise to S5) for the merge decision.

## B. Team-parallel tensor build (`vet_col_team`, default TRUE)

`VetColBuildTeam`: one team per column (Kokkos TeamPolicy, AUTO team size).  Team scratch
holds the column profile (chi, S, the incoming J/H/K sums), the running intensity per ray,
and the intensities of a CHUNK of lc shells.  Per shell, the rays are spread over the
threads (each ray's recurrence stays in one thread per shell step, with the same
arithmetic).  After each chunk, one thread per shell sums the moments in the ray order
r = 0, 1, ... of the one-thread kernel.  Every operation and its order are unchanged, so
the result is BITWISE:
* CPU: team vs one-thread, sp T-S4 n = 32/64/128 and Cartesian Milne, every file identical.
  The 6 vet_col sp gates (ref = one-thread kernel of 9c1a12d7, new = team) are bitwise.
* GPU: T-S4 n = 64 with surface_q, team vs one-thread: 8/8 files bitwise
  (`log.out.11959508`).

Knobs: `vet_col_team_size` (0 = AUTO) and `vet_col_chunk` (0 = min(16, 24 kB scratch
budget)).  Sweep, job 11959410 (`log.out.11959410`), ms per build: AUTO/64/128 threads
x lc 8/16/32 = 3.6-3.8 / 3.3-3.7 / 7.6-8.6.  At lc = 32 the 32 kB scratch halves the
occupancy.  AUTO with lc = 16 is best.  `vet_col_team = false` keeps the old kernel.

### GPU cost (job 11959508, He wedge grid 96 x 128 x 128, 16 blocks on 2 GPUs, 40 steps, same binary, interleaved, 2 repeats)

| arm | wall (s) | build ms/step | build share |
| --- | --- | --- | --- |
| Eddington | 1.842 / 1.837 | | |
| vet_col, one-thread (S5 kernel) | 2.268 / 2.221 | 11.86 / 11.85 | 21 % |
| **vet_col, team (default)** | **1.881 / 1.925** | **3.31 / 3.36** | **7 %** |
| vet_col, team + surface_q | 1.874 / 1.908 | 3.48 / 3.48 | 7 % |

The build is 3.6x faster, and vet_col now costs +2 to +5 % over Eddington (was +22 %).
Its solve is slightly cheaper (27.3 vs 29.5 inner iterations per solve).

## Files

`src/rad_m1/rad_m1_vetcol.cpp` (VcolQuad factored out, with bitwise-identical shell
weights; face tables `vcol_muf`/`vcol_wf`; q in both kernels; `VetColBuildTeam`),
`rad_m1.hpp` (members), `rad_m1.cpp` (options), `rad_m1_implicit.cpp` (the q shadow
lines only).
Scripts: `scripts/` (gate_tbit, gate_sp, gate_gpu, tune_gpu, ts4sq, milne8, rst_sq,
build, snap_new, run).  Inputs: `inp/` (the S5 inputs with the new keys; milne_vc).

## Next

1. Decide vet_col_surface_q default on (recommended).
2. The remaining 3.3 ms: the per-shell team barriers (192 per build) and the serial
   per-shell sums.  A global (column, shell, ray) intensity buffer would allow full-width
   sums at ~1 GB/GPU memory.  Not needed at 7 %.
