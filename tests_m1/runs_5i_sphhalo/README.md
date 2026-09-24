# runs_5i_sphhalo: the fast implicit-M1 halo on the spherical-polar wedge (m1-sphhalo)

Branch m1-sphhalo from rt-integration d9583121.  Run trees, binaries and logs are in
/viper/ptmp2/jinma/sphhalo_0924 (`runs/`, `gpu/`, `gates/`, `*.log`).
ref = `git archive d9583121`, new = this branch.  The scripts and inputs are copied under
`scripts/` and `inp/`.
- CPU: gcc 14 + openmpi 5, Release, MPI.
- GPU: apudev, MI300A gfx942 APU, rocm 6.3, `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1`,
  job 11965648.

## Change

| file | change |
| --- | --- |
| `src/rad_m1/rad_m1_implicit.cpp` | `hmdef` (the implicit_halo_mpi default) no longer excludes `use_spherical_polar`.  On sp it is `false` only for a restart (`restart_run`) whose file lacks the key.  The poles (`use_polar_boundary`), the cubed sphere and SMR/AMR stay excluded.  implicit_halo_overlap and implicit_halo_ovl_faces follow as on the Cartesian mesh: on with halo_mpi and more than one rank, and off for a restart whose file lacks the key. |
| `src/rad_m1/rad_m1_sph.cpp` | `SphericalS1Check` no longer refuses `implicit_halo_mpi = true`, and its message and header are updated. |

**Why nothing else changes.**
- **A wedge has no pole**, so every M1 scratch-halo ghost is a plain same-level copy.
- **The component-role table does not apply.** Design sect. 2 lists the vector pairs of the
  sp exchanges: N2/N3, A2/A3, V2/V3 and DV2/DV3. They act in `bvals_cc.cpp` only under
  `do_pole` (the sign flip) and on cubed-sphere seams.
- **The Krylov vector is a scalar.** `ImplicitHaloMPIInit` already refuses the polar
  boundary, cubed sphere and multilevel meshes.
- **The operators read stored coefficients.** The stencil (`ImplicitStencilOp`) and its
  interior/shell split (`ImplicitStencilOpPart`) read coefficients that the S1/S2 sp
  blocks overwrite when they are built. So the split sees the sp geometry with no sp code.
- **Still refused on sp.** Poles are fatal in `rad_m1.cpp`, SMR is fatal in
  `SphericalS1Check`, and cs is fatal in the constructor.

**New defaults on the sp wedge:**
- implicit_halo_mpi, implicit_halo_overlap and implicit_halo_ovl_faces are on for 2 or more
  ranks wherever implicit_halo_direct is on. implicit_halo_direct is on for eddington and
  vet_col with bicgstab.
- m1/minerbo/kershaw keep halo_direct off by default, so halo_mpi stays off for them.
- Every sp restart written since S1 echoes `implicit_halo_mpi = false` and keeps it.

## Gates

### implicit_op_check on sp (`scripts/opchk.sh`, runs `runs/chk_*`, K = -2, 2 solves each)

- **Wedge:** hewedge, 32 x 12 x 12, r stretched, theta pi/4..3pi/4 and phi 0..pi/2.
- **Geometries** (MeshBlocks, ranks):
  - t2 = 2 x theta on 2 ranks;
  - p2 = 2 x phi on 2 ranks;
  - q2 = 2 x 2 on 2 ranks;
  - q4 = 2 x 2 on 4 ranks;
  - r4 = q4 with theta reflecting instead of periodic;
  - t4 = 4 x theta on 4 ranks (nx2 = 16);
  - f4 = 4 x phi on 4 ranks (nx3 = 16).
- **Closures:**
  - eddington;
  - m1, with implicit_halo_direct and implicit_op_stencil named true (they default off
    for m1);
  - vet_col.
- **Paths checked** (the reference is the production mpi/stencil):
  - halo: mpi and exch;
  - variants: stencil_split_red, overlap, overlap_faces, legacy_od_cache (edd, vc) and
    legacy_7pt, on each halo.

| closure | runs | variant checks | worst dy/row | worst reduction diff | ghosts checked | wrong | FAIL |
|---|---|---|---|---|---|---|---|
| eddington | 7 | 126 | 0 | 2.8e-15 | 96256 | 0 | 0 |
| m1 | 7 | 98 | 0 | 2.9e-15 | 96256 | 0 | 0 |
| vet_col | 7 | 126 | 0 | 2.9e-15 | 96256 | 0 | 0 |

- The split operator (overlap and overlap_faces) equals ImplicitStencilOp exactly, row by row.
- The only differences are in the red = 3 dot products, which are summed in two parts.
- Without an input key, eddington and vet_col select mpi/stencil (the new default).

### halo_mpi on vs off on sp (`runs/w_*`, 4 ranks, 2 x 2 blocks, 10 cycles; `scripts/sp_pairs.py`)

- on = the new default: halo_mpi + overlap + faces.
- off = implicit_halo_mpi = false.
- L = the input tolerances (implicit_tol 1e-10, lin 1e-10). T = tight (1e-12, 1e-11).
- The table gives the final restart state, max relative difference.

| closure | tol | dens | gas E | E_rad | F_rad | NON-CONVERGED |
|---|---|---|---|---|---|---|
| vet_col | L | 0 | 3.1e-12 | 8.1e-14 | 2.7e-13 | 0 |
| vet_col | T | 0 | 3.6e-12 | 7.2e-14 | 2.5e-13 | 0 |
| eddington | L | 0 | 3.2e-12 | 1.0e-14 | 2.2e-14 | 0 |
| eddington | T | 0 | 3.5e-12 | 1.4e-14 | 2.7e-14 | 0 |
| m1 | L, T | bitwise | | | | 0 |

- **Round-off reference:** the same off input on 2 vs 4 ranks, both with the ordinary
  exchange:
  - vet_col L: gas E 3.4e-12, E_rad 7.4e-14, F 2.5e-13;
  - eddington T: gas E 3.8e-12, E_rad 9.6e-15, F 2.1e-14.
- So on vs off sits at the rank-count round-off floor. That is why tightening the tolerance
  does not shrink it.
- The history file prints 6 digits and cannot resolve these levels.

### Restart on sp (vet_col, 4 ranks, all three switches on)

- **rsA:** 20 cycles with a restart every 10.
- **rsB:** restarted from rsA's cycle-10 file and run to cycle 20. Both later restart
  files have identical payloads (bitwise).
- **rsN:** the same restart file with the three implicit_halo_* lines removed from its
  header (an "old" restart).
  - It runs with the ordinary exchange: no "halo_mpi: ON", and it echoes `implicit_halo_mpi
    = 0`, `overlap = 0`, `ovl_faces = 0`.
  - Its final state is identical to rsA's.

### T-S1, T-S4, T-sym, ref vs new (`scripts/gate_sp.sh`, `tabcmp.py`, `cmpdir.py`)

| test | ranks | new default | explicit halo_mpi = false |
|---|---|---|---|
| T-S1 sph_diff n = 64, uniform r and he4 stretch | 1 | bitwise (3/3 files each; one block, halo not in play) | |
| T-S4 sph_atm (m1), n = 64, 4x4 | 1 | bitwise 16/16 | |
| T-S4 sph_atm (m1), 8x8 angular, 4 blocks | 4 | bitwise 16/16 (m1: halo_mpi stays off) | bitwise 16/16 |
| T-S4 sph_atm_vc (vet_col), n = 64, 4x4 | 1 | bitwise 16/16 | |
| T-S4 sph_atm_vc (vet_col), 8x8, 4 blocks | 4 | halo_mpi ON: dE/E 3.3e-11, dF/max abs(F) 1.3e-10 (m1 tabs, 800 cycles, lin_tol 1e-10); transverse F/F_r 1.6e-11 (ref 4.6e-12) | bitwise 16/16 |
| T-sym sph_sym (Eddington), 500 cycles | 4 | halo_mpi ON: dE/E 7.3e-16, dF 8.3e-15; spread(E)/E 1.638e-11 (ref 1.638e-11), F2/F1 2.867e-11 (ref 2.867e-11) | bitwise 18/18 |

Every sp input that names implicit_halo_mpi = false is bitwise against ref.

### Cartesian: new vs ref

- **CPU set** (`scripts/gate_cart.sh`, the runs_5a T-bit set without the sp-hydro case):
  - He slab 1 and 2 ranks, slab new defaults, vet_sc, vimp, box3d_old and box3d_nd
    (halo_mpi on, 4 ranks), explicit, Marshak implicit_x1, Marshak 3-D.
  - All 10 are bitwise.
  - The explicit case (expl2) first showed its two restart files differing by 88 inserted
    zero bytes (new). A rerun of the same new binary was bitwise against ref, and ref vs a
    ref rerun was bitwise too. This is a one-off restart-write glitch under a loaded login
    node, not the code.
- **gates.py** (tests_m1/gates, new binary, /viper/ptmp2/jinma/sphhalo_0924/gates.log):
  GATES: PASS.
- **tst** (in the new snapshot, ATHENAK_M1_DATA = /viper/ptmp2/jinma/faces_0924/m1data):
  test_rad_m1_slab_cpu, test_rad_m1_opcheck_mpicpu and test_rad_m1_restart_mpicpu all
  passed.
- **GPU** (job 11965648): ref vs new is bitwise for all three cases, 9/9 files each:
  - slab 1 GPU;
  - box3d_nd 2 GPUs;
  - box3d_nd + vet_sc 2 GPUs.

## GPU timing: sp He wedge, vet_col, 2 GPUs (job 11965648)

- **Setup:**
  - input `inp/hewedge.athinput`: 96 x 128 x 128, 16 blocks of 96 x 32 x 32, the current
    vet_col defaults (team build);
  - off = `inp/hewedge_off.athinput`, which names implicit_halo_mpi = false;
  - same binary (athena_new_none_gpu, md5 4187aa05), interleaved, 3 repeats each;
  - elapsed time over cycles 10-30 of 40.
- **Profiles:** rocprofv3 traces, `prof.py` window cycles 10-40, rank 0. The files are
  `p_off.r0.sum` and `p_on.r0.sum`.

| arm | ms/cycle (3 repeats) | profiled total | KRY:halo (kernel + idle) | matvec+dot | vec+red | precond |
|---|---|---|---|---|---|---|
| off (ordinary bvals exchange) | 33.48, 33.59, 33.43 | 33.22 | **12.26** (5.77 + 6.49) | 3.05 | 1.80 | 4.57 |
| on (halo_mpi + overlap + faces, new default) | 22.24, 22.34, 22.29 | 22.82 | **1.72** (1.13 + 0.59) | 4.03 | 1.33 | 4.24 |

- The halo drops from 12.3 to 1.7 ms/cycle, and the cycle from 33.5 to 22.3 ms (-33 %, 1.50x).
- The overlap split adds 1.0 ms/cycle of matvec kernel time (StencilOpPart interior + shell).
- VET build (3.4 ms) and hydro (5.2 ms) are unchanged.

## Not done / notes

- The poles stay refused, so the polar ghost flips of the role table are untested here.
- m1/minerbo/kershaw on sp get halo_mpi only when implicit_halo_direct is named true, the
  same rule as on the Cartesian mesh (fdef needs a fixed closure).
