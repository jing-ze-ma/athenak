# tests_m1/gates: closing the M1 test gaps that let the vimp_fold double count through

- **Date:** 2026-09-24, viper. **Branch:** `m1-tests` (from rt-integration e89954e2),
  worktree `/viper/ptmp2/jinma/wt_m1tests`.
- **Run tree:** `/viper/ptmp2/jinma/m1tests_0924` (`base/` = e89954e2, `new/` = this branch,
  `old/` = e0074273 + the check; binaries in `bin/`, logs in `matrix_*`, `gates_*`,
  `bitwise/`, tables in `matrix_*.txt`, `gates_*.txt`, `bitwise_off.txt`).
- **The bug** (tests_m1/runs_4l_sync sect. 1, fixed in 89174004): under
  `implicit_vimp_fold`, `ImplicitStencilOpPart` (the overlapped operator of
  `implicit_halo_overlap`) added the unfolded `M1VimpRow` to the folded stencil. The
  x2/x3 +-1 vimp terms were counted twice.
- **Why it slipped through:**
  1. the overlap on/off gate of runs_4j accepted solver-level differences (KE <= 6e-8);
     the bug was 5.2e-8;
  2. the outer iteration hides an inexact inner operator;
  3. vimp_fold and overlap were gated separately;
  4. there was no direct operator test;
  5. CI runs no M1 test.

This directory closes (1)-(4) and adds the CI tests for (5).

## 1. The operator-equivalence check: `<rad_m1>/implicit_op_check`

`src/rad_m1/rad_m1_opcheck.cpp`, `RadiationM1::ImplicitOpCheck`. Hooks:
- a call at the top of `ImplicitBiCGStab`;
- 9 lines in `ImplicitInit` that read the switch;
- 7 lines in `rad_m1.hpp`;
- one line in `src/CMakeLists.txt`.

**Switch.** `implicit_op_check = K` is read only when named. The default 0 is off, and
then nothing runs. The check covers the first |K| solves.
- K > 0: fatal on a mismatch.
- K < 0: report only.
- `implicit_op_check_tol` (default 1e-12) is the pass threshold.
- AthenaK rejects a command-line key that the input file does not name. An input that
  wants the switch carries `implicit_op_check = 0` (see the inputs here and
  `tst/test_suite/rad_m1/he_slab_m1.athinput`).

**What one check does, before the solve touches its vectors:**
- x = a fixed hash of (gid, k, j, i) in the active cells, and **1e20 in every ghost**.
  (s, rhat) = two more hashes.
- Every halo path the mesh allows fills the ghosts of x:
  - the ordinary exchange (`pbval_kr`);
  - `implicit_halo_direct`;
  - `implicit_halo_mpi`.

  **Each ghost within the operator's reach that has a neighbour block must equal that
  block's active cell exactly.** 576 halo checks and 13.8 M ghosts, 0 wrong.
- Every operator variant the configuration allows is applied, whatever the input
  selects:
  - `stencil`: `ImplicitStencilOp`;
  - `stencil_split_red`;
  - `overlap`: `ImplicitStencilOpPart` interior + shell;
  - `overlap_faces`;
  - `krylov_dev`: `ImplicitOpXD`, halo read in place with the ghosts still poisoned;
  - `krylov_dev_halo`;
  - `legacy_od_cache`;
  - `legacy_7pt`: the 7-point row + `M1VimpRow` + `ImplicitOffDiagOp`, the operator
    before the stored stencil.
- Under implicit_vimp, the stencil is rebuilt with `implicit_vimp_fold` flipped on a
  separate array. Stencil, overlap and device variants are run again (`+fold` / `+nofold`).
- Each y is compared with the first one (the production halo + `stencil`) **per row**:
  max |y - y_ref| / sum_o |a_o x_o|. Round-off is a few 1e-16 in this measure. Against
  max |y|, the old bug was only 4.5e-14, so a global norm would hide it.
- The red = 3 reductions ((y,s), (y,y), (rhat,y)) are compared as well.
- A reference |y| > 1e12 means an unfilled (poisoned) ghost is read with a non-zero
  coefficient.
- Everything is restored afterwards: the whole `iw`, the stencil, `kd` and the flags.

**Combination matrix** (`opcheck_matrix.sh`; new binary; one cycle, K = -2, i.e. the first
two solves):
- Settings covered:
  - vimp on/off;
  - vimp_fold on/off;
  - halo_mpi on/off;
  - overlap on/off, ovl_faces on/off (2 and 4 ranks);
  - krylov_dev on/off (1 rank);
  - 1, 2 and 4 ranks.
- Geometries: 3-D box 84x32x32 (4 blocks) and 2-D slab 84x32, each with closure Eddington
  and with vet_sc + full tensor (non-zero edge coefficients).

| geometry | runs | FAIL | worst dy/row |
|---|---|---|---|
| box, Eddington | 36 | 0 | 6.2e-16 |
| slab, Eddington | 36 | 0 | 4.2e-16 |
| box, vet_sc | 36 | 0 | 8.5e-16 |
| slab, vet_sc | 36 | 0 | 8.5e-16 |

- Exact (0) wherever the arithmetic is the same.
- Every run had 0 unfilled-ghost reads and 0 wrong ghosts.

**Proof on the old bug.** `old/` = e0074273 (before the fix) with this check ported. The port
is the same file minus the lines marked `// OVLF`, since `implicit_halo_ovl_faces` does not
exist there. The same matrix, box + slab, Eddington (`matrix_old.txt`):
- **32 of 72 runs FAIL.** They are exactly the implicit_vimp runs on 2 and 4 ranks, under
  every fold / halo_mpi / overlap setting of the input, because the check forces the
  overlap + fold variant.
- The failing variants are `mpi/overlap` (the input folds) and `mpi+fold/overlap` (the check
  folds), 32 each.
  - dy/row is 1.0e-10 at the first solve and up to 5.3e-8 at the second
    (`matrix_old/box_n2_v1f1_h1o1a0_k0/log.txt`).
  - Every other variant stays <= 6.2e-16 (the unfolded overlap included).
- Runs without vimp, and 1-rank runs, pass: the bug needs both fold and overlap.

## 2. Comparison gates: `gates.py` + `cmp.py`

`python3 gates.py run <athena> <rundir>; python3 gates.py eval <rundir>`.
- Inputs: `box3d_be_x.athinput`, `slab2d_plm_vimp_be_x.athinput` here, with a history row
  every cycle.
- Config: the be levers (vimp_fold, fast_kernels, one_pass 4, predictor order 2).
- Pairs:
  - overlap on/off, ovl_faces on/off, krylov_dev on/off, halo_mpi on/off;
  - 1 vs 2, 1 vs 4, 2 vs 4 ranks.
- Geometries: the box (6 cycles) and the vet_sc slab (12 cycles).
- Each arm runs at the input tolerances (loose: implicit_tol 1e-8, lin_tol 1e-10) and at
  TIGHT = implicit_tol 1e-12, lin_tol 1e-11.
  - lin_tol 1e-12 stalls at round-off: 100 breakdowns and 20 line-Jacobi fallbacks in 30
    box cycles.
  - 1e-11 converges cleanly.

**What the measurements forced on the design** (`runs/tt_*`, `runs/e_*`, and the 30-cycle
arms of the first attempt):
- **Tightening does not shrink the state differences of a 30-cycle run.**
  - Per MeshBlock state after 30 box cycles, overlap on/off: F_rad 2e-10, momentum 8e-8
    of max|rho v|, at both tolerances.
  - The history KE1/KE2 are ~1e-10 loose and ~1e-9 tight.
  - The velocity is the small residual of large balanced forces. Differences of ~1e-11
    in E_rad (the floor that lin_tol and the gas coupling allow) grow to 1e-8 of the
    velocity within tens of cycles. So the gates use short runs and split the measures.
- **The first cycle shrinks.** The new overlap pair goes 9.8e-12 -> 5.2e-13; the old bug stays
  at 2.5e-9 -> 2.5e-9. This is the solver-noise test.
- The bin outputs are single precision (a 1.1e-7 floor), so `cmp.py` reads the double
  state from the restart file (layout in its header).

**Criteria** (cmp.py measures; tight arm unless stated):
1. hst_cons <= 1e-10: mass, tot-E, Etot, F1top/mid/bot, Fres, every row.
2. row1_cons shrinks: tight <= loose/10, or <= 1e-12.
3. hst_dyn <= 3e-8: KE, net momenta (scaled by sqrt(2 M KE)), V1max. This is a coarse
   net: V1max moves 1.5e-8 between 1 and 4 ranks in 6 box cycles.
4. rst_cons <= 1e-9: dens, gas E, E_rad, F_rad (vector) in the active cells.
5. No NON-CONVERGED step. rst_dyn (the momentum field) is reported only; it is 4-15e-8 in
   every pair.

Results (`gates_new.txt`, `gates_old.txt`):

| pair (new) | box: hst_cons, row1 loose->tight, hst_dyn | slab | verdict |
|---|---|---|---|
| overlap on/off | 2.8e-11, 9.8e-12->5.2e-13, 3.4e-9 | 3.2e-11, 2.6e-12->7.4e-13, 9.6e-10 | PASS |
| ovl_faces on/off | 2.4e-11, 7.8e-12->5.5e-13, 8.6e-9 | 2.1e-11, 9.8e-12->7.2e-13, 7.5e-10 | PASS |
| krylov_dev on/off | 0 (bitwise) | 0 | PASS |
| halo_mpi on/off | 0 (bitwise) | 0 | PASS |
| 1 vs 2 ranks | 1.3e-11, 2.4e-12->4.9e-13, 8.3e-9 | 3.3e-11, 2.3e-12->2.6e-13, 1.2e-9 | PASS |
| 1 vs 4 ranks | 2.8e-11, 3.1e-13->5.9e-13, 1.5e-8 | 3.3e-11, 1.9e-12->4.1e-13, 1.3e-9 | PASS |
| 2 vs 4 ranks | 2.8e-11, 2.8e-12->4.7e-13, 1.4e-8 | 5.0e-11, 4.2e-12->3.9e-13, 2.3e-9 | PASS |

**Old binary (e0074273):** the overlap pair **FAILS** on criteria 1 and 2.
- Box: hst_cons 2.5e-9 (Fres), with no shrink (2.5e-9 -> 2.5e-9) and hst_dyn 8.0e-8.
- Slab: 1.6e-9 with no shrink (1.6e-9 -> 1.6e-9).
- All other pairs pass as for the new binary. ovl_faces does not exist there, so f2 = o2.

## 3. CI tests: `tst/test_suite/rad_m1/`

- `m1_common.py` builds a `PROBLEM=box_convection` binary once per session and build type,
  in `tst/build_m1` and `tst/build_m1_mpi`. `conftest.py` removes the builds at the end of
  the session.
- The four data files that are not in git (the He Rosseland and Planck tables, two
  bench/m1_stage2 profiles) are looked for:
  - in `$ATHENAK_M1_DATA`;
  - then at `data/stellar_opac` and `../bench`.

  Without them the tests SKIP, like the dhj correlated-k tests.
- `he_slab_m1.athinput` is the slab input with the data paths set by the test.

| test | what | run time, 2 sessions (excl. build) |
|---|---|---|
| `test_rad_m1_slab_cpu.py` | serial, one block, 12 cycles Eddington + 12 vet_sc at TIGHT; last hst row vs stored REF, rtol 1e-9 (radiation, totals) / 1e-7 (KE, V1max) | 17 s, 27 s |
| `test_rad_m1_opcheck_mpicpu.py` | implicit_op_check = 2 (fatal), 2 ranks vet_sc + 4 ranks Eddington, halo_mpi + overlap + faces + fold; asserts PASS and that the overlap/fold variants ran | 17 s, 24 s |
| `test_rad_m1_restart_mpicpu.py` | 2 ranks, overlap + faces + be levers: 16 cycles vs 8 + restart + 8, rst data and shared hst rows (dt column excepted) bitwise | 25 s, 34 s |

- The box_convection build adds ~3.5 min per build type on the login node, once per
  session. `ATHENAK_M1_KEEP_BUILD=1` keeps the builds for the next session.
- The REF numbers are from e89954e2 + this branch, gcc 14, serial.

## 4. Switch off = bitwise (`bitwise_off.sh`, `bitwise_off.txt`)

The comparison is base e89954e2 against this branch, with the key named = 0 or absent
(`slab2d_nd`).

**Result: hst and the final rst after `<par_end>` are BITWISE in all 7 cases:**
- box 1 rank;
- box 2 ranks, overlap + faces + fold;
- box 4 ranks, overlap + faces + fold;
- slab 1 rank;
- slab 2 ranks, overlap + faces;
- slab 2 ranks, vet_sc;
- slab2d_nd 2 ranks, tlim 200.

## Files

- `opcheck_matrix.sh`: check 1, the matrix.
- `gates.py`, `cmp.py`: check 2.
- `bitwise_off.sh`: the switch-off gate.
- `box3d_be_x.athinput`, `slab2d_plm_vimp_be_x.athinput`: the inputs.

The data paths in these inputs are the viper ones.
