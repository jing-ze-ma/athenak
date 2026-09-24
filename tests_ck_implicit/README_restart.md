# ck-restart: bitwise restarts for the implicit correlated-k RT (T4, c2, cadence), 2026-09-24

Branch `ck-restart` (worktree /viper/ptmp2/jinma/wt_ckrst, from rt-integration e89954e2).
Runs and binaries: /viper/ptmp2/jinma/ckrst_0924. Scripts: `tests_ck_implicit/restart/`.

**Before.** A restart of an implicit-ck run was not a bitwise continuation: c2 alone gave
a max e difference of 3.1e-5 after 4 cycles (cadence/rstgate3_cmp.txt), and the cadence
forced a full call after every restart. T4 was not bitwise either (ck_thk, below).

**After.** N cycles straight = N/2 + restart + N/2, bitwise (rst payload, .hst, bin
dumps), for T4, c2, c2 + every = 4 and c2 + every = 4 + guard, on CPU; and for c2 +
every = 4 on the GPU on the production state (section 4).

## 1. Inventory: what persists from one RT call to the next

| state | who reads it before writing it | treatment |
|---|---|---|
| `ck_thk` (thin/thick verdict of the last apply) | rt_pre_opac of pass 0 (T4 and everything built on it) | **written**: slab THK |
| cadence `ck_cad` slots 0-3: Q0, D, T0, rho0 | every linearised step and the guard | **written**: slabs CAD0-3. Slots 4-5 are per-step scratch. The mask is rewritten by every linear step. |
| `ck_cv` (secant cv the call ended on) | pass 0, only with `ck_impl_cvkeep` | **written** when cvkeep is on: slab CV |
| ck_impl_xstep **stored operator**: kappa rho, coefficient triple, tm factorisation (lP, lG, lC), tpf, the cut, frozen opacity/Planck caches, the frozen beam Qb | a call that re-applies it (ncycle - xs_cyc < xstep) | **rebuilt**. It is O(100) Reals per cell (13 x 88 chains), too large to write. It is a deterministic function of what pass 0 of the storing call read. The store now keeps a device copy of that input: u0, wtemp, bcc0 (MHD) and ck_thk at the call's start, plus the call's bdt (`CkXsSnap`). When the next call can re-apply, the restart writes that input (slabs XS*). The first call after the restart re-runs that pass 0 on it (`CkRstRebuild`), then puts the restarted state and the end-of-call ck_thk back. |
| `ck_impl_xs_cyc` | the xstep reuse decision | **written** per rank (slab XSCYC). It is rank-local under xstep_thr and the guard. |
| xstep store bdt | the rebuild | **written** per rank: slab XSBDT |
| `ck_xsT/ck_xsD` (xstep_thr) | the thr test | rebuilt by the rebuild's store |
| `ck_impl_jac_built` | pass 0 with `jreuse_xc` | **written**: header |
| counters nstore, nreuse, nsweep, cadence nfull/nlin/nguard/fref | verbose lines only | **written**: header, so the logs continue |
| ck_rprev (pred), ck_ep/ck_tp (cvsec), ck_jac, ck_done, ck_dep/ck_seed, prev_res, jac_again, per-call counts | reset or rewritten before being read in every call (pz_, cv0_, CkWarmSeed, the call prologue) | nothing to do |

**Not covered** (documented, not fixed): `ck_impl_warm` (carries ck_dep and ck_ei; it is a
failed lever and refused with the cadence), `ck_impl_reuse_jac = 2`, `ck_impl_glob`,
`ck_impl_aa`.

**A bug found on the way** (it only shows with the restart). When the first implicit call
of a run is a cadence **guard** call, its `ck_done = mask` copy runs before `ck_done` is
allocated. That call then also solved the masked columns, which had already taken the
linearised step. A straight run cannot hit this, because its first call is always a full
call. A restarted run can. Fixed in `CkRstAfterAlloc`.

## 2. File format and backward compatibility

A marked block `CKRST001` (`kCkRstMagic`, `src/utils/two_stream_ck_rst.hpp`), in the same
form as PGENST01 / RTWARM01 / EINTRST1:
- The header (`CkRstHdr`: version, nslab, slab ids, host scalars) comes behind every
  other marked header.
- The nslab per-MeshBlock slabs come behind every other tail slab of each MeshBlock
  record.
- The reader stages the slabs (host `std::vector`s, by slab id). The first RT call
  consumes them: `CkRstAfterAlloc`, `CkRstBeginCall`, `CkRstCadRestore` in
  `src/utils/two_stream_ck_rst_state.hpp`.

Behaviour by case:
- **ck_implicit off**: no block is written, and the file is byte-identical to before
  (gate D, semi).
- **Old restart** (no block): it loads, prints one WARNING and takes the old path. That
  path is ck_thk all thin, no stored operator, and a cadence full call first. It is
  bitwise with the old binary (gate O).
- **Old binaries** cannot read a file that carries the block. The same holds for EINTRST1.

Extra size of each restart:
- T4: 1 slab.
- c2: 3 slabs, plus 7 (hydro) or 10 (MHD) when the next call re-applies.
- c2 + every = 4: 7 slabs. It never re-applies: full calls are 4 cycles apart, which is
  at least xstep = 2.

This is against about 13 slabs of hydro state.

Code:
- `src/utils/two_stream_ck_rst.hpp` (new, light): the format and the staging.
- `src/utils/two_stream_ck_rst_state.hpp` (new): collect, restore, snapshot, rebuild.
- `src/utils/two_stream_rt.hpp`: 4 hook lines plus the include.
- `src/outputs/restart.cpp`, `outputs.hpp`: write the block.
- `src/pgen/pgen.cpp`: read it.
- `two_stream_column_ck.hpp`: comment only.

No change to `mhd_tasks.cpp`, and none to the ck kernel storage.

## 3. CPU gates (`restart/gate.sh` -> /viper/ptmp2/jinma/ckrst_0924/gate/gate.out)

Setup:
- Binaries (git archive snapshots, gcc/14 openmpi/5.0, Release): `athena.cpu.base` =
  e89954e2, `athena.cpu.r2` = this branch.
- Input: `restart/wp_rst.athinput`. This is the well-posed full-hydro transient
  (cadence/wp_cad.athinput) with `data_format = %.17e` on the history.
- Each run writes hst and bin every cycle, and rst at exit.
- `restart/cmp.py` compares the rst bytes after `<par_end>`, the hst rows matched by
  time, and the bin dumps matched by content.

| gate | arms | result |
|---|---|---|
| R: 8 straight vs 4 + 4 and vs 3 + 5; 6 straight vs 3 + 3 | semi, t4 (per RK stage), c2, c2 + every 4, c2 + every 4 + thr 0.02 | **all BITWISE** (rst, hst, bin). The c2 3 + 3 and 3 + 5 restarts land on an xstep re-apply and exercise the rebuild ("ck restart: rebuilt the ck_impl_xstep operator stored on cycle 2"). The c2e4 3 + x restarts land on a linearised step. |
| D: base vs new, 6 cycles straight | same | **BITWISE** (hst, bin; semi also rst bytes) |
| O: a base-written (old) restart read by base and by new | same | **BITWISE** (hst, bin); the WARNING is printed |

Before the guard fix (binary r1), the c2e4g 3 + x restarts differed (eint 8.6e-3 in the
24 masked columns on the first step). Everything else was already bitwise on r1.

## 4. GPU gate (`restart/gpugate.sub`, apudev, 2 ranks, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1)

Setup:
- Binary: `athena.gpu.r2` (same source as athena.cpu.r2, ROCm 6.3, gfx942), md5
  417aae040aaaddb7dfbfa3df292c415a.
- Restart: bench/cs_hyd4_prod/rst/dhj.00160.rst, read in place: t = 2.43883e7 s,
  **rotation 80.0**, cycle 1256387.
- Input: `restart/prod_rst.athinput` (cadence/prod_cad.athinput with the history at
  %.17e).
- Runs: 1660 cycles straight (t = 2.44186e7 s, 0.099 rotation) against 830 + restart
  + 830 (about 0.05 rotation each).
- Compared: the final rst payload, the hst every 1500 s, and the bin dump.

| job | arm | result |
|---|---|---|
| 11956954 | c2 + every = 4 | rst payload, 10 hst rows, bin: **BITWISE** |
| 11956955 | c2 (every = 1) | rst payload, 10 hst rows, bin: **BITWISE** |

Logs: /viper/ptmp2/jinma/ckrst_0924/gpugate/log.out.1195695{4,5}, and the runs are in
`gpugate/<arm>/{s,a,r}`.

## 5. On the merge with rt-integration 9c1a12d7 (ck-scratch, mhd-split), and MHD

Merge commit 8dbab6bb. CPU binaries `athena.cpu.m1` (the merge) and `athena.cpu.m2`
(ccb9712d, the fix below).

**MHD gate** (`restart/mhdgate.sh`, 12 ranks):
- Setup: bench/cs_mhd_prod4/rst/dhj.00127.rst, read in place: t = 1.93151e7 s,
  **rotation 63.3**, cycle 1192186. Input `restart/prod4_mhd_rst.athinput` (the prod4
  input plus the ck keys). Arm: c2 + every = 4. Runs: 6 cycles straight vs 3 + restart
  + 3.
- m1 result: **DIFFER**. The three cycles after the restart are all linearised steps,
  so no implicit call allocated `ck_thk`. The restart written at their end then carried
  no ck block at all, and the staged state was lost.
- Fix (ccb9712d): `CkRstCollect` passes on whatever is still staged: THK, CAD0-3, CV and
  the XS snapshot.
- m2 result: rst payload **BITWISE** (mhdgate/gate.out).

**Hydro gates on m2** (gate_m2/gate.out):
- R: wpc 4 + 4 and 3 + 5, and wpf 3 + 3, for semi, c2, c2e4 and c2e4g: all BITWISE.
- D (against e89954e2) and O (old restart): BITWISE.

GPU: athena.gpu.r2 predates the ccb9712d fix. The fix only matters for a restart written
before the first implicit call after a restart. It does not affect the GPU gates, nor the
ab3 chain, whose segment-1 restart comes after 25000 cycles.
