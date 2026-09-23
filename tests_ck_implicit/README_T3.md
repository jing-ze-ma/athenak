# T2 + T3: frozen operator and linear re-apply kernel for implicit tm (2026-09-23)

Phases T2 and T3 of `DESIGN_tm.md`. Not committed. The patch is
`bench/impl_t3_0923/t3.patch` (T1 + T2 + T3 against HEAD `eb060b1b`, 1057 lines). It touches
`src/utils/two_stream_rt.hpp`, `src/utils/two_stream_column_ck.hpp` and
`src/pgen/deep_hot_jupiter_rt.cpp`. Everything was built in `bench/impl_t3_0923/src_new`, a
`git archive HEAD` snapshot with `impl_t1_0922/t1.patch` applied (it applied cleanly).

## Verdict

* **Off path: unchanged.** On CPU, production tm + sph + beam is bitwise against HEAD, and
  so is tm without the beam. On GPU, all 118 T1 chain kernels have identical instruction,
  scratch, VGPR, spill and private-segment counts.
* **T2 is exact.** It is bitwise against T1 in both Newton modes.
* **T3 agrees with T2 to round-off.** On the same `B`, Src agrees to 4e-15 and Fb to 1.3e-15,
  and Em is identical. Every call takes the same number of passes.
* **Convergence (restart `rst/dhj.00135.rst`, `arat = 1e30`, maxit 8, rj 1, seed 2).**
  After a 15-cycle restart transient: 0 of 260 calls fail to converge, the mean is 6.20
  passes, and the gap is max 7.85e-7 / median 4.96e-7. **The 1e-6 target is met.**
* **Cost, GPU, same binary:** T1 is **4.58x** the semi-implicit step, T2 **3.20x** and T3
  **2.73x**. The linear kernel costs **0.114** of a production sweep, and a whole linear
  pass costs 0.30. The design estimate was 1.75x. The gap to it comes from the two full
  passes per call (46.8 ms) and from `ck_impl_tri` (1.95 ms per pass), not from the new
  kernel (section 4).

## 1. The change

**T2.** `launch_ck_form` now instantiates `FOP = 1` for `FRM = 1`, with and without `JAC`
(4 more kernels per tier; the `JAC` ones only up to tier 264). The tm body already honoured
`ckfus`. The FATAL on `ck_impl_frozen_op` + tm is gone.

New rule: under `FOP` the cut `icut` is frozen for passes > 0, because `rt_pre_cut` is
skipped. The reason is that a cut which moved down by one cell would read a `kappa rho`
that this call never stored. Without `FOP` nothing changes.

**T3.** `problem/ck_impl_lin` (default false) needs `ck_implicit`, `ck_sweep_form = 1`,
`frozen_op`, `frozen_cof` and `!refresh_kappa`; anything else is a FATAL. It works in three
steps:

* On the storing pass, `ck_lin_build` fills two packed arrays from what that pass stored:
  * `lP`, 9 slots per (cell, chain): `e0`, `cin`, `cout`, `R`, `1/(1+R beta)`, `dslv/dB`,
    `dsuu/dB`, `dt_l/dtc` and `2 wfc/mu kappa rho`;
  * `lG`, 4 slots per cell: `beta`, `rat`, `fsc`, `1/dz`.
* Every pass that does not assemble the Jacobian runs `rt_chain_ck_lin1` instead of the
  chain kernel. It is one thread per chain: the Sc forward substitution plus the ray back
  substitution. It does no table look-up, no exponential, no divide and no beam.
* `rt_chain_ck_lin_sum` then sums the block's chains into `Src_g`/`Fb_g`/`Em_g`.

`ck_impl_lin_thr = 4` selects the block-threaded variant instead. `ck_impl_lin_check = 1`
runs the chain kernel on the same `B` and prints the difference. `ck_impl_debug = -1`
prints the report line from every rank (from `stall_debug.patch`).

**Memory:** +9 x 446 MB for `lP` and +2 x 446 MB for the partials per rank at production,
plus T2's 4 x 446 MB.

## 2. Correctness gates (CPU, `bench/impl_t3_0923/gate/`, 20 cycles, serial)

| gate | arms | result |
| --- | --- | --- |
| (1) off path | HEAD `athena.cpu.head` vs `athena.cpu.final`: production tm + sph + beam, and tm without beam | **BITWISE**, 25 files each (`off_final.log`) |
| impl unchanged | T1 vs new: four-pass implicit; tm implicit with rj 1 + seed 2 | BITWISE (`final.log`) |
| T2 | T1 tm (rj 1, seed 2) vs + `frozen_op`; T1 tm (rj 0, seed 0) vs + `frozen_op` | **BITWISE** both |
| (2) T3 | `lin_check`, 212 passes | max\|lin-chain\|/max\|chain\|: Src 3.5e-15, Fb 8.9e-16, Em 0 |
| (2) T3 | T2 vs T3, 40 calls | pass counts identical in every call; res agrees to 1.3e-7 relative |
| (2) T3 | T2 vs T3, restart payload after 20 cycles | density/energy words 2.6e-12; near-zero momentum words 9e-9 (`rstcmp2.py`) |
| (2) threading | thr 1 vs thr 4 | final states identical |

(1b) GPU code object: `counts2.v3.txt` against `impl_t1_0922/new.counts2.txt`.
* **All 118 T1 kernels are identical**, the production one (FOP=F, IMP=F) included; T1 had
  already shown that one identical to HEAD.
* 14 FOP=T kernels are new.
* `rt_chain_ck_lin1` has 89 VGPR, 0 VGPR spills, 4 SGPR spills and 1104 B private (Sc
  column). `lin_sum` has 62 VGPR.
* A first version with 5 separate Views ran at 128 VGPR with 81 VGPR spills (`meta.t3.txt`).
  Packing into `lP`/`lG` removed them.

(2) GPU: `lin_check`, both ranks, 793 passes: Src 4.0e-15, Fb 1.3e-15, Em 0
(`gpu/j3/k_chk1`, job 11944461).

## 3. Restart choice and convergence (gate 3)

**Scan of prod3.** 143 bins and 568 restarts were scanned. The statistic is osc =
\|T_i - (T_{i-1}+T_{i+1})/2\|/T above 1e-3 bar, with T taken from the general-EOS table
(`bench/impl_t3_0923/scan/oddeven_scan.py`, `scan_rst_*.txt`).
* The detector fires on the known column (rst 567, gid 13 corner: osc 1.82).
* Only rst 0 is fully clean. From rot 0.5 on, 180-340 columns sit at osc > 0.1 (fronts near
  1e-6 and 1e-8 bar), and the corner column spikes above 1 intermittently (rot 38, 88, ...,
  283).
* **Choice: `rst/dhj.00135.rst` (rot 67.5, ncycle 1153131).** It is the latest restart with
  no column above 0.3 (max osc 0.290, 220 columns above 0.1). Rst 567 has 9 columns above
  0.3 and max 1.82.

**Convergence** (GPU, job 11944460, `athena.gpu.v3`, 150 cycles = 300 calls, both ranks;
`ana.py` logic):

| window | calls | non-conv | passes (mean, max over ranks) | res max | gap max / median |
| --- | --- | --- | --- | --- | --- |
| all | 300 | 20 | 6.39 | 2.7e-7 | 9.81e-7 / 4.96e-7 |
| after cycle +20 | 260 | **0** | 6.20 | 1.0e-8 | **7.85e-7 / 4.96e-7** |

* All 20 non-converged calls are in the first 15 cycles: the switch from prod3's RT to the
  prod4 form.
* T1, T2 and T3 give the same table.
* Semi-implicit gap (`rt_desum`, 32 lines): mean -10.2 %, range -6.4 to -14.4 %.
* `ck_impl_seed = 0` (store and Jacobian in one pass) does **not** converge: 153 of 260
  calls fail, gap max 5.0e-6. Seed 2 stays.

## 4. Cost (GPU, apudev, 2 ranks / 2 GPUs, `cpu time used`, 150 cycles from rst 00135)

Job 11944460, one binary (`athena.gpu.v3`), interleaved s, t1, t2, t3, t3s, t3s, t3, t2,
t1, s:

| arm | r1 / r2 [s] | cycles/s | x semi |
| --- | --- | --- | --- |
| semi-implicit | 11.29 / 11.54 | 13.14 | 1.00 |
| T1 (implicit tm) | 52.68 / 52.02 | 2.87 | 4.58 |
| T2 (+ frozen_op) | 36.53 / 36.53 | 4.11 | 3.20 |
| **T3 (+ lin, thr 1)** | **30.83 / 31.51** | **4.81** | **2.73** |

Earlier binaries, same arms:
* Job 11944353: s 12.18 / 12.28, T1 56.80 / 53.95, T2 39.14 / 38.67.
  * Unpacked thr 1: 36.05 / 36.02.
  * thr 4: 37.54 / 37.45.
* Job 11944227 (thr 4, unpacked): s 11.32 / 13.26, T1 53.19 / 51.46, T2 36.21 / 38.21,
  T3 35.53 / 35.85.

**Per pass** (rocprofv3, rank 0, 40 cycles = 80 calls; `gpu/prof.py`; T1 from job
11944260, production and T2 from 11944354, T3 from 11944461):

| kernel | ms per launch |
| --- | --- |
| production chain kernel (semi-implicit) | 17.45 |
| T1 pass / T1 Jacobian pass | 15.26 / 31.62 |
| T2 re-apply pass (FOP), derived | ~6.7 |
| T2/T3 storing pass (with beam) / Jacobian pass (FOP) | 22.6 / 24.2 |
| **T3 `rt_chain_ck_lin1` + `lin_sum`** | **1.46 + 0.53 = 1.99 (0.114 of a sweep)** |
| T3 thr 4, packed | 2.88 |
| `ck_lin_build` (once per call) | 2.93 |
| per pass, any arm: `ck_impl_tri` / apply / `ck_impl_res` / `rt_pre_opac` | 1.95 / 0.57 / 0.50 / 0.18 |

**Where a T3 call goes** (rank 0, all kernels 89.8 ms, against 26.6 ms semi-implicit):

| item | ms per call |
| --- | --- |
| storing pass + Jacobian pass | 46.8 |
| `ck_impl_tri` (6.3 per call) | 12.2 |
| linear passes (5.05 per call) | 10.0 |
| apply + res | 7.5 |
| build | 2.9 |

The next levers are two, both outside T3:
1. Assemble the Jacobian from the stored factorisation, instead of in a full sweep.
2. Coalesce `ck_impl_tri`, whose `(m,k,j,i)` arrays are read with lanes on `j`.

## 5. Files (`bench/impl_t3_0923/`)

* `t3.patch`, `build.sh`, `src_new/` (snapshot and builds), `src_head/`.
* Binaries: `athena.cpu.{head,v3,final}`, `athena.gpu.v3` (timing binary; the final source
  differs from it only in comments and one split string literal).
* Code-object tables: `counts2.{t1,t3,t3b,v3}.txt`, `meta.*.txt`, `dis2.sh`.
* `gate/`: `gate.sh`, `cmpall.py`, `cmprel.py`, `rstcmp*.py` and the logs.
* `gpu/`: `submit{,2,3}.sh`, `submit_prof{,2,3}.sh`, `smoke.sh`, `ana.py`, `prof.py`,
  `log.out.*`, and the arm directories `j2/`, `j3/`.
* `scan/`: the odd-even scan.
