# T4 vs T3 vs semi-implicit from prod4 and hyd4 restarts (2026-09-23)

This redoes the GPU comparison of `README_T4.md`, but starts from the newest production
restarts instead of prod3's `rst/dhj.00135.rst`. The user asked for this because the prod3
restarts carry odd-even columns. No code was changed and nothing was committed.

There are two states:

* **Primary: hyd4.** `bench/cs_hyd4_prod/rst/dhj.00023.rst` is at rot 11.50 (t = 3.5075e6 s,
  ncycle 158989, dt 19.22 s). It runs with the hyd4 input. This is the more developed state,
  and its dt is near the production dt.
* **MHD check: prod4.** `bench/cs_mhd_prod4/rst/dhj.00003.rst` is at rot 1.50
  (t = 4.5750e5 s, ncycle 118114, dt 2.93 s). It runs with the prod4 input. **This state is
  early:** prod4 is still spinning up from its cold start, and dt is about 3 s against
  13-20 s at steady state.

Both restart files were copied into `bench/impl_t4p4_0923/`, and the md5 of each copy matches
its production file. Nothing was written into either production directory.

## Verdict

* **T4 stays the cheaper and faster-converging implicit arm on both production states.**
  * hyd4: T4 costs **1.77x** the semi-implicit step and T3 **2.48x**, so T4 is 0.71 of T3.
  * prod4: T4 costs **1.56x** and T3 **2.05x**, so T4 is 0.76 of T3.
  * For reference, from prod3 rst 135 (README_T4) the numbers were 1.70x and 2.46x.
* **Convergence after the first 20 cycles (260 calls): 0 non-converged calls in either arm
  or state.**
  * Passes per call: T4 4.50 vs T3 6.16 (hyd4), and T4 4.75 vs T3 6.00 (prod4).
  * The gap is smaller for T4 than for T3 in both states.
* **The restarts are not clean either.** Section 2 has the scan.
  * Both restarts flag more kinked columns than prod3 rst 135 (osc > 0.1: 1337 hyd4 and
    1359 prod4, against 220).
  * No column carries a long odd-even train.
  * prod4 has one ragged patch at the top of gid 2 (k 1-5, j 9-15, 3e-9 .. 2e-8 bar),
    where T jumps between the 200 K floor and about 1000 K.
* **T4 vs T3 after 150 cycles: large local differences remain. They sit in the kinked or
  ragged cells, and the tight references show that tolerance sets them, not the scheme.**
  * **hyd4:**
    * T4 vs T3: max 6.9e-2, rms 3.1e-4.
    * The worst cell is gid 1 k14 j0 i62 at 1.5e-6 bar, next to hyd4's worst kinked
      column. There T4 equals both tight references to 1e-5, and T3 is 7.4 % off.
    * Most other cells above 1e-2 lie above 1e-8 bar. There T3 differs even from its own
      repeat by up to 4.0e-2 (GPU atomics; T4 is bitwise reproducible). The two tight
      runs differ there by 2.9e-2.
  * **prod4:**
    * T4 vs T3: max 2.0e-1, rms 7.3e-4. 75 of the 81 cells above 1e-2 are in the gid 2
      patch.
    * The two tight references agree to 1.1e-3 (0 cells above 1e-2).
    * T3 and T4 are **equally** far from them: max 1.18 and 1.26, 173 and 168 patch cells
      above 1e-2, both at gid 2 k3 j13 i82.
    * In that patch the 1e-8 tolerance answer differs from the 1e-10 answer by O(1),
      whichever Jacobian is used.

## 1. Setup

* **Binary: `athena.gpu.v4`** (md5 `85aab163...`), copied from `bench/impl_t4_0923`. It was
  built from `d14f96de` plus `t4.patch`. The md5s of its three dhj sources
  (`two_stream_column_ck.hpp`, `two_stream_rt.hpp`, `deep_hot_jupiter_rt.cpp`, in
  `src_md5.v4.txt`) equal `git show HEAD:<file>` at HEAD `b35f6c66` / `95595afc`. Between
  `d14f96de` and HEAD, the only other source changes are in `src/rad_m1/`, which dhj does
  not construct. No rebuild was needed.
* **Job layout:** GPU apudev, 2 ranks / 2 GPUs on 1 node, with `HSA_XNACK=1` and
  `HSA_NO_SCRATCH_RECLAIM=1`. Each run is 150 cycles, and the cost is `cpu time used`.
* **Arm order:** `s t3 t4 | t4 t3 s | s t3 t4`, then the tight references `t3x t4x` once
  each.
* **Jobs:** 11945062 (hyd) and 11945063 (mhd). Jobs 11945053/54 are void: the `STATE`
  environment variable did not reach the job, so it is now passed as an argument.
* **Inputs:** `gpu/hyd.athinput` and `gpu/mhd.athinput` are the production inputs,
  byte-for-byte, with T4's appended `<problem>` block (`gpu/impl_block.txt`). The block
  declares the `ck_impl_*` keys: `maxit 8`, `reuse_jac 1`, `seed 2`, `tol = dtol = 1e-8`,
  `rt_cell_report` every 20 cycles.

| arm | switches (on top of production tm + sph + beam) |
| --- | --- |
| `s` (semi-implicit, production) | none (`ck_impl_debug=-2` only) |
| `t3` | `ck_implicit`, `ck_impl_arat=1e30`, `ck_impl_frozen_op`, `ck_impl_lin`, `lin_thr=1` |
| `t4` | `t3` + `ck_impl_fuse` + `ck_impl_jac_lin` + `ck_impl_cvsec` (no `lw`, as briefed) |
| `t3x`, `t4x` | as `t3` / `t4` with `ck_impl_tol = dtol = 1e-10`, `maxit = 20` (as README_T4) |

## 2. Odd-even scan of the two restarts

The statistic is T3's: osc = |T_i - (T_{i-1}+T_{i+1})/2| / T above 1e-3 bar, per column.
The T3 tool cannot read these restarts directly, because prod4 and hyd4 restarts carry 3
extra cache arrays per MeshBlock (112.6 MB vs 82.2 MB). Two scripts handle them:

* `scan_inv.py` reuses the T3 tool's potential, geometry, EOS inversion and metric on the
  new record layout. Check: on prod3 rst 135 it reproduces T3's line exactly
  (`maxosc 0.290`, `n>0.1 220`).
* `scan4.py` uses the cached T instead. It gives the same result: 1.131 vs 1.134.

| restart | max osc | cols > 0.1 | cols > 0.3 | s3 / s4 (3- / 4-cell alternating runs) | worst column |
| --- | --- | --- | --- | --- | --- |
| prod3 rst 135 (T3/T4's state) | 0.290 | 220 | 0 | 276 / 221 | gid 4 k0 j0 i69, 5.3e-6 bar |
| **hyd4 rst 23 (rot 11.5)** | 0.681 | 1337 | 20 | 320 / 86 | gid 1 k14 j1 i60, 1.1e-6 bar, T 661 K |
| **prod4 rst 3 (rot 1.5)** | 1.134 | 1359 | 113 | 710 / 314 | gid 2 k2 j11 i86, 2.5e-9 bar, T 200 K (floor) |

By pressure band (`scan_band.txt`), columns > 0.1 / > 0.3:

| band [bar] | hyd4 | prod4 |
| --- | --- | --- |
| 1e-5 .. 1e-3 | 312 / 3 | 716 / 70 |
| 1e-7 .. 1e-5 | 302 / 20 | 118 / 30 |
| 1e-9 .. 1e-7 | 983 / 0 | 725 / 13 |
| < 1e-9 | 0 / 0 | 23 / 2 |

**Neither restart is free of radial T kinks.**

* Both flag more columns than prod3 rst 135, which T3 picked as its cleanest usable
  restart.
* No column carries a long odd-even train. The longest alternating run is 5 pairs
  (6 cells) in prod4 and 4 pairs in hyd4 (`scan_prof.txt`). These are ragged fronts, not
  prod3's persistent corner column.
* **prod4 has one concentrated patch:** gid 2, k 1-4, j 10-14. It holds the 8 worst
  columns of the whole mesh, at about 3e-9 .. 2e-8 bar. There T jumps between 200 K (the
  floor) and 1000 K from cell to cell.
* In hyd4 the kinks are spread thinly; the worst sit in gid 1 near 1e-6 bar.

## 3. Cost (GPU, `cpu time used`, 150 cycles, 3 interleaved repeats)

| arm | hyd4 r1 / r2 / r3 [s] | x semi | prod4 r1 / r2 / r3 [s] | x semi |
| --- | --- | --- | --- | --- |
| semi-implicit | 10.00 / 9.77 / 9.73 | 1.00 | 10.60 / 10.86 / 10.89 | 1.00 |
| T3 | 24.17 / 24.52 / 24.56 | 2.48 | 21.97 / 22.13 / 22.19 | 2.05 |
| **T4** | **17.44 / 17.40 / 17.39** | **1.77** | **16.67 / 16.84 / 17.08** | **1.56** |
| T3x (tol 1e-10, 1 run) | 35.84 | 3.64 | 24.92 | 2.31 |
| T4x (tol 1e-10, 1 run) | 18.75 | 1.91 | 17.57 | 1.63 |

The repeats spread by 1-3 %. Sources: `gpu/log.out.11945062` (hyd4), `gpu/log.out.11945063`
(prod4), `gpu/ana5_{hyd,mhd}.txt`.

## 4. Convergence and gap (after the first 20 cycles: 260 calls, both ranks)

* Script: `gpu/ana5.py`, which applies `ana4.py`'s logic.
* Gap:
  * implicit arms: |ckdesum|, max and median over the calls;
  * semi-implicit: `rt_desum rel`.

| state / arm | passes mean / max | non-conv | res max | gap max / median | median res by pass |
| --- | --- | --- | --- | --- | --- |
| hyd4 T3 | 6.16 / 8 | 0 | 1.0e-8 | 3.81e-6 / 1.67e-6 | 1.1e-2 1.9e-3 5.4e-5 1.8e-6 7.2e-8 1.4e-8 1.0e-8 1.0e-8 |
| **hyd4 T4** | **4.50 / 5** | **0** | 1.0e-8 | **1.14e-6 / 6.03e-7** | 1.1e-2 1.9e-3 5.6e-6 1.4e-8 1.0e-8 |
| hyd4 T3x | 13.50 / 20 | 130 (maxit 20) | 2.8e-10 | 3.3e-9 / 1.7e-9 | crawls, 0.7 per pass below 1e-8 |
| hyd4 T4x | 5.70 / 7 | 0 | 1.0e-10 | 1.25e-8 / 6.3e-9 | |
| prod4 T3 | 6.00 / 7 | 0 | 1.0e-8 | 1.34e-5 / 9.85e-6 | 2.6e-3 9.3e-5 4.8e-7 2.7e-8 8.4e-8 3.7e-8 9.9e-9 |
| **prod4 T4** | **4.75 / 6** | **0** | 1.0e-8 | **7.22e-6 / 4.73e-6** | 2.6e-3 9.3e-5 4.8e-7 1.0e-8 1.0e-8 9.9e-9 |
| prod4 T3x | 9.02 / 15 | 0 | 1.0e-10 | 1.62e-8 / 1.58e-8 | |
| prod4 T4x | 6.55 / 15 | 0 | 1.0e-10 | 1.29e-7 / 6.99e-8 | |

* Repeats r2 and r3 give the same passes, non-converged counts and gap max for both T3 and
  T4.
* The whole window, including the first 20 cycles, is also clean: T3 on hyd4 has 0
  non-converged calls in 300. Starting from a prod4-form restart there is no prod3 -> prod4
  transient, unlike README_T3/T4.
* **Semi-implicit gap** (`rt_desum rel`, after +20 cycles):
  * hyd4: 24 lines, mean -10.4 %, range -6.5 to -15.0 %.
  * prod4: 28 lines, mean -1.4 %, range -0.9 to -2.0 %.
* The prod4 implicit gaps (median 5e-6 for T4, 1e-5 for T3) are above the 1e-6 target of
  README_T3. On hyd4, T4's median meets it (6.0e-7) and T3's does not (1.7e-6).

## 5. T4 vs T3 after 150 cycles, and the tight references

* Script: `gpu/tcmp.py` (output in `gpu/tcmp.txt`).
* Quantity: relative difference of the cached T on all active cells, final restart of each
  arm.

| pair | hyd4 max / rms / cells > 1e-2 | prod4 max / rms / cells > 1e-2 |
| --- | --- | --- |
| **T3 vs T4** | 6.88e-2 / 3.12e-4 / 164 | 2.02e-1 / 7.30e-4 / 81 |
| T3 r1 vs r2, r1 vs r3 | 4.04e-2, 3.74e-2 / ~2e-4 / 64, 81 | 4.6e-11 / 2e-13 / 0 |
| T4 r1 vs r2 | 0 (bitwise restart) | 0 (bitwise restart) |
| T3x vs T4x (tight vs tight) | 2.86e-2 / 1.97e-4 / 81 | 1.09e-3 / 1.96e-6 / 0 |
| T3 vs T3x | 6.88e-2 / 2.67e-4 / 115 | 1.18 / 3.00e-3 / 178 |
| T4 vs T4x | 4.02e-2 / 2.95e-4 / 158 | 1.26 / 3.06e-3 / 172 |
| T3 vs T4x | 6.88e-2 / 2.77e-4 / 119 | 1.18 / 3.00e-3 / 178 |
| T4 vs T3x | 4.02e-2 / 2.77e-4 / 132 | 1.26 / 3.06e-3 / 172 |

**Where the differences are.**

**hyd4:**
* The T3 vs T4 maximum is at gid 1 k14 j0 i62, 1.5e-6 bar. That is next to the worst
  kinked column of the hyd4 restart (gid 1 k14 j1 i60, 1.1e-6 bar).
* T there, in code units:

  | T3 | T4 | T3x | T4x |
  | --- | --- | --- | --- |
  | 7.236e10 | 6.738e10 | 6.738e10 | 6.738e10 |

  **At the one large mid-atmosphere difference, T4 is right and T3 is the outlier.**
* 133 of T3 vs T4's 164 cells above 1e-2 sit above 1e-8 bar, 117 of them in gid 0-3.
  There the T3 repeats differ by up to 4e-2 and the two tight runs by 2.9e-2, so these
  cells amplify round-off.
* Away from the top and gid 0-3, the rms is:
  * T3 vs T4: 1.1e-4
  * T3 vs T3x: 9.4e-5
  * T4 vs T4x: 1.4e-4
  * T3x vs T4x: 1.5e-5

**prod4:**
* 75 of T3 vs T4's 81 cells above 1e-2 lie in gid 2, k 1-4, j 11-14, i 69-91,
  1.7e-9 .. 2e-7 bar. This is the ragged patch the scan found at the restart.
* The two tight runs agree there. T3 and T4 both miss them by O(1):
  * worst cell gid 2 k3 j13 i82: T3 1.72e10, T4 1.66e10, T3x = T4x = 3.76e10;
  * 173 and 168 cells above 1e-2 in the patch.
* **Neither is closer.** The difference is set by `tol = 1e-8` in a front that the 1e-8
  and 1e-10 answers resolve differently. It is not set by the T3 vs T4 Jacobian.
* Outside the patch and the top, the rms is:
  * T3 vs T4: 8.2e-5
  * T3 vs T3x: 4.4e-5
  * T4 vs T4x: 9.4e-5

**Answer to the open item of README_T4.**
* A large local T4 vs T3 difference still appears without the prod3 odd-even column. It
  again appears where the restart carries a kinked or ragged front: hyd4 near 1e-6 bar, and
  prod4 in the gid 2 top patch.
* The tight references show that it is the tolerance-level answer in those cells that is
  ill-determined. Neither scheme has a bias:
  * In hyd4, T4 matches the references at the worst cell.
  * In prod4, both schemes miss them equally.
* Elsewhere, T4 and T3 agree to about 1e-4 rms, at the level of each one's own tolerance
  error.

## 6. Files (`bench/impl_t4p4_0923/`)

* `dhj.00003.rst` and `hyd.00023.rst`: copies of the production restarts (md5 checked).
* `athena.gpu.v4`: the binary.
* Scan:
  * scripts `scan4.py` (decoder and cached-T scan), `scan_inv.py` (the T3 tool with EOS
    inversion), `scan_band.py`;
  * outputs `scan4.txt`, `scan_inv.txt`, `scan_band.txt`, `scan_prof.txt`,
    `*.Tinv.npy`.
* `gpu/`:
  * scripts `submit.sh`, `hyd.athinput`, `mhd.athinput`, `impl_block.txt`, `ana5.py`
    (uses `ana4.py`), `tcmp.py`;
  * outputs `ana5_{hyd,mhd}.txt`, `tcmp.txt`, `log.out.11945062/3`;
  * run directories `hyd/`, `mhd/`: `{s,t3,t4}_r{1,2,3}`, `t3x_r1`, `t4x_r1` (2.5 GB,
    including the final restarts).
