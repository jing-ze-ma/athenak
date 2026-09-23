# T4: making the implicit correlated-k call cheaper on the GPU (2026-09-23)

Phase T4, following `README_T3.md`. Nothing is committed. The patch is
`bench/impl_t4_0923/t4.patch` (792 lines). It applies to `d14f96de` and to the current
HEAD `a6202226` (`git apply --check`). It touches `src/utils/two_stream_column_ck.hpp`,
`src/utils/two_stream_rt.hpp` and `src/pgen/deep_hot_jupiter_rt.cpp`.

It was built in `bench/impl_t4_0923/src_new`, a `git archive` snapshot of `d14f96de` with
the patch applied (`src_head` is the same snapshot without it). The timing binary is
`athena.gpu.v4` (md5 `85aab163...`); `src_md5.v4.txt` holds the md5s of the three source
files it was built from. The final CPU binary is `athena.cpu.v5`, built from the same
source.

## Verdict

* **T4 costs 1.70x the semi-implicit step. T3 costs 2.46x.** Same binary, 3 interleaved
  repeats (job 11945029).
  * semi-implicit: 10.83 / 10.72 / 10.72 s
  * T3: 26.38 / 26.47 / 26.40 s
  * T4: 18.20 / 18.49 / 18.31 s
  * T4 is 0.69 of T3. The design goal was about 1.75x; this meets it.
* **T4 is three switches on top of T3.** Each is default off, and each was measured on
  its own (section 3):
  * `problem/ck_impl_fuse = true`: the residual test and the tridiagonal step run as one
    team-per-column kernel.
  * `problem/ck_impl_jac_lin = true`: the Jacobian is built from the stored
    factorisation. There is no JAC chain pass any more.
  * `problem/ck_impl_cvsec = true`: the heat capacity in the Jacobian rows is a secant.
  * `ck_impl_lw = true` was also on in the T4 arm. It changes nothing (lever d).
* **Gate 1, implicit off: the production path is unchanged.**
  * CPU: bitwise against HEAD in the production configuration (tm + sph + beam) and in tm
    without the beam.
  * CPU: T3's own implicit path is bitwise against HEAD.
  * GPU: all 133 chain-kernel entries have the same instruction, scratch, VGPR, spill and
    private-segment counts as T3.
* **Gate 2, the converged answer: met.**
  * Non-converged calls after cycle +20: 0 of 260.
  * Gap: max 6.02e-7 and median 4.84e-7. T3 has max 7.85e-7 and median 4.96e-7.
  * T4 and T3 converge to the same fixed point. After 1 or 10 post-transient cycles, T4
    is closer to a tightly converged run (tol 1e-10) than T3 is (section 2).
* **The pass count fell from 6.20 to 4.50 per call.** The lever is the secant heat
  capacity:
  * T3's Newton converges linearly, about 0.03 per pass, and then crawls near 1e-8.
  * With the secant cv the median residual goes 8.7e-4 -> 9.6e-7 -> 8.9e-9.

## 1. The change (all default off; T3 is recovered to the bit with every switch off)

| knob | what it does |
| --- | --- |
| `ck_impl_fuse` | `CkImplStep` becomes one `TeamPolicy` kernel with one 64-lane team per column. The residual, the rows, `CkThinSolve`, the caps and the fallback run in parallel over the cells, reading coalesced along i. Only the two Thomas recurrences run on one lane, out of LDS. The arithmetic and its order are T3's, so the state is bitwise T3's; only the two diagnostic sums (slots 4 and 5) are summed in a different order. There is one host read-back per pass instead of two. Refused with `ck_impl_debug > 0`. |
| `ck_impl_jac_lin` | `rt_chain_ck_jlin` runs the JAC recurrences of the tm body term for term (dSc/dB carried up, dd^+/dB and the reflection scalar W carried down) on what `ck_lin_build` stored in `lP`/`lG`. It uses one thread per chain and writes per-chain partials, which `ck_jlin_sum` adds in chain order with no atomics. The window parked at each face lives in the same partial arrays, so the kernel has no private column and no scratch. The residual of the Jacobian pass comes from the linear kernel. The result is deterministic: T4 repeats are bitwise. Needs `ck_impl_lin` with `lin_thr = 1`. |
| `ck_impl_cvsec` | The Jacobian rows use cv = (e_p - e_{p-1})/(T_p - T_{p-1}), cell by cell, from pass 1 on. The pass 0 -> 1 seed step gives a clean secant. The previous cv is kept wherever \|dT\| < 1e-7 T or the slope falls outside 1/50 .. 50 of e/T. T3 uses cv = e/T, which `README_T1_stall.md` measured at 5x off in the dissociating top. This is the grey mode-3 lesson: a linear term missing from the Jacobian. Needs `ck_impl_fuse`. |
| `ck_impl_jac0` | Builds the Jacobian on pass 0 (at e^n) instead of pass 1. This is lever (a). **Rejected.** |
| `ck_impl_jneg` | Keeps the negative per-chain off-diagonal parts and drops only a negative net entry. **Rejected.** |
| `ck_impl_lw` | Launches the implicit-only kernels light-weight (the functor goes as a kernel argument). **No effect.** |
| `ck_impl_debug = -2` | Adds `hist=res/dstep/active,...` for every pass to the per-rank report line. |

Memory:
* `jac_lin` adds 1 x 446 MB per rank (`ck_lpj`).
* `cvsec` adds 3 cell arrays.
* The fused kernel uses 6 x 8 x n1 B of LDS per team.

## 2. Gates

**(1) Implicit off: CPU bitwise.**
* Source: `bench/impl_t4_0923/gate/gate.sh`; logs `off_v5.log` and `gates_v5.log`.
* Setup: 20 cycles, serial; `.hst`, `.bin` payloads and `.rst` compared.
* Reference: `athena.cpu.head`, built from `d14f96de`. New: `athena.cpu.v5`.

| gate | result |
| --- | --- |
| production tm + sph + beam | BITWISE (25 files) |
| tm, no beam | BITWISE (25 files) |
| T3 (implicit tm, rj 1, seed 2, fop, lin): HEAD vs new | BITWISE |
| T3 vs T3 + `fuse` (with and without `lw`) | BITWISE; pass counts and residuals identical |

**(1b) GPU code object.**
* `counts2.v4.txt` (`dis2.sh dhj.v4.o`) is identical to T3's `counts2.v3.txt` on all
  133 chain-kernel entries, production included.
* `rt_chain_ck_jlin` uses 113 VGPR, 0 spills and 0 private bytes.
* The fused kernel uses 92 VGPR, 0 VGPR spills and 0 private bytes.

**CPU cold start** (`dhj_ck_implicit`, 40 calls; `cv_v5.log`, `gates_v3.log`):

| arm | passes / call | non-converged | notes |
| --- | --- | --- | --- |
| T3 | 7.30 | 22 of 40 | |
| `jac_lin` | 7.30 | 22 of 40 | Same pass count in every call as T3; res agrees to 1.3e-7. It matches T3's JAC to round-off. |
| `jac_lin + cvsec` | 5.50 | 4 of 40 | |
| `jneg` | 7.28 | 22 of 40 | |
| `jac0` | 7.90 | 25 of 40 | |

The `jac_lin` jac0 arm is bitwise equal to the chain-JAC jac0 arm, which confirms the new
assembly a second way.

**(2) Convergence on the production state.**
* GPU, `rst/dhj.00135.rst` with the prod4 `<problem>` block, as in T3.
* Window: cycles +20 to +150, 260 calls, both ranks.
* Script: `gpu/ana4.py`; runs `gpu/j4` (job 11944960) and `gpu/j6` (job 11945029).

| arm | passes mean / max | non-conv | res max | gap max / median | median res by pass |
| --- | --- | --- | --- | --- | --- |
| T3 | 6.20 / 8 | 0 | 1.0e-8 | 7.85e-7 / 4.96e-7 | 6.3e-3 8.7e-4 2.3e-5 6.9e-7 2.6e-8 1.5e-8 1.0e-8 1.0e-8 |
| + fuse | 6.20 / 8 | 0 | 1.0e-8 | 7.85e-7 / 4.96e-7 | as T3 |
| + fuse + jac_lin | 6.20 / 8 | 0 | 1.0e-8 | 7.85e-7 / 4.96e-7 | as T3 |
| **T4** (+ cvsec) | **4.50 / 5** | **0** | 1.0e-8 | **6.02e-7 / 4.84e-7** | 6.3e-3 8.7e-4 9.6e-7 8.9e-9 1.0e-8 |
| T4 + reuse_jac 0 | 4.50 / 5 | 0 | 1.0e-8 | 6.03e-7 / 4.84e-7 | as T4 |
| T4 + jneg | 4.50 / 5 | 0 | 1.0e-8 | **1.16e-6** / 9.19e-7 | 6.3e-3 8.7e-4 1.0e-6 1.4e-8 1.0e-8 |
| T3 + jac0 (chain) | 7.15 / 8 | **168** | 6.2e-5 | 4.75e-6 / 7.1e-7 | stalls at ~1e-6 |

The per-pass contraction of T3 is 0.02-0.05, and 0.23-0.93 in the tail. That is a
linearly converging iteration. Stage 1 of every step needs 7-8 passes and stage 2 needs
5.

**(2) The answer.**
* Source: `gpu/rstdec.py`, which decodes the restart records in double precision: e, rho,
  and the EOS caches T and p on the active cells.
* Start: T3's own state after 150 cycles (`gpu/rst_t3_150`), i.e. after the transient.
* Jobs: 11945010 (`h1`, 1 cycle) and 11945011 (`h10`, 10 cycles).
* The x arms use `tol = dtol = 1e-10` and `maxit = 20`.

| pair | 1 cycle: max rel T | 10 cycles: max rel T | 10 cycles: cells with dT/T > 1e-4 |
| --- | --- | --- | --- |
| T3 vs T3 repeat | 0 | - | - |
| T3 vs T3x (T3's own tolerance error) | 1.77e-4 | 4.22e-4 | 54 |
| T4 vs T4x (T4's own tolerance error) | 6.95e-5 | 1.57e-4 | 54 |
| T3 vs T4 | 1.48e-4 | 3.50e-4 | 115 |
| T3x vs T4x (same fixed point) | 1.06e-5 | 2.71e-5 | 0 |

T4 differs from T3 by less than T3 differs from its own converged answer. The two
tight runs agree to 3e-5.

Tolerance-level differences grow with time. After 150 cycles from rst 00135:
* T4 vs T3 is 5.6e-1 in T in 5 cells, with 1054 cells above 1e-4.
* All the large differences sit in one patch: MeshBlock 13, k ~ 7, j ~ 7, i ~ 103-118.
* T3 against its own repeat stays at 1e-8, because the GPU atomics perturb it only at
  round-off.
* T4 runs are bitwise reproducible.

The first cycle after the restart is the prod3 -> prod4 transient, and there both hit the
8-pass cap:
* T3's first call ends at res 7.4e-8.
* T4's first call ends at res 4.1e-7.
* Job 11945007 (`g1`).
* This window is excluded by the brief's "after the first 20 cycles".

## 3. Cost, lever by lever

* GPU apudev, 2 ranks / 2 GPUs (1 node), `cpu time used`, 150 cycles from rst 00135.
* `HSA_NO_SCRATCH_RECLAIM=1` is set in every job (coordinator, 09-23).
* One binary per job, with arms interleaved forward and then reversed.

| arm | job 11944960 (v3), r1 / r2 [s] | x semi | job 11944895 (v1*), r1 / r2 [s] |
| --- | --- | --- | --- |
| semi-implicit | 12.08 / 10.48 | 1.00 | 11.38 / 10.70 |
| T3 | 25.67 / 25.39 | 2.26 | 25.55 / 25.86 |
| + fuse | 22.19 / 21.78 | 1.95 | 22.26 / 21.87 |
| + fuse + lw | - | - | 22.05 / 23.01 |
| + fuse + jac_lin | 20.09 / 20.44 | 1.80 | - |
| **+ fuse + jac_lin + cvsec (T4)** | **18.26 / 18.24** | **1.62** | - |
| T4 without lw | 18.20 / 18.47 | 1.63 | - |
| T4 + reuse_jac 0 | 20.87 / 20.93 | 1.85 | - |
| T4 + jneg | 18.11 / 18.75 | 1.63 | - |
| T3 + jac0 (lever a) | - | - | 25.23 / 24.88 |
| coloured probes, 3 / 5 colours | - | - | 21.34 / 21.31, 22.23 / 22.41 |

\* v1 is `c02ddd34` plus the early patch, without the `d14f96de` bvals change. Its
probe-Jacobian arms (M applied to dB/dT on every 3rd or 5th cell) were dropped from the
patch: folding the far entries onto the three kept ones leaves 130 and 50 of 260 calls
non-converged.

Repeats of the final pair:
* Job 11944978 (v3): s 10.33 / 10.83 / 10.44, T3 25.67 / 25.61 / 25.63, T4 18.43 /
  18.29 / 18.22. That is 1.74x and 2.43x.
* Job 11945029 (v4, the final source): 1.70x and 2.46x.
* The spread of the semi-implicit repeats (10.3-12.1 s) is the largest error bar.

**Per call, rank 0** (rocprofv3, 40 cycles, `gpu/pf4` job 11945028; semi-implicit from
`gpu/pf`, job 11944979). All kernels:

| | semi | T3 | T4 |
| --- | --- | --- | --- |
| all kernels | 26.8 | 76.6 | 51.1 |
| chain store pass | 17.8 (production) | 17.1 | **22.0** |
| JAC chain pass | - | 18.2 | - |
| `ck_impl_tri` + `ck_impl_res` | - | 11.4 + 3.3 | fused 1.4 |
| `lin1` + `lin_sum` | - | 7.3 + 2.6 | 6.7 + 2.5 |
| `jlin` | - | - | 2.7 |
| apply | - | 4.0 | 3.1 |
| build | - | 2.8 | 2.8 |
| sweeps per call (incl. transient) | 1 | 7.05 | 5.86 |

### Open: the store pass is 5 ms slower in every arm without the JAC chain pass

* The same kernel with the same dispatch runs at 17.1 ms (T3, fuse) and 22.0 ms (every
  `jac_lin` arm) on rank 0, and at 18.8 vs 24.6 ms on rank 1 (`gpu/pf5`, job 11945037).
* T3's first call is also 22 ms; from call 2 on it is 17 ms.
* It is not the jlin scratch: v3's jlin had 3.3 kB of private memory per thread and v4's
  has none, and both are slow.
* It is not the ROCr scratch limit: `HSA_SCRATCH_SINGLE_LIMIT=16 GB` changes nothing
  (`gpu/pf6`, job 11945044).
* No other kernel overlaps it, and its predecessors are identical.
* It is worth about 1.5 s of the 18.3 s. The cause is not found.

## 4. The levers the brief listed

* **(a) Merging the store pass and the Jacobian pass: no.** They are not at the same
  state: the store is at e^n and the Jacobian at the seeded iterate. A Jacobian built at
  e^n (`jac0`) stalls at ~1e-6: 168 of 260 calls do not converge, and it saves only 0.4
  s. What the brief wanted from it is delivered by `jac_lin`, which drops the second
  full sweep: 18.2 ms becomes 1.4 (lin1) + 0.5 (sum) + 2.7 (jlin) ms.
* **(b) Pass count: linear, not quadratic.**
  * The missing term is cv. `cvsec` takes the pass count from 6.20 to 4.50.
  * Adding the negative off-diagonal parts: no gain, and a worse gap.
  * Refreshing the Jacobian every pass: same passes, +2.6 s.
  * Coloured probes: worse.
  * Warm start: not re-measured. Phase 3 rejected it because it perturbs non-converged
    calls, and with 0 non-converged calls and 4.5 passes (seed, 2 Newton passes,
    confirm) there is at most one pass left to win.
  * Colskip is already on. After the secant cv, 2232-3072 columns are still active on
    pass 3, so there is little left for it to skip.
* **(c) `ck_impl_tri`: 11.4 + 3.3 ms -> 1.4 ms per call** with the fused team kernel.
  * The serial Thomas out of LDS was enough; PCR was not needed.
  * Batching over g-points does not apply, because the matrix is summed over chains.
* **(d) Launch overhead and host syncs: small once `HSA_NO_SCRATCH_RECLAIM=1` is set.**
  * The ~160 us gaps before constant-memory launches in the T3 trace were the scratch
    re-attach, not the functor copy.
  * `lw` changes nothing (18.20 / 18.47 vs 18.26 / 18.24 s).
  * RT launch gaps are now ~20 us each, about 1.7 ms per T4 call.
  * The fused step has one read-back per pass instead of two.
* **(e) Kernels over 32 kB: none.**
  * `meta.v4.txt` (all 256 kernels of `dhj.v4.o`): 176 constant-memory launches, 80
    local, 0 global-memory.
  * No RT or implicit kernel reaches the > 32 kB path.

## 5. Files (`bench/impl_t4_0923/`)

* The patch: `t4.patch`.
* Build: `build.sh`, `src_new/` (build directories deleted: inode quota), `src_head/`,
  `src_md5.v4.txt`.
* Binaries:
  * `athena.gpu.v4`: final; timing and profiles.
  * `athena.gpu.v3`: jlin with private arrays; job 11944960/78.
  * `athena.gpu.v1`: `c02ddd34` base, with probes.
  * `athena.cpu.v5`: final; gates.
  * `athena.cpu.head`.
* Code-object tables: `counts2.v{1,3,4}.txt`, `meta.v*.txt`, `dis.v4.txt`.
* `gate/`: `gate.sh`, `cmpall.py`, `cmprel.py`, `rstcmp2.py`, the logs and the run
  directories (bin/rst deleted).
* `gpu/`:
  * scripts: `submit.sh`, `submit_prof.sh`, `arms.sh`, `ana4.py`, `rstdec.py`,
    `bincmp.py`, `prof.py`;
  * run directories `j1 j4 j5 j6 pf pf4 pf5 pf6 g1 h1 h10`.
* Logs: `log.out.<job>` in `gpu/`, and for jobs 11944895, 11944960, 11944978 and
  11945029 in the bench root.
* The v2 GPU binary was built while its source was being edited. It crashed and was
  deleted, and job 11944936 is void.
