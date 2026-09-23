# Newton globalisation for the implicit correlated-k column solve (T4), 2026-09-24

Switch: `problem/ck_impl_glob = none | ls | ls_sub` (default `none`), branch `ck-glob`.
Code: `src/utils/two_stream_column_ck.hpp` (fused step, `CkImplStep`), the pass loop in
`src/utils/two_stream_rt.hpp`, key parsing in `src/pgen/deep_hot_jupiter_rt.cpp`. It needs
`ck_impl_fuse` (fatal otherwise). Extra keys: `ck_impl_ls_ntry` (4), `ck_impl_ls_c` (1e-4),
`ck_impl_sub_max` (8).

## What it does

* **ls**: a backtracking line search for each column. The merit is
  phi = ||R_i/(e*_i + eps e*_max)||_2, with e* = e^n fixed over the call.
  * The first trial is T4's capped step (ck_impl_dtmax / demax), unchanged.
  * The next pass's sweep evaluates the trial. If phi > (1 - c alpha) phi_old, the column
    is pulled back to half its trial and takes no Newton step that pass. There is no new
    Jacobian and no extra sweep for the other columns. After `ls_ntry` halvings the trial
    is accepted.
  * Seed steps and fallback steps are not line-searched.
* **ls_sub**: ls, plus sub-stepping for each column.
  * A column still above tol after `maxit` residual evaluations of its current (sub-)step
    is put back to e^n (u0 -= dep). It is then redone as 2, 4, ... `sub_max` backward-Euler
    sub-steps of bdt/2^L, each with its own maxit budget.
  * Each sub-step starts with the `ck_impl_seed = 2` guess.
  * The pass after a sub-step seed rebuilds the chord Jacobian for every live column
    (`ck_impl_jac_again`). Without that rebuild, and without the sub-step seed, sub-steps
    as short as 6 s did not converge (the chord J was built at the whole step's seeded
    state).
  * The source of finished sub-steps is carried in `ck_sacc`. The ckdesum gap therefore
    compares sum (e - e^n) dx with sum (sacc + h S) dx.
  * A column that fails at level sub_max is left at its last iterate and counted
    (`subfail`, NOT-CONVERGED).
  * Pass cap: maxit (1 + 2 + .. + sub_max) + 2 (1 + .. ) + 2 = 152 at the defaults.
* Everything is per column (one team per column), with no host sync per trial. The
  counters are atomics on the diagnostic array only, so the state is deterministic.
* Report line: `lsrej`, `subrst` and `subfail` (totals over the call).

## Gates

The CPU gates use the well-posed suite (wellposed/ plus wellposed_hooks.patch, which is
not committed). Binaries: `athena.cpu.base` = HEAD d80c84ca + hooks, and `athena.cpu.glob`
= ck-glob + hooks. Runs are under `/viper/ptmp2/jinma/ckglob_0924`. Scripts are in
`tests_ck_implicit/glob/`: `run_glob.sh G1|A1|A2|D|Dref|E`, `ana_glob.py`, `ana_gpu.py` and
`submit_gpu.sh`. Arms: `t4` (as wellposed/arms.sh, tol 1e-8, maxit 8), `ls`, `sub`
(ls_sub, sub_max 8) and `s32` (ls_sub, sub_max 32).

| gate | result |
| --- | --- |
| 1 none bitwise | **PASS.** CPU: t4 and semi at dt 20 / 2000, and D-H2 at 2000 (G1b.log): .hst, bin and rst are bitwise after the parameter text, and rec.txt is identical. GPU, production hyd4 restart, 150 cycles: t4 on glob vs base binary, 3 files BITWISE (job 11953105). |
| 2 A1 / A2 | A1: sub removes every non-converged call at 2000 s (200 -> 0), with an identical state (max 7.076e-3 vs 7.075e-3 against T**). ls changes nothing. A2 at 2000 s (one call from the IC): error 0.53 -> 0.086 (sub), still 1 non-converged. |
| 3 D, H2, 200 / 2000 s | 200 s: **PASS** for ls_sub (0/100 non-converged, T 3229-4406 K, error vs the dt=20 reference 1.4e-4). 2000 s: **FAIL at the briefed sub_max = 8** (100/100 non-converged, T collapses to 32 K). PASS with sub_max = 32 (2/100, T 3114-4406 K, error 4.8e-3 / 7.4e-4). |
| 4 E symmetry, 2000 s | **PASS** for ls_sub. Column spread 2.3 / 18.7 (t4, tm / tmb) -> 8.5e-11 / 5.5e-11 (sub) and 2.2e-12 / 1.5e-12 (s32). ls alone does not help (23 / 18). |
| 5 20 s unchanged | **PASS.** rec.txt identical for t4 / ls / sub / s32 in A1, A2, D and E at 20 s (D t4 vs sub: all output files bitwise). No trial was rejected, so the passes are the same. A1 at 200 s is identical too. |
| 6 GPU cost | ls and ls_sub cost the same as T4 at the production dt. They are bitwise T4 there (no rejection, no sub-step). |

### A1 (from T** = the final restart of wellposed_0923 a1_t4x_20, 200 calls; the A0b restart is gone)

```
arm     dt    max@200   mean@200   passes  maxp nonconv  resfinal    |bud0|  lsrej subrst  sfail
t4      20  2.033e-04  4.693e-05     3.00     3       0  2.83e-10  6.50e-05      0      0      0
ls/sub/s32 20: identical to t4
t4     200  7.202e-04  9.930e-05     4.00     5       0  5.99e-09  7.13e-05      0      0      0
ls/sub/s32 200: identical to t4
t4    2000  7.075e-03  6.108e-04     8.00     8     200  4.48e-08  5.69e-05      0      0      0
ls    2000  7.075e-03  6.108e-04     8.00     8     200  4.48e-08  5.69e-05      0      0      0
sub   2000  7.076e-03  6.108e-04    18.09    35       0  3.60e-08  3.71e-05      0   6432      0
s32   2000  (= sub)
```

At 2000 s, T4 crawls linearly: the residual is 4.5e-8 at maxit against tol 1e-8. ls_sub
restarts every column at level 1 on every call (subrst 6432 = 32 per call). That doubles the
passes (18.1 vs 8) for the same state.

### A2 (vs t4x dt=4, t = 2000 s)

```
arm     dt        max       mean   passes  maxp nonconv  lsrej subrst  sfail
t4/ls/sub/s32 20: 1.176e-04 1.980e-05, 4.02 passes, 0 nonconv (identical)
t4     200  2.267e-03  3.433e-04     7.10     8       4      0      0      0
ls     200  2.151e-03  3.301e-04     7.30     8       4    128      0      0
sub    200  1.730e-03  2.487e-04    16.90    64       0     64    192      0
t4    2000  5.302e-01  3.811e-02     8.00     8       1      0      0      0
ls    2000  5.137e-01  2.963e-02     8.00     8       1     32      0      0
sub   2000  8.573e-02  2.434e-02    37.00    37       1    224     96     32
s32   2000  1.324e-01  4.130e-02    55.00    55       1    224    160     32
```

### D (eos_h2 true, mu0 1, 100 calls). err@t1 / err@t2: max rel. T error vs t4 at dt 20 at the same time (calls 10 / 100 of the arm)

```
arm     dt   passes  maxp nonconv     err@t1     err@t2  Trange@100  lsrej subrst  sfail
t4/ls/sub/s32 20: 4.10 passes, max 6, 0 nonconv, T 3312-4407 (identical)
t4     200     8.00     8     100  5.669e-01  9.836e-01    72- 4745      0      0      0
ls     200     6.18     8       9  2.263e-03  1.414e-04  3229- 4406    448      0      0
sub    200     7.52    69       0  4.715e-04  1.354e-04  3229- 4406     64    256      0
t4    2000     8.00     8     100  6.045e-01  9.835e-01    72- 4994      0      0      0
ls    2000     8.00     8     100  6.133e-01  9.928e-01    32- 4740   6968      0      0
sub   2000    44.16    47     100  6.322e-01  9.928e-01    32- 4566  41144   9600   3200
s32   2000    19.66   193       2  4.781e-03  7.442e-04  3114- 4406    384   1856     64
```

* ckdesum gap, |max| / median over the calls:
  * D 200: t4 3.3e-2 / 4.3e-3, sub 7.6e-6 / 5.2e-7.
  * D 2000: s32 2.9e-4 / 5.7e-7.
  * A1 2000: t4 2.9e-5 / 3.3e-6, sub 2.8e-6 / 1.1e-6.
* Note: the rec.txt `resfinal` / `bud0` measure the WHOLE-step backward-Euler residual of
  the final state. A sub-stepped column does not satisfy it by construction. For ls_sub,
  the convergence measures are the solver's own sub-step residual and the gap above.

### Why ls alone does not converge: the diagnosis (per-pass history, ck_impl_debug=-2, G1/gdh_*)

* At 200 s with H2, T4 converges only LINEARLY, even with a fresh Jacobian on every pass
  (reuse_jac = 0).
  * Rate per pass: ~0.87 (chord) and ~0.6 (fresh).
  * Passes to tol 1e-8: 25-30 (chord) and 12-23 (fresh).
* Linear convergence means the tridiagonal Jacobian is inexact. It truncates the
  non-local two-stream coupling and drops the negative per-chain parts. The contraction
  factor grows with the step length h.
* At 2000 s the iteration enters a period-3 limit cycle, with steps alternating at the
  0.5 dtmax cap and the 0.333 demax bound.
* Backtracking breaks the cycle, but it cannot speed up a linear rate. The inexact
  direction is often not a descent direction either (thousands of rejections per run at
  2000 s).
* Shortening h (sub-steps) is what restores convergence.
* In the D-H2 transient at 2000 s, 250-s sub-steps are still too long. 62.5-s sub-steps
  (sub_max 32) work.

## GPU cost (job 11953105; apudev; 2 ranks / 2 GPUs; HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1; 150 cycles)

* Restart: `bench/cs_hyd4_prod/rst/dhj.00113.rst`, read in place. Rotation 56.50,
  t = 1.72325e7 s, ncycle 873580, dt 18.7 s, md5 353655f0...
* Input: hyd4 production input plus the T4 `<problem>` block and the glob keys
  (`gpu/hyd.athinput`).
* Binary: `athena.gpu.glob` (md5 ec58abdb...). The last run uses the base binary
  (dc3b3d5c...).
* Arm order: s t4 ls sub | sub ls t4 s | s t4 ls sub. `cpu time used`.

```
arm            cpu r1/r2/r3 [s]  x semi   passes  maxp nonconv  lsrej subrst  sfail
s               10.00/9.96/9.91    1.00
t4            17.60/18.52/17.39    1.79     4.50     5       0      0      0      0
ls            17.21/17.59/17.45    1.75     4.50     5       0      0      0      0
sub           17.26/17.37/17.66    1.75     4.50     5       0      0      0      0
```

* ls_r1 and sub_r1 are bitwise t4_r1 (3 files). So is t4base_r1 on the HEAD binary, and
  so is t4_r2.
* At the production dt the globalisation is never triggered and costs nothing measurable.
  The ~2 % difference is within the repeat spread.
* A first submission (11952750) and a debug job (11952869) crashed with a GPU memory
  access fault. The cause was the glob GPU binary: its pgen object had been compiled from
  a source snapshot that was modified in the middle of the compile. A clean rebuild (touch
  + make) fixed it. Those jobs are void.

## Open

* **Gate 3 at 2000 s fails with the briefed sub_max = 8.** sub_max = 32 passes, at 19.7
  mean passes per call (maxp 193) in the transient.
* **The escalation discards progress.** A column that is converging linearly and slowly
  (A1 at 2000 s: residual 4.5e-8 at maxit) is restarted from e^n, which doubles its
  passes.
  * Cheaper alternative 1: restart only the failing sub-step from its start state.
  * Cheaper alternative 2: grant extra passes while phi keeps contracting.
* **The root cause is the inexact tridiagonal Jacobian.** Any improvement there (non-local
  coupling, keeping the negative parts) acts on the rate itself.
