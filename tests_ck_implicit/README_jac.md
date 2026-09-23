# Cheaper escalation and a better Jacobian for the implicit ck column (T4), 2026-09-23

Branch `ck-jac` (from rt-integration be02c647). Code: `src/utils/two_stream_column_ck.hpp`
(the fused step `CkImplStep`), the pass loop in `src/utils/two_stream_rt.hpp`, key parsing
in `src/pgen/deep_hot_jupiter_rt.cpp`. Scripts: `tests_ck_implicit/jac/`. Both switches
default off, and off is bitwise (gate 1).

| key | default | meaning |
| --- | --- | --- |
| `ck_impl_esc` | 0 | 1 = cheaper escalation of `ck_impl_glob = ls_sub` (fatal without it) |
| `ck_impl_esc_rho` | 0.9 | largest per-pass merit ratio that counts as "still converging" |
| `ck_impl_esc_extra` | 8 | extra residual evaluations a converging (sub-)step may take |
| `ck_impl_aa` | 0 | Anderson depth m (0 = off, at most 8; needs `ck_impl_fuse`) |
| `ck_impl_aa_rst` | true | restart the Anderson history when the weighted residual grew |

## Deliverable 1: ck_impl_esc (ls_sub escalation)

* **Extra passes.** A (sub-)step that has used its `maxit` residual evaluations normally
  escalates. It keeps iterating instead when both of these hold:
  * the merit fell by `esc_rho` or better over the last pass;
  * the passes to tol projected at that rate, log(tol/r)/log(rate), fit in
    `maxit + esc_extra`.
* **Local restart.** A sub-step that still fails is redone on its own, from its start state
  e^n + `ck_sacc`, as two sub-steps of half its length. `dep` is reset to `sacc`. ls_sub, by
  contrast, puts the whole column back to e^n and redoes every sub-step at the next level.
* **Coarsening.** A sub-step that converged in <= maxit/2 steps hands the column back to the
  coarser level, where the finished sub-steps end on that level's grid.
* **Bookkeeping.** Level and index stay in `ck_lsc` (slots 3 and 4). The pass cap grows by
  `esc_extra` x (1 + 2 + .. + sub_max). The report line adds `escx` (extra column-passes
  granted) and `coarse`.

## Deliverable 2: ck_impl_aa (Anderson / DIIS acceleration: a Jacobian-free Krylov step)

* **Why a solver change, not a better matrix.** The tridiagonal chord J leaves out two
  things: the non-local two-stream coupling and the negative per-chain parts. The dense
  exact J is not built. Instead, each column keeps the last m differences of its iterates
  (dX) and of its backward-Euler residuals (dR). The residual comes from the ordinary sweep
  every pass already runs, so there is **no extra sweep**.
* **The step.**
  * gamma = argmin || W (r - dR gamma) ||_2, where W is the merit weights
    1/(e*_i + eps e*_max).
  * e_new = e - dX gamma + P^{-1}(-(r - dR gamma)), where P is the chord tridiagonal of
    that pass.
* **What this is.** On a linear problem with a fixed P this is GMRES preconditioned by the
  tridiagonal. The Jacobian's missing parts are sampled through the sweep: the sweep is
  the operator and the tridiagonal is the preconditioner. Because the least squares is on
  the true residual, a P that changes from pass to pass (`cvsec`) does not invalidate the
  history.
* **Least squares.** It uses a regularised Cholesky of the m x m Gram matrix (lambda =
  1e-10 max diag) on one lane, in team scratch. It is skipped if the matrix is not
  positive definite.
* **Thin rows.** Rows whose step is the per-cell thin solve are neither in the least squares
  nor mixed.
* **History restarts:**
  * pass 0;
  * every new sub-step (a new backward-Euler residual);
  * a fallback step;
  * whenever the weighted residual norm grew (`aa_rst`).
  A seed step's pair is never stored.
* **Report line:** `aa` counts the accelerated column-steps.
* **The negative per-chain parts (`ck_impl_jneg`, already in HEAD) are NOT safe** (A2 at
  200 s, dt 200, 10 calls):
  * jn alone: 10/10 non-converged, against 4/10 for t4.
  * jn + aa: 3/10, the same as aa alone.
  This is left out.

## Gates (CPU: well-posed suite + wellposed_hooks.patch; runs /viper/ptmp2/jinma/ckjac_0923)

Binaries: `athena.cpu.base` = be02c647 + hooks, and `athena.cpu.jac` = ck-jac + hooks.
Scripts: `jac/run_jac.sh G1|A1|A2|D|E [arms]` and `jac/ana_jac.py A1|A2|D|E` (tables), plus
`ana_jac.py H <run> <call>` for per-pass histories. The input is `wpj.athinput`, which is
wp.athinput with the glob, esc and aa keys (`jac_keys.athinput`).

The arms all use T4 with tol 1e-8, maxit 8 and ck_impl_debug -2:

| arm | switches |
| --- | --- |
| t4 | none added |
| sub | ls_sub, sub_max 32 (the current default) |
| esc | sub + esc 1 |
| aa | t4 + aa 4 |
| aas | sub + aa 4 |
| aae | sub + esc 1 + aa 4 |

| gate | result |
| --- | --- |
| 1 off bitwise | **PASS.** CPU: t4, sub and semi at 20 / 2000 s, plus D-H2 at 2000 s (G1): base vs jac binary, 5 files BITWISE and rec.txt identical, all 9 pairs. GPU: see below. |
| 2 dt 20 s | **PASS.** esc: identical to t4 (A1, A2, D, E). aa / aas / aae: within tolerance with the same passes: A1 3.00, A2 4.02, D 4.10, E 4.05 / 4.10. A2 error vs t4x 1.175e-4 vs 1.176e-4. D T vs t4 at call 100: 3.1e-8. E column spread 8.9e-13 vs 9.8e-13. |
| 3 A1 / D / E at 200-2000 s | Passes per call fall at every point, and non-converged falls to 0 wherever aae runs (tables). |
| 4 A2 at 2000 s (one call from the IC) | Accuracy up 2.5x, passes up 2.2x (see "Cost of the finest level" below). |

### A1 (from T** = wellposed_0923 a1_t4x_20 rst 3, 200 calls; errors vs T**)

```
arm     dt    max@200   mean@200   passes  maxp nonconv  resfinal    |bud0| subrst   escx     aa
t4..aae 20: identical state 2.033e-04 / 4.693e-05, 3.00 passes, 0 nonconv
t4..aae 200: 7.202e-04 / 9.930e-05, 4.00 passes, 0 nonconv (aa rows: resfinal 5.79e-9 vs 5.99e-9)
t4     2000  7.075e-03  6.108e-04     8.00     8     200  4.48e-08  5.69e-05      0      0      0
sub    2000  7.076e-03  6.108e-04    18.09    35       0  3.60e-08  3.71e-05   6432      0      0
esc    2000  7.075e-03  6.108e-04    14.71    23       0  3.87e-08  3.70e-05    464  19856      0
aa     2000  7.075e-03  6.108e-04     8.00     8     196  1.94e-08  8.41e-06      0      0  37000
aas    2000  7.075e-03  6.108e-04    10.05    18       0  3.84e-08  3.78e-05    208      0  37736
aae    2000  7.075e-03  6.108e-04     9.48    18       0  3.85e-08  3.78e-05     72    256  37488
```

The A1 cost case of README_glob (ls_sub 18.1 passes/call) falls to 14.7 (esc) and to
9.5 (aae), with the same state and 0 non-converged.

### D (eos_h2 true, mu0 1, 100 calls; err vs t4 at dt 20 at the same t)

```
arm     dt   passes  maxp nonconv  resfinal     err@t1     err@t2  Trange@100 subrst  sfail   escx     aa
t4      20     4.10     6       0  9.69e-09        nan  0.000e+00  3312- 4407      0      0      0      0
aa/aas/aae 20: 4.10 passes, max 6, 0 nonconv, err@t2 3.143e-08
t4     200     8.00     8     100  4.56e-03  5.669e-01  9.836e-01    72- 4745      0      0      0      0
sub    200     7.52    69       0  2.73e-04  4.715e-04  1.354e-04  3229- 4406    256      0      0      0
esc    200     6.96    58       0  2.30e-04  6.656e-04  1.382e-04  3229- 4406     96      0    672      0
aa     200     5.26     8       5  1.85e-05  8.267e-04  1.405e-04  3229- 4406      0      0      0   7392
aas    200     5.77    43       0  2.30e-04  6.661e-04  1.382e-04  3229- 4406     96      0      0   8320
aae    200     5.56    31       0  1.48e-04  7.632e-04  1.395e-04  3229- 4406     32      0    320   8032
t4    2000     8.00     8     100  7.75e-02  6.045e-01  9.835e-01    72- 4994      0      0      0      0
sub   2000    19.66   193       2  3.53e-02  4.781e-03  7.442e-04  3114- 4406   1856     64      0      0
esc   2000    12.65   179       1  3.53e-02  3.152e-03  4.103e-04  3114- 4406    448     32   2304      0
aa    2000     8.00     8     100  8.12e-02  6.043e-01  6.956e-01  1333- 5654      0      0      0  14464
aas   2000    11.74   201       0  9.40e-03  1.157e-03  6.673e-05  3114- 4406    448      0      0  22176
aae   2000     9.43   135       0  8.98e-03  1.325e-03  6.723e-05  3114- 4406    224      0    480  18400
```

D at 2000 s is the transient where ls_sub needed sub_max 32. Compared with ls_sub, aae
has:

* 2.1x fewer passes (19.7 -> 9.4);
* no non-converged columns (2 -> 0);
* 11x smaller error against the dt = 20 s reference at call 100 (7.4e-4 -> 6.7e-5).

Plain aa (no sub-steps) cannot rescue the H2 transient at 200 or 2000 s: the step must
still be shortened there.

### E (symmetry, identical columns, 20 calls; column spread max rel T)

```
run             passes nonconv  spread    | run              passes nonconv  spread
e_tm_2000_t4      8.00      20  2.35e+00  | e_tmb_2000_t4      8.00      20  5.47e+00
e_tm_2000_sub    41.95       1  2.21e-12  | e_tmb_2000_sub    45.00       1  1.51e-12
e_tm_2000_esc    21.95       0  1.53e-12  | e_tmb_2000_esc    22.65       0  9.58e-13
e_tm_2000_aa      8.00      20  2.43e+00  | e_tmb_2000_aa      8.00      20  3.70e+00
e_tm_2000_aas    23.70       0  2.17e-12  | e_tmb_2000_aas    24.70       0  1.21e-12
e_tm_2000_aae    14.85       0  1.96e-12  | e_tmb_2000_aae    15.10       0  1.55e-12
dt 20: every arm 4.05 (tm) / 4.10 (tmb) passes, 0 nonconv, spread <= 9.8e-13
```

### A2 (from the IC to t = 2000 s, vs t4x dt = 4)

```
arm     dt        max       mean   passes  maxp nonconv subrst  sfail   escx     aa
t4/sub/esc 20: 1.176e-04 1.980e-05 4.02; aa/aas/aae 20: 1.175e-04 1.979e-05 4.02
t4     200  2.267e-03  3.433e-04     7.10     8       4      0      0      0      0
sub    200  1.730e-03  2.487e-04    16.90    64       0    192      0      0      0
esc    200  2.070e-03  3.068e-04    12.80    51       0     64      0    448      0
aa     200  2.279e-03  3.443e-04     7.10     8       3      0      0      0   1408
aas    200  1.953e-03  2.860e-04    12.00    42       0     96      0      0   2336
aae    200  2.162e-03  3.217e-04     9.20    27       0     32      0     64   1824
t4    2000  5.302e-01  3.811e-02     8.00     8       1      0      0      0      0
sub   2000  1.324e-01  4.130e-02    55.00    55       1    160     32      0      0
esc   2000  5.307e-02  1.004e-02   174.00   174       0    160      0    256      0
aa    2000  4.977e-01  3.047e-02     8.00     8       1      0      0      0    128
aas   2000  5.306e-02  1.006e-02   188.00   188       0    160      0      0   3296
aae   2000  5.308e-02  1.001e-02   123.00   123       0    128      0    128   2464
```

**Cost of the finest level.** A2 at 2000 s is a single call from the IC. ls_sub gives up at
the finest level: 55 passes, 1 non-converged, error 0.13. esc / aae finish the column
instead: 0 non-converged, error 0.053, but 123-188 passes, because most 62.5-s sub-steps of
the call run at the finest level. That is the price of converging there.

### Per-pass residual history (max rel. residual, ck_impl_debug = -2; `ana_jac.py H`)

```
A1 2000 t4   call 100: 1.49e-04 1.48e-04 7.81e-06 3.50e-07 2.25e-07 1.42e-07 8.82e-08 5.46e-08
A1 2000 aa   call 100: 1.49e-04 1.48e-04 7.81e-06 3.46e-07 2.11e-07 1.21e-07 4.03e-08 1.16e-08
D  200  sub  call 50:  8.11e-04 7.24e-04 1.87e-05 5.80e-07 1.94e-08 7.27e-10
D  200  aae  call 50:  8.11e-04 7.24e-04 1.87e-05 2.57e-07 2.34e-09
D  200  t4   call 50:  1.78e-03 1.90e-03 2.83e-03 8.82e-04 4.40e-04 1.74e-03 1.50e-03 3.99e-04
A2 200  t4   call 2:   6.14e-03 5.01e-03 2.53e-04 1.03e-05 5.35e-06 4.49e-06 3.65e-06 3.06e-06
A2 200  aa   call 2:   6.14e-03 5.01e-03 2.53e-04 1.32e-06 4.04e-07 2.35e-07 6.51e-08 1.46e-08
```

* **Superlinear only where the column is near linear.**
  * D at 200 s: the ratio per pass after pass 3 is 0.014, then 0.009, against 0.03 for the
    chord.
  * A2 at 200 s, the same call: 1.5e-8 against 3.1e-6 after 8 passes (200x).
* **A1 at 2000 s the rate improves but stays linear:** 0.6 per pass (t4) against 0.3-0.5
  (aa). m = 4 secant directions do not span the error of the dense coupling within maxit.
  aa = 8 and aa_rst = false gave the same A2 result as aa = 4 (sweep in X/).

## GPU cost (production restart; apudev; 2 ranks / 2 GPUs; HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1; 150 cycles)

**PENDING at commit time.** Jobs 11953388 (SET=A: s t4 aa aae, plus t4 with the base
binary) and 11953389 (SET=B: s t4 sub esc, plus sub with the base binary) were submitted
with `jac/submit_gpu.sh`. Details:

* Restart: `bench/cs_hyd4_prod/rst/dhj.00127.rst`, read in place. Rotation 63.50,
  t = 1.93675e7 s, ncycle 986896, dt 19.47 s, md5 9fb8a77f...
* Input: `gpu/hyd.athinput`, which is the hyd4 production input plus the T4 `<problem>`
  block of README_glob plus `jac_keys.athinput`.
* Binaries: `athena.gpu.jac` and `athena.gpu.base` (be02c647 + hooks).
* Arm order: interleaved r1 / reversed r2 / r3.

Table: `python3 tests_ck_implicit/jac/ana_gpu.py`. GPU gate 1: the `cmpdir.py` line at the
end of each log (`/viper/ptmp2/jinma/ckjac_0923/gpu/log.out.<job>`) must read BITWISE.

## Open

* The finest-level cost in the one-call transient (A2 at 2000 s) could be reduced by a
  coarser `sub_max` with esc, or by coarsening sooner. Neither was tried.
* The dense two-stream Jacobian (built from the stored factorisation) would be needed for
  true Newton at 2000 s in the steady state (A1). Anderson recovers only part of it.
