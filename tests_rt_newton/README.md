# rt_newton on the deep-hot-Jupiter correlated-k run (2026-09-23)

Question: should `problem/rt_newton` (Newton refinement of the semi-implicit two-stream
apply, `src/utils/two_stream_rt.hpp` ~515 and ~5575) become DEFAULT ON?

## Verdict

| user | verdict | why |
| --- | --- | --- |
| **dhj (ck, prod4)** | **flipping `rt_newton` alone: NEUTRAL, bit for bit.** Flipping `rt_semi_lin = false` + `rt_newton = true`: **HURTS on the day side, helps at night: keep the legacy default** | see below |
| box_convection (He / B star) | **not affected** by any header-default change | the pgen reads all three switches with its own defaults (`box_convection.cpp:1826-1828`: semi_lin false, newton true, rescue_eq true) |
| red_giant | **not affected** | the pgen reads all three (`red_giant.cpp:1682-1684`: default = the mode-3 flag; `:1958-1960` ck path: legacy), and every production input sets them explicitly |

**Recommendation: do not flip anything for dhj.** Keep the legacy step
(`rt_semi_lin = true`, `rt_newton = false`, `rt_rescue_eq = false`) as the dhj default.
Also leave the header defaults alone. The only pgen that inherits them is
`deep_hot_jupiter_rt`, and there the flip is either inert (`rt_newton` alone) or measurably
worse on the day side (with `rt_semi_lin = false`). Pgens that want the new step already
select it themselves.

### Why rt_newton alone does nothing on dhj

* `deep_hot_jupiter_rt.cpp` (HEAD d837c355) never reads `rt_semi_lin`, `rt_newton` or
  `rt_rescue_eq`. An input key `problem/rt_newton = true` in a dhj input is ignored.
* Even if it is read, the Newton branch sits in the `else` of `if (semilin)`
  (two_stream_rt.hpp:5579), and `rt_semi_lin` defaults TRUE (:484). Under the legacy
  linearised step (`lambda = 4E/e`) the Newton code is never reached.
* Measured: arm **b** (semi_lin true, newton true) is **bitwise identical** to arm **a** in
  every variable of the final dump and in the .hst. The .bin files differ only in the
  embedded input text.

## Setup

* Binary: `bench/rtnewton_0923/athena.gpu`. It is `git archive` HEAD d837c355, plus the T1
  implicit patch (`bench/impl_t1_0922/t1.patch`), plus `bench/rtnewton_0923/pgen.patch`
  (three `GetOrAddBoolean` reads in the dhj pgen, defaults = the header defaults). Built
  with HIP gfx942 and MPI (`build.sh`). Nothing is committed.
* **Default-path gate:** arm **g** runs the T1 binary without the pgen reads. It is
  **bitwise identical** to arm **a** (all variables and the hst), so the patch is inert on
  the default path. Arm a_r2 is also bitwise equal to a, and N_r2 to N, so GPU runs are
  deterministic here.
* Job **11943604** (apudev, 1 node, 2 ranks / 2 GPUs, `submit.sh`). Restart
  `bench/cs_mhd_prod3/rst/dhj.00567.rst` (rot 283, cycle 4429831). Input `rtn.athinput` is
  `impl_t1_0922/gpu/impl.athinput` (the prod4 `<problem>` block plus diagnostics) with the
  three keys and `ck_impl_arat` added. Each arm runs 400 cycles (about 6.5e3 s, 0.02 rot).
  Arms were run in the order a, N, L, b, c, NR, a_r2, N_r2, g.

| arm | rt_semi_lin | rt_newton | rt_rescue_eq | other |
| --- | --- | --- | --- | --- |
| a (= prod4) | true | false | false | |
| L ("b" in the correction brief) | false | false | false | closed form `e_eq = e (A/E)^(1/4)` |
| N ("c") | false | true | false | |
| NR ("d") | false | true | true | |
| b | true | true | false | proves the flag is inert |
| c = reference | - | - | - | `ck_implicit = true`, `ck_impl_arat = 1e30`, maxit 8, reuse_jac 1, seed 2 |

## Numbers

Counters and timing come from `analysis/logs.py` (output in `analysis/logs.txt`) and the
`cpu time used` line in each run.log. The gap is the `rt_desum` rel,
`sum(de dx)/sum(src dt dx) - 1`, over every report line.

| arm | gap mean | gap max | rt_de_max clips | efix / rescue_eq / rescue_floor | eos_efloor | eos_tfloor | dt end [s] | cycles/s (loop) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| a | -9.43 % | -14.8 % | 0 | 0 / 0 / 0 | 1270943 | 171888 | 16.182 | 13.46 (a_r2 13.60) |
| L | -9.37 % | -14.7 % | 0 | 0 / 0 / 0 | 1275024 | 175600 | 16.182 | 13.43 |
| N | -3.91 % | -7.8 % | 0 | 0 / 0 / 0 | 1289355 | 198686 | 16.137 | 13.24 (N_r2 13.16) |
| NR | = N bitwise | | 0 | 0 / 0 / 0 | = N | = N | 16.137 | 11.88 (outlier, same arithmetic as N) |
| c | ck gap: max 1.8e-6, median 2.0e-7; 19 of 800 calls not converged | | - | - | 1302936 | 223165 | 15.691 | 2.97 |

* No `rt_de_max` clip and no positivity rescue fired in any arm. `rt_rescue_eq` therefore
  has nothing to act on, and NR is bitwise equal to N.
* Cost of Newton: 13.46 / 13.60 (a) against 13.24 / 13.16 (N) cycles/s, i.e. **+2.4 %
  loop time**. The implicit reference costs **4.5x**.
* A smaller rt_desum gap does not mean a more accurate answer. The gap only measures how
  far the applied `de` is from the explicit `src dt`. The accuracy metric is T(p) against c.

### T(p) against the implicit reference

`analysis/tp_vs_ref.py`, output in `analysis/tp_vs_ref.txt`. dT = arm - c, in kelvin.
Regions: day `mu > 0.1`, terminator `|mu| <= 0.1`, night `mu < -0.1`.

| p [bar] | region | a: mean / rms / max | L: mean / rms / max | N (= NR): mean / rms / max |
| --- | --- | --- | --- | --- |
| 1e-6 | day | -0.6 / 0.9 / 3.8 | -0.6 / 0.9 / 3.9 | -4.1 / 4.2 / 11.4 |
| 1e-6 | term | -2.2 / 3.0 / 9.5 | -2.3 / 3.1 / 10.1 | 4.1 / 9.9 / 30.5 |
| 1e-6 | night | 1.6 / 13.9 / 113 | 0.8 / 13.8 / 127 | 9.3 / 17.5 / 160 |
| 1e-5 | day | 8.8 / 9.9 / 26.0 | 9.8 / 11.1 / 29.4 | **21.8 / 26.4 / 71.7** |
| 1e-5 | term | 5.0 / 5.4 / 11.7 | 5.9 / 6.2 / 13.6 | **23.5 / 25.2 / 55.6** |
| 1e-5 | night | 2.6 / 20.9 / 331 | 2.1 / 19.3 / 331 | 10.0 / 17.3 / 148 |
| 1e-4 | day | -2.9 / 3.1 / 6.8 | -3.2 / 3.4 / 7.3 | **-19.1 / 21.7 / 36.8** |
| 1e-4 | term | -2.8 / 3.3 / 6.7 | -3.3 / 3.8 / 7.1 | 1.1 / 19.2 / 56.8 |
| 1e-4 | night | 15.1 / 30.8 / 226 | 13.9 / 28.8 / 226 | 15.4 / 20.0 / 130 |
| 1e-3 | day | -6.7 / 7.0 / 9.9 | -7.0 / 7.3 / 10.2 | **-33.6 / 36.4 / 54.0** |
| 1e-3 | term | -0.3 / 5.0 / 36.3 | -0.6 / 4.9 / 35.6 | -0.7 / 6.9 / 18.1 |
| 1e-3 | night | 38.3 / 44.0 / 176 | 36.0 / 41.5 / 174 | 18.6 / 22.0 / 207 |
| 1e-2 | day | 4.0 / 7.6 / 20.9 | 3.7 / 7.5 / 20.6 | 1.8 / 3.2 / 9.2 |
| 1e-2 | term | 25.4 / 25.7 / 32.9 | 25.1 / 25.3 / 32.4 | 13.8 / 14.2 / 21.4 |
| 1e-2 | night | 28.5 / 29.6 / 58.7 | 27.9 / 29.0 / 58.8 | 21.2 / 22.0 / 40.3 |

The file also records where each maximum sits (lat/lon). Two cases matter:

* The legacy arm's night-side maxima (113-331 K, 1e-6 to 1e-4 bar) sit at lat 9-15, lon
  -145 to -165, i.e. a few night-side columns. Whether they are the odd-even columns that
  README_T1_stall.md finds in the restart was not checked. N lowers the 1e-5 and 1e-4 bar maxima (331 -> 148 K, 226 -> 130 K)
  and raises those at 1e-6 and 1e-3 bar (113 -> 160 K, 176 -> 207 K).
* N's day-side errors are broad rather than one column: the rms is close to the mean.

Reading:
* **L is about the same as a.** Replacing the linearisation with the closed-form fixed
  point changes the rms by at most about 2 K.
* **Newton (N) is 3-7x worse than legacy on the day side at 1e-5 to 1e-3 bar**
  (rms 26 / 22 / 36 K against 10 / 3 / 7 K). It is also worse at the terminator at 1e-5 and
  1e-4 bar.
* **Newton is better at night at 1e-4 to 1e-2 bar and at the deep terminator**
  (rms 20 / 22 / 22 K against 31 / 44 / 30 K).
* Newton moves the floor counters toward the reference (tfloor 171888 for a, 198686 for N,
  223165 for c).

Caveat: this is a 400-cycle (0.02 rot) window after switching the scheme on a legacy
restart, so it measures the transient departure from the converged solve, not a
long-term climate.

## Why Newton hurts on the day side: the T^4 emission scaling

The Newton residual is `F(de) = de - A dt + E_old (T_new/T_old)^4 dt` (:5586). Newton now
takes T(e) and c_v from the EOS, but the cell's emission still scales as T^4. At the
day-side top (README_T1_stall.md) the true slopes are:

* d ln T/d ln e = 5.4
* d ln E/d ln T = 9.8 (Wien-side bands)
* so d ln E/d ln e = 53

The three schemes assume the following d ln E/d ln e:

* legacy: 4 (`lambda = 4E/e`)
* Newton: 4 x 5.4 = 22
* truth: 53

Newton therefore gets the heat capacity right while its fixed point is still
`E ~ T^4 = A`. It relaxes toward an equilibrium T that differs from the true one,
`T (A/E)^(1/9.8)`, by a factor 2.45 in ln T. The legacy step never targets a fixed point,
since it is only a damped explicit step, so it suffers less from the wrong exponent.
Making the EOS "more correct" in one factor while the emission slope stays wrong makes
the day side worse.

**What a band-summed E(T) would need.** Where the sweep kernel forms `Em_g` per band, it
would also form the band-summed logarithmic slope
`s = d ln sum_b(kappa_P,b(T) rho B_b(T)) / d ln T`. That can be done analytically from
dB_b/dT and the table's d kappa_b/dT, or by one extra evaluation at T(1+eps). Store `s` in
one extra (m,k,j,i) array, and in the Newton use `r = (T1/T0)^s` and `dr/de = s r/T1 dT/de`
in place of the 4. This is one array plus one more Planck/opacity evaluation per band and
cell in the sweep, so the per-iteration cost stays unchanged. Re-evaluating the band sums
exactly at each Newton iterate would cost `nband x ng` table lookups per iteration and is
not needed. This is untested. It is the change that would make `rt_newton` physically
consistent on dhj, and it should be gated against arm c as done here.

## What flipping the defaults would change bitwise

* **Header `rt_newton` default true (alone):** no change in any pgen or input.
  * box_convection and red_giant read it with their own defaults.
  * dhj never reaches the branch while `rt_semi_lin = true`, as arm b shows.
* **Header `rt_semi_lin` default false (with or without `rt_newton`):** changes every
  `deep_hot_jupiter_rt` run, and nothing else. None of these dhj inputs sets it:
  * `inputs/production/deep_hot_jupiter_cs_prod4.athinput`
  * `inputs/production/deep_hot_jupiter_cs_hyd4.athinput`
  * `inputs/tests/dhj_ck_implicit.athinput`
  * `inputs/tests/dhj_ck_spherical.athinput`
  * `inputs/mhd/deep_hot_jupiter_rt_{eos,ideal_xe,blowup}.athinput`
  * the regression tests built with `PROBLEM=deep_hot_jupiter_rt`:
    `tst/test_suite/rad/test_rad_dhj_ck_{cpu,mpicpu}.py`,
    `test_rad_dhj_srclim_cpu.py` (via `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`) and
    `tst/test_suite/nr/test_nr_geneos_repro_gpu.py`

  These tests were not run here. They check physical tolerances and run-to-run
  reproducibility rather than stored dumps, so whether they would still pass is untested.
  Also affected: the running prod4 production, if it were rebuilt.
* A red-giant input that omits the switches keeps its pgen default. The memory note from
  09-11 records that the six opt-ins must be set in every red-giant input from pin10 on.

## Box / red-giant evidence (from memory notes, nothing run)

* `rt-source-dt-forcing.md` (09-14, He box `hestar_fecz/conv_fix`, clip experiment): "rt_semi_lin / rt_newton off shift 10 %, same
  dt spread". These are not the pump, and the box already uses the new step by default.
* `red-giant-explicit-conduction-explosion.md`: `rt_newton=false stalls 1.51497e6
  (T15's time)`, i.e. the red giant's T15 failure is independent of rt_newton.
* two_stream_rt.hpp ~5600 (RG_v4): the bare Newton diverged in thin radiation-dominated
  red-giant cells and was safeguarded with a bisection. The red giant runs Newton ON by
  input.
* `red-giant-rt-bface-switch-and-head-regression.md` (09-11): these defaults were set to
  the OLD step on purpose, so that ck dhj stays bit-identical to HEAD.

## Files

`/viper/u2/jinma/ATHENAK/bench/rtnewton_0923/`:
* `pgen.patch`, `build.sh`, `athena.gpu`, `HEAD.txt`
* `rtn.athinput`, `submit.sh`, `log.out.11943604`
* arm directories `a L N NR b c a_r2 N_r2 g`
* `analysis/{logs.py,logs.txt,tp_vs_ref.py,tp_vs_ref.txt}`
