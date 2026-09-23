# ck-fast: five cost levers for the implicit correlated-k RT (T4), 2026-09-23/24

Branch `ck-fast` (worktree /viper/ptmp2/jinma/wt_ckfast, from rt-integration 76899483).
The goal was to bring the T4 RT down to the cost of the hydro. The starting point is
tests_ck_implicit/prof_0923/README.md: on 2 GPUs, hydro takes 15.4 ms/cycle and T4 RT
takes 110 ms/cycle.

Every lever has its own `<problem>` switch. All default off. With all of them off, the
code is bitwise the code without them (section 2). Runs are in
/viper/ptmp2/jinma/ckfast_0923. Scripts are in `tests_ck_implicit/fast/` (GPU job
scripts in `fast/gpu/`).

**What "T4" means here.** It is the prof_0923 arm t8: ck_implicit, arat 1e30, frozen_op,
lin, lin_thr 1, fuse, jac_lin, cvsec, tol = dtol = 1e-8, maxit 8. It does not use
reuse_jac or seed. The well-posed README's "t4" also had reuse_jac 1 and seed 2.

**Result.** The best measured combination reaches **RT = 1.6x hydro** (c8: 40.6 ms/cycle
against 15.3 for hydro). That is 4.4x cheaper than T4, which is 7.2x hydro. **The target
of 1x hydro is not reached.** Section 7 says why and what is left.

## 1. The switches

| lever | switch(es) | what it does | needs |
|---|---|---|---|
| 1 | `ck_impl_once` (existed) | Runs the whole ck RT once per cycle, after the last RK stage, with the full dt. It works unchanged with T4/lin/fuse; no fix was needed. | ck_implicit, not rt_rad_force |
| 2 | `ck_impl_xstep = k`, `ck_impl_xstep_thr` | Keeps the frozen_op operator across calls: kappa rho, the coefficient triple, the tm factorisation, the cut and the frozen beam Qb. A call that re-applies it starts with a frozen linear pass: no cut, opacity, chain sweep or lin_build. A new operator is stored when none exists yet, when k cycles have passed since the last store, or when thr > 0 and some ck cell has moved in T or rho by more than thr since the store. For the last test, T is measured after rt_pre_tp of pass 0. | frozen_op |
| 3 | `ck_impl_jreuse = rho`, `ck_impl_jreuse_act` (0.25), `ck_impl_jreuse_xc` | Chord Jacobian. It is built on pass 0 and rebuilt after a pass in which the max residual contracted by less than rho, or in which fewer than act x columns were still active. xc carries the Jacobian across calls that re-apply an operator. **Do not use xc** (section 4). | reuse_jac 0, glob none |
| 4 | `ck_impl_pred`, `ck_impl_pred_fac` (0.5), `ck_impl_pred_chk`; `ck_impl_warm` + `ck_impl_warm_step`; `ck_impl_cvkeep` | pred removes the confirmation pass. A column is marked converged right after its step when r_p^2/r_{p-1} <= fac x tol and no cap was hit. pred_chk is a diagnostic only: one extra sweep gives the final state's own residual and gap. warm seeds each call with the previous call's increment; warm_step forces the pass-0 step. cvkeep builds the pass-0 rows with the previous call's secant cv. | pred: fuse; cvkeep: cvsec |
| 5 | `ck_impl_nosync` | Removes all per-pass device allocations. Each 1-element dummy View used to cost a hipMallocAsync, a zero fill and 2 hipDeviceSynchronize; they are now cached (`CkDum<V>`). The host mirror is cached. Scalar fills are stream-ordered. The rt_apply clip count goes to a device scalar. The clamp and efix read-backs and rt_pre_geom run once per call. rt_apply returns early under ck_implicit when no diagnostic is on. | - |

## 2. Gates: bitwise when off (CPU; GPU for nosync)

CPU binaries were built from git archive snapshots:
- `athena.cpu.base`: 76899483 plus `wellposed/wellposed_hooks.patch`.
- `athena.cpu.fast6`: 2b3e6917.

The same hooks are committed on ck-fast (be70b99b).

| comparison | result |
|---|---|
| base vs fast (all levers off), well-posed A2 (100 calls), arms t4 and semi | rec.txt bitwise; rst data bitwise (the parameter header differs only by the new keys) |
| base vs fast, full hydro with the ck RT per RK stage (6 cycles from the wp IC), T4 | rst data and .hst bitwise |
| t4 vs t4 + nosync, CPU (A2 and full hydro) | bitwise |
| t4 (GPU binary v1 = be70b99b) vs t4 + nosync (v4 = aa656a75), production A/B 0.1 rotation | final rst data bitwise |

These were repeated on every CPU binary (fast1, fast4, fast5, fast6).

## 3. Well-posed suite (CPU, dt = 20 s, frozen hydro)

`python3 fast/wp_ana.py /viper/ptmp2/jinma/ckfast_0923/wp1 /viper/ptmp2/jinma/ckfast_0923/wp9`,
plus wp4/wp5/wp6/wp7 for the older variants.

Tests:
- **A1**: 200 calls from the A0b steady state. T** is the t4x arm.
- **A2**: transient from the IC, 100 calls, against the t4x dt = 4 s reference.
- **D**: eos_h2, mu0 = 1, against t4x of the same binary.
- **C**: final-state residual and budget.

**Lever 1 is not testable here.** The harness already applies the RT once per call.

"c2/c4/c8" is the final combination of levers 2-5: xstep k + jreuse 0.2 (act 0.25, no
xc) + pred fac 0.5 + nosync (+ pred_chk).

| arm | A1 max vs T** | A2 max vs ref (vs t4) | A2 passes | D max@10 / @100 vs t4x | D passes | C: A2 resfinal, bud/dtS |
|---|---|---|---|---|---|---|
| t4 | 3.4e-8 | 1.18e-4 (0) | 4.02 | 1.4e-4 / 1.6e-6 | 3.37 | 8.6e-9, 9.6e-6 |
| t4x | 0 | - | - | 0 / 0 | 4.89 | - |
| x2 | 3.3e-8 | 1.24e-4 (3.6e-5) | 4.03 | 5.6e-3 / 2.2e-5 | 3.40 | 9.7e-9, 1.0e-5 |
| x4 | 3.3e-8 | 1.79e-4 (1.1e-4) | 4.04 | 1.8e-2 / 6.9e-5 | 3.45 | 8.6e-9, 8.6e-6 |
| x8 | 3.3e-8 | 3.80e-4 (3.0e-4) | 4.05 | 3.7e-2 / 1.6e-4 | 3.44 | 9.1e-9, 1.1e-5 |
| x4 thr 0.01 | 3.3e-8 | 1.71e-4 (9.3e-5) | 4.02 | 1.4e-4 / 5.3e-5 | 3.37 | 8.6e-9, 9.6e-6 |
| j (0.2) | 3.4e-8 | 1.18e-4 (6.0e-8) | 4.02 | 7.4e-5 / 1.6e-6 | 3.39 | 8.3e-9, 5.5e-6 |
| j5 (0.5) | 3.4e-8 | (5.2e-8) | 4.05 | 8.5e-5 / 1.6e-6, **1 non-conv** | 3.40 | 7.4e-9 |
| p5 (pred 0.5) | 3.4e-8 | 1.18e-4 (1.8e-8) | 3.01 | 3.3e-4 / 1.6e-6 | 3.02 | 1.2e-8, 1.2e-5 |
| p1 (pred 1.0) | 3.4e-8 | (1.8e-8) | 3.01 | 3.7e-4 / 1.6e-6 | 3.00 | 1.2e-8 |
| pred + cvkeep, fac 1 | - | - | - | **1.9e-3**, residual 1.1e-7 | 2.05 | - |
| w (warm) | **2.3e-5, drift 6.7e-6/call** (t4: 3.6e-7) | 1.30e-4 (7.0e-5) | 3.07 | **5.6e-3** / 4.5e-5 | 2.90 | - |
| ws (warm_step) | **1.2e-5** | as w | 3.07 | **5.6e-3** | 2.90 | - |
| k (cvkeep) | 3.4e-8 | 1.18e-4 (5.6e-6) | 3.45 | 1.2e-4 / 3.8e-6 | 3.48 | 1.0e-8, 1.1e-5 |
| n (nosync) | bitwise t4 | bitwise | | | | |
| c2 | 3.3e-8 | 1.24e-4 (3.6e-5) | 3.01 | 5.6e-3 / 2.2e-5 | 3.11 | 3.3e-8, 1.8e-5 |
| c4 | 3.3e-8 | 1.79e-4 (1.1e-4) | 3.01 | 1.8e-2 / 6.9e-5 | 3.14 | 3.3e-8, 1.8e-5 |
| c8 | 3.2e-8 | 3.80e-4 (3.0e-4) | 3.01 | 3.7e-2 / 1.6e-4 | 3.12 | 3.3e-8, 1.8e-5 |

What the table shows:
- xstep error grows with k in the transients (A2, D). A1 (steady) is unaffected.
- A threshold of 0.01 removes the D transient error. In production, however, it fires on
  every call (section 4).
- jreuse and pred stay within the solver tolerance of T4, with pred slowly degrading as
  fac grows.
- **warm fails**: the A1 drift is 20x T4's and the D error is 40x. The residual norm
  |R|/(e + 1e-3 e_max) is loose in the thin top, and a seeded state passes it.
- cvkeep combined with pred mispredicts, so fac 0.5 is used and cvkeep is left out.

## 4. Production A/B (GPU)

Setup:
- Restart: bench/cs_hyd4_prod/rst/dhj.00138.rst, read in place: t = 2.1045e7 s, **rotation 69.0**.
- Input: fast/gpu/prod_fast.athinput.
- Each arm runs to t0 + 3.05e4 s (0.1 rotation, about 1665 cycles).
- `python3 fast/ab_ana.py /viper/ptmp2/jinma/ckfast_0923/ab t4` (v5 combos: `.../ab_v5`).
- Jobs: 11954121/122 (v1: t4 t4x o x2 x4 x8 j), 11954181 (v4: n p ws), 11954701 (v5: c2 c4 c8 p5).
- Columns: relative T difference to t4, as rms by region.

**The night side is chaotic.** The t4x arm, which differs only by solver tolerance, is
already at rms 7.8e-3 there. The day side and the kink counts are the meaningful
comparison.

| arm | max day | rms day | rms night | day rms by band (>1e-3, 1e-3..1e-5, 1e-5..1e-7, <1e-7 bar) | mean / max \|ckdesum\| | passes/call | non-conv | kinked cols >0.1 (3 bands) |
|---|---|---|---|---|---|---|---|---|
| t4 | 0 | 0 | 0 | 0 | 3.4e-6 / 5.8e-6 | 4.00 | 0 | 523 374 737 |
| t4x (tol 1e-10) | 6.7e-4 | 3.4e-6 | 7.8e-3 | 1.9e-6 2.9e-6 4.6e-6 7.8e-6 | 2.9e-8 / 2.4e-7 | 5.85 | 0 | 518 372 741 |
| o (lever 1) | 2.0e-2 | **3.0e-4** | 1.3e-2 | 7.0e-5 3.2e-4 2.2e-4 1.1e-3 | 1.2e-6 / 1.5e-6 | 4.00 | 0 | 501 369 710 |
| x2 | 1.3e-2 | 6.7e-5 | 8.2e-3 | 4.9e-5 9.7e-5 6.2e-5 8.0e-6 | 3.4e-6 / 5.9e-6 | 4.00 | 0 | 517 380 736 |
| x4 | 1.1e-2 | 6.7e-5 | 1.0e-2 | 4.4e-5 9.9e-5 6.7e-5 2.2e-5 | 3.5e-6 / 5.9e-6 | 4.23 | 0 | 515 380 **771** |
| x8 | 3.9e-3 | 6.4e-5 | 1.2e-2 | 3.4e-5 9.9e-5 6.5e-5 8.3e-6 | 3.5e-6 / 6.0e-6 | 4.38 | 0 | 518 374 **821** |
| j (0.2) | 1.1e-3 | 2.3e-6 | 4.6e-3 | 2.9e-6 1.1e-6 1.9e-6 2.8e-6 | 3.4e-6 / 5.8e-6 | 4.00 | 0 | 522 375 742 |
| p (pred fac 0.1) | 5.2e-4 | 4.6e-6 | 3.9e-3 | 1.8e-6 4.2e-6 4.6e-6 1.6e-5 | 3.4e-6 / 5.8e-6 | 3.50 | 0 | 519 376 738 |
| p5 (pred fac 0.5) | 3.3e-4 | 5.2e-6 | 3.9e-3 | 1.3e-6 3.2e-6 5.8e-6 2.0e-5 | 3.4e-6 / 9.5e-6 | 3.05 | 0 | 522 374 745 |
| n (nosync) | bitwise t4 | 0 | 0 | 0 | as t4 | 4.00 | 0 | as t4 |
| ws (warm+step) | 5.2e-3 | **1.0e-4** | 1.3e-2 | 6.0e-5 1.7e-4 6.9e-5 1.6e-4 | 3.5e-6 / 6.4e-6 | 4.00 | 0 | 494 377 732 |
| c2 (v5) | 2.0e-2 | 3.0e-4 | 1.4e-2 | 7.0e-5 3.2e-4 2.2e-4 1.1e-3 | 1.2e-6 / 4.9e-6 | 3.07 | 0 | 499 378 716 |
| c4 (v5) | 2.0e-2 | 2.9e-4 | 1.7e-2 | 7.0e-5 3.2e-4 2.2e-4 1.1e-3 | 1.0e-6 / 5.1e-6 | **4.71** | 0 | 498 376 738 |
| c8 (v5) | 2.0e-2 | 3.0e-4 | 2.1e-2 | 7.1e-5 3.2e-4 2.2e-4 1.2e-3 | 3.3e-6 / 4.7e-5 | **5.68** | **103** | 490 372 **791** |

The combination arms (c2/c4/c8) all include lever 1 (once), and their day-side error is
that of lever 1.

Failures kept as evidence:
- **c-arms with jreuse_xc.** With v4 (ab_old/), c4 needed 6.92 passes and had 1
  non-converged call; c8 needed 6.29 passes and had 209 non-converged calls.
- **jreuse together with xstep and once.** diag job 11955056
  (/viper/ptmp2/jinma/ckfast_0923/diag/v5b, 449 cycles), passes per call:

  | arm | passes per call |
  |---|---|
  | once + x4 | 4.37 |
  | + pred | 3.29 |
  | + jreuse 0.2 | about 7 |

  A few columns (20 to 280) keep a slow chord tail (contraction about 0.5). The global
  max-residual trigger does not see them. jreuse_act (2b3e6917) was added for this.
- **xstep_thr.** The max per-cycle relative T/rho change over the ck cells of a rank is
  about 0.08 (median over calls; diag v5 ck8t), so a threshold of 0.01 refreshes on every
  call. The c4t arm cost 60.8 ms/cycle.

## 5. Cost (GPU, apudev, 2 ranks, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1)

Method:
- Same restart as the A/B.
- Each arm runs 170 cycles (t0 + 3000 s).
- Arms are interleaved forward / reversed / forward.
- ms/cycle is taken from the elapsed= lines, and the table gives the median of 3.
- RT = arm - h (the h arm also drops usrsource, 0.9 ms).
- Command: `python3 /viper/ptmp2/jinma/ckprof_0923/ana_time.py /viper/ptmp2/jinma/ckfast_0923/time/<tag>`.

| arm | job (binary) | ms/cycle | RT ms/cycle | RT x hydro |
|---|---|---|---|---|
| h (hydro only) | 11954179 (v4) | 15.33 | - | - |
| t4 | 11954179 (v4) | 126.30 | 111.0 | 7.24 |
| o (lever 1) | 11954179 | 74.01 | 58.7 | 3.83 |
| x4 (lever 2, per stage) | 11954179 | 81.85 | 66.5 | 4.34 |
| x8 | 11954179 | 80.69 | 65.4 | 4.26 |
| j 0.2 (lever 3) | 11954179 / 11955219 (v6, act) | 110.03 / 110.98 | 94.7 / 95.6 | 6.18 / 6.24 |
| j5 0.5 | 11954179 | 109.14 | 93.8 | 6.12 |
| p fac 0.1 (lever 4) | 11954179 | 115.66 | 100.3 | 6.55 |
| p5 fac 0.5 | 11954700 (v5) | 114.60 | 99.2 | 6.46 |
| ws (warm+step) | 11954179 | 126.93 | 111.6 | 7.28 |
| n (lever 5) | 11954179 | 120.86 | 105.5 | 6.88 |
| once+pred+nosync (dn) | 11955116 (v5) | 65.86 | 50.5 | 3.29 |
| d4 / d8 (once + x + pred + nosync) | 11955116 | 46.41 / 45.35 | 31.0 / 30.0 | 2.02 / 1.95 |
| **c2** (all five, k 2) | 11955219 (v6) | 46.54 | 31.2 | 2.03 |
| **c4** | 11955219 (v6) | 42.27 | 26.9 | 1.76 |
| **c8** | 11955219 (v6) | 40.59 | 25.3 | **1.65** |

The earlier v1 profile (prof/v1/c4, v4 prof/v4/c8) shows where the remaining cost is.
The host gap is gone: c4 wall is 40.4 against 35.3 ms of GPU busy, while hydro alone has
4.3 ms. The rest is kernel time:
- A frozen pass costs 2.8 ms: lin1 1.31, lin_sum 0.53, apply 0.46, fused 0.24,
  pre_tp 0.17, pre_opac 0.07.
- A Jacobian build (jlin + sum) costs 1.8 ms.
- A store sweep costs 23 ms plus 3 ms of lin_build.

lin1 is memory-bound: about 5 GB per pass, near HBM peak.

## 6. Choice of k (lever 2)

Evidence by k:
- **Well-posed A2** (error vs t4): k = 2: 3.6e-5; 4: 1.1e-4; 8: 3.0e-4. T4's own error
  against the dt = 4 s reference is 1.2e-4.
- **Well-posed D, call 10**: 5.6e-3 / 1.8e-2 / 3.7e-2.
- **Production day rms**: 6.7e-5 / 6.7e-5 / 6.4e-5, flat in k. Kinked columns in the
  < 1e-7 bar band: 736 / 771 / 821 (t4: 737).
- **Combination** (v5): k = 8 leaves 103 calls non-converged, and each refresh after 8
  stale cycles jolts the solve (dstep 0.3-0.4, 8 passes; diag v5 c8).

**Chosen: k = 2** for accuracy. It is the only k with no kink increase and a transient
error below T4's own. k = 4 is the cost option (c4 42.3 against c2 46.5 ms/cycle). It
stays pending on the v6 A/B, which tests whether jreuse_act removes its 4.71 passes per
call.

## 7. HANDOVER (2026-09-24)

**Commits on ck-fast.** be70b99b, 075b5296, 8a719220, aa656a75, 22d12d91, 2b3e6917 (code),
then the README/scripts commit. The final binary is 2b3e6917:
- CPU: /viper/ptmp2/jinma/ckfast_0923/athena.cpu.fast6
- GPU: athena.gpu.v6, md5 f46dc3b0e341634edae92fb68269cd60

**Per-lever status.**

| lever | built | off-gate | accuracy vs T4-every-stage | cost (RT x hydro, T4 = 7.24) |
|---|---|---|---|---|
| 1 once | existed, works with T4 | n/a (existing switch) | production day rms 3.0e-4 (T4 tolerance floor 3.4e-6), kinks not worse, gap smaller | 3.83 |
| 2 xstep | yes | bitwise | transient error grows with k (A2, D); production day 6.5e-5, kinks +5% (k4) / +11% (k8); thr useless in production | 4.34 (k4, per stage) |
| 3 jreuse | yes (+act, +xc) | bitwise | within tolerance (A2 6e-8, production day 2.3e-6); xc fails (non-conv) | 6.2 |
| 4 pred | yes | bitwise | within tolerance at fac 0.5 (production day 5.2e-6, D@10 3.3e-4 vs 1.4e-4) | 6.46 |
| 4 warm / warm_step / cvkeep | yes | bitwise | warm, warm_step FAIL (A1 drift, D, production day 1e-4); cvkeep OK alone, mispredicts with pred | no gain |
| 5 nosync | yes | bitwise (CPU and GPU) | bitwise | 6.88 |
| combination c2 / c4 / c8 | | | c2: 0 non-conv, 3.07 passes; c4 / c8 pending v6 A/B | 2.03 / 1.76 / 1.65 |

**Merged into rt-integration (2026-09-24).** All levers default off. Gates on the merge
(ck-fast 21d13e2b = b01b50d3 + rt-integration f4e25cac), base = rt-integration f4e25cac:
- CPU (+ wellposed hooks): A2 t4/semi rec.txt and rst, full hydro T4 6 cycles rst + .hst: bitwise
  (/viper/ptmp2/jinma/ckfast_0924/cpugate.out).
- GPU (apudev job 11955755, T4, 0.1 rotation from cs_hyd4_prod dhj.00159.rst, rotation 79.5):
  final dhj.00160.rst data and .hst bitwise.
- A/B 11955220 verdict: c2 passes; c4 marginal (4.56 passes/call); c8 fails (103 non-converged).

**Pending job (done, verdict above).** 11955220 (apudev, v6, AB_EXTRA pred_chk) runs the A/B arms c2 c4 c8 into
/viper/ptmp2/jinma/ckfast_0923/ab/{c2,c4,c8}. Analyse it with:

    cd /viper/ptmp2/jinma/wt_ckfast/tests_ck_implicit/fast
    python3 ab_ana.py /viper/ptmp2/jinma/ckfast_0923/ab t4
    for a in c2 c4 c8; do echo $a $(grep -o 'passes=[0-9]*' /viper/ptmp2/jinma/ckfast_0923/ab/$a/run.log | sort | uniq -c | tr '\n' ' ') nc=$(grep -c NOT-CONV /viper/ptmp2/jinma/ckfast_0923/ab/$a/run.log); done

Pass criteria: 0 non-converged calls; passes about 3-3.5; day rms about 3e-4, which is
lever 1's error; kinks in the < 1e-7 bar band not above t4's 737 by more than 5%.

**Next step, to go from 1.65x toward 1x.** The remaining cost is kernel time:
- about 3.1 passes x 2.8 ms
- the Jacobian at 1.8 ms per build
- the store at 26 ms / k

In order:
1. Fuse lin_sum and the ck part of rt_apply into one kernel: about 0.8 ms/pass. The apply
   re-sums 22 block partials that lin_sum has just written.
2. Store lP in FP32 for lin1, which is memory-bound: about 0.6 ms/pass. This needs an
   accuracy gate.
3. A per-column refresh criterion for xstep (a column-local T/rho change) in place of
   the rank-wide max, so that k can grow without the refresh jolt.
4. Build the chord Jacobian only for the columns that are still active and slow
   (jreuse_act does this partly).
