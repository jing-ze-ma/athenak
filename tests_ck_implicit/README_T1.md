# T1: `ck_implicit` on the tm sweep (`ck_sweep_form = 1`), production form (2026-09-22)

Phase T1 of `DESIGN_tm.md` (T0 said GO: `T0_RESULTS.md`). Not committed. The patch is
`bench/impl_t1_0922/t1.patch` against HEAD `92d36db3`. It touches
`src/utils/two_stream_rt.hpp` and one comment in `src/pgen/deep_hot_jupiter_rt.cpp`.
Everything was built in `bench/impl_t1_0922/src_new`, a `git archive HEAD` snapshot.

## Verdict

* **Off path: unchanged.** CPU is bitwise in all 5 gates. The production GPU kernel has
  identical instruction count, scratch, VGPR, spill and private-segment figures at every tier.
* **The Jacobian is correct.** The nearest-neighbour entries match a finite-difference
  probe of the tm recurrence to 1e-7 relative. The Newton behaves like the four-pass
  arm's Newton on the same restart.
* **Gate (c) is NOT met, by either form.** From the rot-283 restart, with `maxit = 8`,
  `reuse_jac = 1` and `seed = 2`, none of the 600 calls converges. The residual stalls
  at ~4e-6 and `|ckdesum|` has a median of 5e-6, against the 1e-6 target. The four-pass
  arm (`ck_sweep_form = 0`) stalls in the same way (section 3). The stall comes from the
  Newton driver at this state and dt, not from the tm assembly.
* **Cost (GPU, same binary): implicit-tm is 5.5x the semi-implicit step** (section 4).
  The four-pass implicit arm is 7.6x. The design's estimate for T1 was ~4.4x.
* The tm FATAL is lifted for `ck_sweep_form = 1` only. sd with `ck_implicit`, and
  `ck_impl_frozen_op` with tm, are still refused. **Recommendation: keep the lift in the
  patch but do not commit it** until the stall (section 3) is understood, because the
  convergence gate has not passed.

## 1. The change

* A new **compile-time tag `JAC`** on `rt_chain_ck` (7th tag; `count2.py` prints it in
  the `IMP` column). It is true only in the FRM = 1 instantiation, and only on a pass
  that builds the Jacobian (`ckjacp_`). The production kernel, and passes that reuse the
  Jacobian, run the JAC = 0 code. That code is the tm sweep with no change at all.
* **What JAC computes.** At frozen opacity the tm deposit of cell i is exactly
  `Src_i = (wfc/dz_i)[rat_i D_i - D_{i+1}]`, where
  `D_f = (1+b)(Sc_f - (1-R_f) d_f^+)/(1+R_f b)` and `d_f^+` is the down ray arriving at
  face f. `dR/dB = 0`.
  * Pass 1 carries `dSc_f/dB_{f-2,f-1,f}` through the same Moebius and half-layer maps,
    and parks them at every face. That is 3 Reals per face and chain of private memory,
    in the JAC kernel only.
  * Pass 2 carries `dd^+/dB_{i-1,i,i+1}`. Everything deeper enters through a single
    scalar, `W_f = T_f^2 c_f (alpha_{f+1} W_{f+1} + beta_{f+1})`, because a `B_j` with
    j <= f-2 can reach `d_f^+` only by reflection of `Sc_f`.
  * The nearest-neighbour entries are therefore **exact**, reflection included. Entries
    two or more cells away stay in the residual. Negative per-chain off-diagonal
    contributions (O(beta) reflection) are dropped, which keeps the M-matrix. The
    coefficients come from the `Cc0/Cci/Cco` cache, so there is no extra `expm1`.
* JAC is instantiated for the radial tiers 72, 136 and 264 only. At tier 520 the parked
  window pushes the stack frame to 167888 B, over the 131056 B limit, and the HIP build
  failed. `ck_implicit` + tm with n1 > 264 is therefore a startup FATAL.
  Production is tier 264.
* FATAL at the first RT call: `ck_implicit` with `ck_sweep_form = 2` (sd), and with
  `ck_impl_frozen_op` under tm (T2 is not built). The pgen default is unchanged: an input
  with `ck_implicit` still defaults to form 0, so the tm path must be asked for explicitly.

## 2. Gates: correctness

**(a) Bitwise off path, CPU** (`bench/impl_t1_0922/gate/gate.sh`, `gate.log`, `md5.txt`).
20 cycles; the `.hst` file, the bin dumps and the clean-exit rst (25 files per arm).
REF = `bench/tags_imp_leg_0922/athena.cpu.head`, built from HEAD 3b6e0e52. Its `src/`
differs from 92d36db3 only in `mesh/build_tree.cpp`, which is on the restart path and is
not exercised here. `ck_sweep_form` was added as a key to the gate inputs so that it can
be overridden.

| gate | settings | result |
| --- | --- | --- |
| a_prod_tm | ck_spherical, ck_beam_sph, ck_sweep_form = 1, implicit off (**production**) | BITWISE |
| b_imp_4p | ck_implicit, spherical + beam, ck_sweep_form = 0 | BITWISE |
| c_imp_pp | ck_implicit, plane-parallel | BITWISE |
| d_imp_4p_rj | b + reuse_jac = 1, seed = 2 | BITWISE |
| e_sph_tm_nobeam | spherical, tm, no beam | BITWISE |

**(a') GPU code object** (`count2.py` over `dhj.new.o`; HEAD numbers from
`tags_imp_leg_0922/head.counts2.txt`, whose pgen sources are identical to HEAD's):

| SPH BSP CCH FRM | JAC | tier 72 / 136 / 264 / 520: instr, scratch, vgpr, spill, priv |
| --- | --- | --- |
| T T 2 1 (**production**) HEAD | - | 25682/1528/128/0/17040, 25682/1528/128/0/31376, 25659/1512/128/0/60016, 25545/1457/128/0/117376 |
| T T 2 1 (**production**) T1 | F | **identical on all four tiers** |
| T T 2 1 | T | 29401/2114/128/0/24528, 29401/2114/128/0/45008, 29374/2117/128/0/85984, (not instantiated) |

`new.counts2.txt` lists all 118 chain kernels: 112 plus the 6 JAC ones.

**(b) Finite-difference probe of the assembly** (`bench/impl_t1_0922/fdcheck/fd_tm.py`,
`fd_tm.out`). The script is a line-by-line Python transcription of tm pass 1 and pass 2
plus the JAC algorithm. It runs on random 40-cell columns with a BFace-blended weak cell
and beta up to 0.05. Because the sweep is affine, the central difference is exact.
Result: max |J - FD|/|FD| over the nearest-neighbour entries is 2e-11 to 1.1e-7 across 4
trials. The probe checks the derivation. The C++ transcription is checked only by the
Newton behaviour below: this is a model check, not a probe of the compiled kernel.

**CPU side note** (`gate/cv_tm`, `gate/cv_fp`): `dhj_ck_implicit` cold start, spherical +
beam, maxit 8. tm and four-pass behave alike. The first stage of each cycle hits the cap
in both (res ~6e-7 for tm and ~5e-7 for four-pass at cycle 3), and the second stage
converges in 6-7 passes in both.

## 3. Gates: convergence and gap, GPU (job 11942987, `bench/impl_t1_0922/gpu/`)

Restart `cs_mhd_prod3/rst/dhj.00567.rst`, prod4 `<problem>` form
(`tm_prof_growth/prof_tm` input), plus the keys in `impl.athinput`. 2 ranks on apudev.
Summaries come from `ana.py` over the rank-0 `### ck_implicit` lines.

| arm | form | knobs | calls | passes | res max / median | \|ckdesum\| max / median | capped | non-converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| i1 | tm | maxit 8, rj 1, seed 2 | 600 | 8.00 (cap) | 1.12e-3 / 8.4e-6 | 1.7e-4 / 5.0e-6 | 14 | 600 |
| f1 | four-pass | same | 300 | 8.00 (cap) | 1.15e-3 / 9.4e-6 | 1.5e-4 / 3.9e-6 | 15 | 300 |
| i0 | tm | maxit 8, **rj 0, seed 0** (plain Newton) | 300 | 8.00 (cap) | 1.14e-3 / 1.1e-5 | 1.7e-4 / 3.7e-6 | 14 | 300 |

**Pending: job 11943097** (`gpu/submit2.sh`, `log.out.11943097`). It runs arms `i20`
and `f20` (maxit 20, tm and four-pass) to test whether the stall clears with more
passes. Summarise them with `python3 ana.py i20 f20`.

* From the second cycle on, the last-pass residual of i1 sits at 3.6-5.8e-6, with
  `dstep` ~0.11 (the largest per-cell step of the last pass). This is a **stall, not a
  slow approach**: some cells keep taking 11 % steps. The four-pass arm does the same.
  Plain Newton with the Jacobian rebuilt every pass (i0) stalls in the same way, so the
  cause is not the chord reuse or the seed. The phase-4 figures (5.67 passes, gap
  9.7e-7) were not reproduced by the four-pass arm at this state either.
* **Sweep-to-gas gap.** Semi-implicit (arm s, same binary and input, `rt_desum` every 20
  cycles, 60 lines): mean -9.5 %, range -5.8 to -14.8 %. T0 measured -8.2 % over its
  settled window. The implicit tm arm, even without converging, has `|ckdesum|` with a
  median of 5.0e-6 and a max of 1.7e-4. That is ~4 orders better than semi-implicit, but
  it does not meet the <= 1e-6 target.

## 4. Cost, GPU only (`cpu time used`, 2 ranks, same restart)

Job 11942987 had one repeat per arm for the implicit pair. **Pending: job 11943097**
repeats the s/i1 pair interleaved x2 at 200 cycles (`s_r1 i1_r1 s_r2 i1_r2`, `grep
"cpu time used" log.out.11943097`). Its first arm, s_r1, gave 15.03 s.

| arm | binary | setting | cycles | cpu time used | per 300 cycles | vs semi-implicit |
| --- | --- | --- | --- | --- | --- | --- |
| h1 / h2 | HEAD 3b6e0e52 | production, implicit off | 300 | 24.20 / 21.82 s | | |
| n1 / n2 | T1 | production, implicit off | 300 | 22.29 / 22.17 s | | |
| s | T1 | implicit off (+ rt reports) | 300 | 22.10 s | 22.10 | 1.00 |
| i1 | T1 | implicit tm, rj 1, seed 2 | 300 | 120.91 s | 120.91 | **5.47x** |
| f1 | T1 | implicit four-pass, rj 1, seed 2 | 150 | 83.86 s | 167.7 | 7.59x |

* **Production is unchanged.** HEAD 24.20 / 21.82 s against T1 22.29 / 22.17 s,
  interleaved. The spread between repeats (h1 against h2: 11 %) is larger than the
  difference between binaries, and the code object is identical anyway (a').
* The implicit arms run into the cap on every call (8 sweeps per stage), so 5.5x is the
  cost at the cap. A converging call would cost less.

## 5. Files

`bench/impl_t1_0922/`: `t1.patch`, `build.sh`, `src_new/` (the snapshot and its build
dirs), `athena.cpu.new`, `athena.gpu.new` (md5 d069a7c6...), `dhj.new.o`,
`new.counts2.txt`, `gate/`, `fdcheck/`, and `gpu/` (`submit.sh`, `submit2.sh`,
`log.out.*`, one directory per arm, `ana.py`).
