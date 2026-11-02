# NOTE 2026-11-02 Caltech -> viper: prod-1102 = 4e86ab10f (2nd-order ck keys), please run the HIP bitwise gate

prod-1102 moved by fast-forward from 1f6a938e9 to **4e86ab10f** (branch ck2ndp-1102). It is prod 1f6a938e9 plus the 2nd-order ck cadence keys. Every new key is read only when named; with the key absent the code is bitwise the old code. All new code lives in src/utils/ck_2nd.hpp, with hooks in:
- two_stream_rt_driver.cpp, two_stream_rt_pass.cpp, deep_hot_jupiter_rt_src.cpp;
- the ck restart block (slab ids 26..30);
- driver.cpp: <time>/dt_fixed now also replaces a restart's dt (test key only).

## Keys (USER 11-02: r1m launches with these, cadence 4)
<problem>:
- ck_cad_probe = true: the cadence model is refreshed by a BE probe on U^n every ck_impl_every cycles; both hesdirk2 stages are the linearised step.
- ck_cad_drift = true.
- ck_cad_cont = true, ck_cad_cont_len = 0.25: a continuous model across refreshes.
- ck_impl_xstep = 0.
- ck_vef_dt = V: the VEF is refreshed at the probe that passes t0 + n V, a stateless schedule. r1m uses 672; hydro p1 uses 1500 (~150 dt, today's VEF cost).
- Test-only keys: ck_cad_dt, ck_cad_t0, ck_xstep_dt, ck_cad_drift2, ck_cad_lin_iter, ck_cad_probe_fac / _dt / _tolinc / _twice / _nofsal.

Launch deck: Caltech prod_1102/r1m_launch_v2.athinput (= r1m_launch_gate_es + the keys above).

## Results (H200)
- Frozen-hydro order ladder, double bins, T ratios (1e-6..1e-3 / 1e-3..1 bar), 17 = second order:
  - cadence-4 hesdirk2 today: 4.7 / 4.1, first order;
  - new keys, VEF frozen: 18.7 / 14.2;
  - new keys, VEF on schedule V = 128 s: 14.7 / 10.7.
  - X1 error: 0.48 K / 0.032 K, against 5.0 K / 0.33 K today.
  - The VEF schedule leaves a dt-independent lag bias (0.2-0.3 K / 0.03-0.05 K at V = 128 s). Today's VEF refresh every 150 cycles is first order.
- Cost:
  - hydro W121 p1: 26.9-27.2 ms/cycle, against 27.1 for the split and 26.5 for hesdirk2 + probe + seed;
  - r1m with emf_smooth: 56.4 against 51.6-52.3 (+8 %, accepted by USER);
  - counters: ncap 0, nonconv 0, nfall 0.

## Caltech gates (CUDA, --fmad=false)
- Keys absent vs prod1102m_1f6a938e:
  - hydro W121, 150 cycles: rst + hydro.hst identical;
  - MHD r1m es deck, 40 cycles: rst + mhd.hst identical.
- New keys on the r1m launch deck (es): 300 straight vs 151+149 (split on an hst sync cycle): rst data identical, all 29 post-restart mhd.hst rows identical.
- r1m v2 deck, 100 cycles, 1 vs 2 MPI ranks (two H200, 02-03 GPUs 1+3): rst identical. mhd.hst differs only in the reduction round-off (mass 1.8e-16).
- Accretor cs96 stream restart 830-850 (2 GPU), keys absent, ry_per_accretor 4e86ab10 vs afa17af2: rstcmp SAME (ndiff 0), rc 0, no Kokkos finalize error.

## Request
Please run your HIP bitwise gate on prod-1102 4e86ab10f:
- keys absent, and the r1m keys above;
- in particular the 150+150 continuation with the new keys.
