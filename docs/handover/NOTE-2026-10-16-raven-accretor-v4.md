# NOTE 2026-10-16 raven: accretor v4 (HANDOVER 01fdd1c3) - INTERIM (stopped on USER credit limit)

Code accretor-bigdom-1015 @ 99adf52d (v4 code 5b0d4bc1). Raven dir /raven/ptmp/jinma/accbig_1016 (v4/); the full
RESUME (jobs, commands, open checks) is at the end of /raven/ptmp/jinma/accbig_1016/NOTE.md.

## Builds (nofma, login node; interactive was ~6 h out)
5b0d4bc1: md5 ac371900ebe237ba0f3980a361669478. Test binary 5b0d4bc1 + thinclos-1016 (02ddf595, 34b1b5f8
cherry-picked, local branch accv4-thinclos-test 2e8f2bfe, not pushed): md5 4f0dbbfc74b86c7951d7a7a90a55e4c2.

## 1. Regression gate: PASS (BITWISE), job 31228556, 4 A100, 20 cycles
- v3 deck + dt_report, donor_irr = false: efe33e6a == 5b0d4bc1 == thinclos test binary (key absent).
- Old decks s5 (G1, fresh) and s5 + stream + ROL2U (G2): orbit1 fdb18517 == 5b0d4bc1.
Hst (hydro, user), bins (hydro_w, m1, whole files incl. the parameter header) and the cycle/dt(/dt_report) lines
are identical in every pair.

## 2. Start-up check: PASS
gdark in-plane T_loc 7.78 (alpha 0.1 deg) / 16.77 (5.1) / 21.33 (15.1) / 24.19 (30.1) / 26.15 (60.1) / 26.22 (90.1) /
25.50 (120.1) / 24.19 kK (179.9). Design values: 7.8/16.9/21.4/24.2/26.2/26.2/25.5/24.2. Pole 27.6 kK, A_lobe 2133.4,
T_mean 25933.5 K. s_geo 0.635505 (0.636). Target L_d/(4 pi D^2) = 152225 (1.52e5).

## 3. Calibration: s = 0.22 (used); recommend s ~0.30 next round
- At r 9.6 (the TASK's radius) the metric is ill-posed: the donor light is absorbed and re-emitted in the atmosphere
  there. -(dF_r) is 6.8e3 at s = 1, which would give s = 22. Used instead: the facing (|phi| < 10 deg) incident flux just
  above the atmosphere, at r 10.0 / 10.4.
- s = 1 (gdark 0.25; v4/s1, job 31228557): 6.67e5 / 7.42e5, i.e. 4.4 / 4.9x the target (corrected for D: s = 0.236 /
  0.220). Static dt 7.03e-6 flat even at s = 1; atmosphere r 9.1-9.6 max |v| 115 km/s at t 0.00946 (no irradiation:
  47; v3 uniform: 228).
- s = 0.22 (deck v4/run/cal.athinput; static v4/cal/static, job 31229071): 1.10e5 / 7.96e4 = **0.73 / 0.52x the
  target**. The response is nonlinear (~s^1.5, gas coupling), so **s ~0.30** would hit the target. Not switched mid-test.
- Static at s = 0.22: **dt 7.03e-6 flat to t 0.00946**, 0 retries, 0 NC. Atmosphere max |v| 11.7 / 25.5 / 46.1 / 49.6
  km/s at t 0.0024 / 0.0047 / 0.0071 / 0.0095 (no irradiation: 11.7 / 22.1 / 39.5 / 46.5). Envelope vs orbit1 at t
  0.00118: rho > 0.01 median relative difference 6e-7 (max 1.2e-2). Lobe at floor with v = 0 (sink_zero_v works).
  Ltop 5.1e4 (orbit1 4.5e4, v3 2.0e5). Eirr rate 2.6e5 per time unit (v3: 7.4e6).

## 4. v4 stream phase through impact: RUNNING (not evaluated)
gpudev links 31229257 (done), 31229258 (running), 31229259 (queued), and gpu 31229101 (race), dir v4/cal/stream.
**At t 0.0185 (cycle 3700): dt 3.36e-6, 0 retries.** dt drops after the switch-on even at s = 0.22. The cause is not yet
located (the v3 drop was the atmosphere; here the static phase was clean). Open checks: see RESUME.

## 5. Fix 3 (vet_col lateral beam): INCOMPLETE
- thinclos m1 override (s = 1) on Raven: dt collapse by cycle 30 (eint <= 0 from cycle 5, Picard NC resid 1e16). The
  debug M1 thin-layer closure is not usable as is.
- Arms queued on the gpu partition, 31229100: vet_col without irradiation, Eddington with and without irradiation, M1
  without irradiation, M1 at s = 0.22, all to t 0.00237. Compare with the vet_col s = 1 arm (v4/s1 m1 bin at t 0.00237).
- **BUG (MI300A, not debugged):** the accretor pgen built for Viper MI300A (rocm 7.2, nofma, md5 dc0815b0) dies on
  cycle 1 in EVERY arm with a HIP illegal memory access (HSA_STATUS_ERROR_MEMORY_APERTURE_VIOLATION in a Kokkos
  deep copy). That includes plain vet_col with donor_irr off, dt_report off and dt_retry off (Viper jobs
  12246122/12246124/12246144).
