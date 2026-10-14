# NOTE caltech -> raven: He giant LARGE (x1max 700 Rsun) gate: A PASS, B PASS, C (BG) FAIL on the photosphere band -> production NOT submitted

Reply to TASK-2026-10-14-caltech-hegiant-large.md (921c86ff). Dir /resnick/groups/carnegie_poc/jingze/hegiant_large_1014.
Your 200-Rsun chain (hgF4w) and hgGuard untouched.

## Build / files
- binary athena_gpu_he_star_m1_d505d21f_nofma md5 1440c05320a7bef88c4ad960019b9840 (no new build)
- deck hegiant_large.athinput.in 005be1fc -> after @FILES@ substitution hegiant_large.athinput 09a3b446c4c267898f1b081910ab52d7
- IC ic_giant_own_x800.txt ac69934e7d17951d708c4e917f95a46a (head 077ddfd6, tail 242826ad, ddc72692); tables c1ef063c / b699c60d
- gate job 4515555 (1 node x 4 H200, hpc-sm-02-10, UCX_TLS=^xpmem), 4 ranks x 1 MeshBlock 1697x64x64

## A. smoke C (24 cycles): PASS
- rc 0; FATAL 0, NaN 0, DIVERGED 0, NON-CONVERGED 0; Picard mean 4.94 max 11; dt 20.0 every cycle
- s/cycle steady 0.8-0.9 s (cycles 3-24 incl. the 2 output spikes: rst write at cycle 12 +37 s, final rst+dumps +71 s -> 5.98 s mean)
- GPU peak 77775 MiB per H200 (of 141 GB)
- `<hydro>/tfloor_kelvin = 1200 K`; `sponge rate = 0.001 above r = 4.38291e+13` (r_top 4.8699e+13)
- vet_gd halo: `band-clamped lateral reads (ALL ranks, all sweeps) 0`; vet_col_lat clamped reads 0
## B. restart R from A/rst/hegiant.00001.rst to 24: PASS (BITWISE)
- final hydro_w / m1 / m1_vet dumps data SAME; user.hst and hydro.hst rows after cycle 12: 12 of 12 verbatim in A; Picard 3.20/4, NC 0

## C. background BG (he_seed 0, -t 01:55:00): FAIL (photosphere band does not match E0)
- rc 0, reached cycle 6619, t = 1.310e5 s; FATAL/NaN/DIVERGED 0, NC 0, Picard 4.45/11; GPU peak 77201 MiB; ~1.0 s/cycle with dumps every 1000 s
- budget: (L_in + P_esrc) x 12 = 6.7551e6 Lsun: PASS. Mdot_top = 0 in every row: PASS
- floor bands: 1.50-3.23 R_ph |<v_r>| 0.08, 3.23-6.00 0.02, 6.00-r0 0.01, >r0 0.00 km/s; max|drho/rho| <= 0.002: PASS
- bands 0-1.5 R_ph at t = 1.0e4 s (km/s / max|drho/rho|), this run vs E0:
  - 0.00-0.50: +0.12 / 0.001 vs +0.29 / 0.002   (x0.4)
  - 0.99-1.01: +1.15 / 0.330 vs +2.69 / 0.510   (x0.43) <- FAIL
  - 1.01-1.03: -0.04 / 0.003 vs -0.11 / 0.008
  - 1.03-1.10: -0.17 / 0.025 vs -0.17 / 0.035
  - 1.10-1.50: -0.15 / 0     vs -0.15 / 0       (match)
- same deficit vs E0 at every time: t 1e3 +0.19 vs +0.27, 2e3 +0.36 vs +0.73, 5e3 +0.70 vs +1.89; at 2.3e4 +1.92 vs E0 +5.2.
- cross-check with OUR 897-cell FIN4 run (same binary and FIN4 physics; seeded) at equal t, band 0.99-1.01:
  t 2.16e4: 897 +5.16 / 1.11 vs large +1.70 / 0.51; t 4.32e4: 897 +7.84 / 9.98 vs large +5.6 / 2.3.
  0-0.5 band: 897 +0.40 vs large +0.22. So the 897 FIN4 run agrees with E0 and the 1697 grid is ~2-3x slower at the photosphere
  and ~2x lower in the deep band in the first ~0.5 d; by 1.3e5 s the photosphere band reaches +12 km/s, 1.03-1.10 +40 km/s
  (the envelope moves into the floor gas there, drho/rho ~6e5 is expected for floor cells).
- Candidates (not tested here): the new grid near R_ph/base (dr/r 0.8-1.3x of fit2_897), tfloor 1.2e3 / eos_logt_min 3.0,
  or the radiation field of the larger thin top (marshak at 700 vs 200 Rsun). A 1697-grid BG with tfloor 3.0e3/logt_min 3.4, or
  an 897 BG with tfloor 1.2e3, would separate them (each ~1 h on 4 H200). Tell us which; we run it on request.
- full table: docs/handover/hegiant-large-1014/bg_C_caltech.txt

## Production: NOT submitted (gate rule). Ready to go on your word:
- link script hegiant_large_1014/link_hgl.sh (copy of link_hg4w.sh, job name hgLw, 7-d, --exclude=hpc-sm-01-09), run/run.cfg
  (UCX_TLS=^xpmem, tlim 1e9, prune keep 3), plan link 1 + 2 afterany
- rmax guard rmax_guard_large.sh (RLIM 525, expansion partition, 7 d, + 1 afterany follow-on); its numbers come from a minimal
  reducer guard/reduce_rmax.py on the NEW run's dumps (same R_ph/r12 definition as the page; checked on the 897 run: t 18.50 d
  R_ph 86.991, r12 100.179 = the page). Stops only hgLw jobs.
