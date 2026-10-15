# NOTE caltech -> raven: He giant LARGE gate on FIN5.1 (reply to aa27fbd6): A PASS, B PASS, C stopped by the USER at t ~ 9.7e4 s; same E0 deficit as FIN4

Dir /resnick/groups/carnegie_poc/jingze/hegiant_large_1014/gate51 (old gate/ untouched). Binary athena_gpu_he_star_m1_21c4839d_nofma
md5 f06db8e6c49a5014db07577289196d8b (build + smoke: NOTE-2026-10-15-caltech-fin51-build.md, 6fee27be). 1 node x 4 H200, 4 ranks x
1 MeshBlock 1697x64x64, UCX_TLS=^xpmem. USER 10-15: the He giant LARGE run is DROPPED (GPUs freed); C was stopped before its end;
no production was or will be submitted.

## Deck
hegiant_large_fin51.athinput **md5 8ee79d2c4367919706dca0cdf496fbd6** = 09a3b446 + `dt_retry = true` under `<time>`, with the inline
`# ...` comments stripped (408 keys parse identical, `<comment>` block untouched). Reason: with the comments the FIN5.1 rst header is
45 282 B, over the 11 x 4096 = 45 056 B limit of ParameterInput::LoadFromFile -> every restart FATALs ("<par_end> is not found in the
first 40KBytes", seen in the smoke; details in the build NOTE sect. 2b). Stripped: 39 194 B. (Commented version kept as
hegiant_large_fin51_withcomments.athinput, 8a490843.)

## A. smoke (job 4542525, 24 cycles, keys as the FIN4 gate): PASS
- rc 0; FATAL 0, NaN 0, DIVERGED 0, NON-CONVERGED 0; Picard mean 3.81 max 9; dt 20.0 every cycle
- **NO `ic_rad: no fused twin intensities` line** (A, B and C: 0)
- s/cycle steady 0.86-0.91 (FIN4 0.8-0.9); output cycles 12 (+31 s) and 24 (+70 s)
- **GPU peak 56 059 MiB** per H200 (FIN4 gate: 77 775 MiB, -28 %), incl. the dt_retry checkpoint 7.17 GB on rank 0
  (24 checkpoints, 0 retries)
- `tfloor_kelvin = 1200 K`; `sponge rate = 0.001 above r = 4.38291e+13`; `band-clamped lateral reads (ALL ranks, all sweeps) 0`;
  vet_col_lat clamped 0; out-of-slab main 0 twin 0; dense band saves 0
## B. restart from A/rst/hegiant.00001.rst to 24: PASS (BITWISE)
- final hydro_w / m1 / m1_vet dumps SAME; user.hst + hydro.hst 12 of 12 rows verbatim in A; Picard 3.20/4, NC 0

## C. background (he_seed 0, job 4542526): stopped by the user at cycle 4800, t = 9.67e4 s (of ~1.3e5 planned)
- up to the stop: FATAL/NaN/DIVERGED 0, NC 0, DT_RETRY 0, no-fused-twin 0; GPU peak 56 005 MiB; Mdot_top 0 in every row
- floor bands at 9.6e4 s: 1.50-3.23 R_ph -0.08 km/s, 3.23-6.00 -0.02, 6.00-r0 -0.01, >r0 -0.00; max|drho/rho| <= 0.002: PASS
- bands 0-1.5 R_ph (km/s / max|drho/rho|), FIN5.1 C vs E0 vs the FIN4 C (75437e36):

| t (s) | band | FIN5.1 C | E0 | FIN4 C |
|---|---|---|---|---|
| 1.0e3 | 0.99-1.01 | +0.19 / 0.020 | +0.27 / 0.024 | +0.19 / 0.020 |
| 2.0e3 | 0.99-1.01 | +0.36 / 0.045 | +0.73 / 0.078 | +0.36 / 0.045 |
| 5.0e3 | 0.99-1.01 | +0.70 / 0.143 | +1.89 / 0.276 | +0.70 / 0.143 |
| 1.0e4 | 0.00-0.50 | +0.12 / 0.001 | +0.29 / 0.002 | +0.12 / 0.001 |
| 1.0e4 | 0.99-1.01 | **+1.15 / 0.330** | **+2.69 / 0.510** | +1.15 / 0.330 |
| 1.0e4 | 1.03-1.10 | -0.17 / 0.025 | -0.17 / 0.035 | -0.17 / 0.025 |
| 1.0e4 | 1.10-1.50 | -0.15 / 0 | -0.15 / 0 | -0.15 / 0 |
| 2.3e4 | 0.99-1.01 | +1.92 / 0.557 | (E0 ends 1.09e4) | +1.92 / 0.557 |
| 4.3e4 | 0.99-1.01 | +5.39 / 2.120 | - | +5.38 / 2.131 |
| 9.6e4 | 0.99-1.01 | +12.82 / 67.2 | - | - |

- **Verdict: the photosphere band still FAILS vs E0 (x0.43 at 1e4 s), and FIN5.1 C equals the FIN4 C to the printed digits up to
  4.3e4 s.** So the FIN4 IC overflow (skipped twin J correction) is NOT the cause of the E0 mismatch; the candidates of 75437e36
  (the 1697 grid near R_ph/base, tfloor 1.2e3 / logt_min 3.0, the larger thin top) remain. Not pursued (run dropped).
- analysis: bg_large.py (eada8978) as in 75437e36: `BASE=hegiant LOGF='run*.log' VISPY=<FIN5.1 source>/vis/python python3
  bg_large.py 4.30368e12 4.38291e13 gate51/C 1` (the repo checkout's bin_convert.py does not read bin format 1.2, so VISPY points
  to the FIN5.1 build source). Full table: docs/handover/hegiant-large-1014/bg_C51_caltech.txt.

## Production: NOT submitted (run dropped by the user)
Prepared only, in case it is revived: run/run.cfg -> BIN 21c4839d (md5 f06db8e6), IN hegiant_large_fin51.athinput (8ee79d2c),
TLIM/tlim 1e9, UCX_TLS=^xpmem (FIN4 copy run/run.cfg.fin4_d505d21f); link_hgl.sh unchanged apart from a header line (7-d, hgLw,
--exclude=hpc-sm-01-09); rmax_guard_large.sh unchanged (RLIM 525). run/ holds no rst or logs.

## 200-Rsun FIN4 He giant fresh start (restated from 17c4957b)
hegiant_1014/fin4/run.log (job 4489484, chain hgF4) DOES contain, once:
`<rad_m1> ic_rad: no fused twin intensities: J is the raw formal solution (ray noise not removed)`.
So that FIN4 production began from the degraded radiation IC (no twin J correction). Nothing was changed on that chain.

## Not done
- C2/C3 FIN4 BGs (4522986/4522987) left to the coordinator, as instructed.
