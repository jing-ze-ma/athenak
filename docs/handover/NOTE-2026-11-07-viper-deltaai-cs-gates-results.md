# NOTE 2026-11-07 viper (on DeltaAI): results of TASK-2026-11-07-deltaai-cs-defaults-gates (hgb + rsb)

Run by viper through the USER's ssh master (jobs csdai_*, account bivj-dtai-gh, ghx4-interactive). Facts only; no default decided here.

## Summary
- **G1 restart (BASE+r8-cs+T): BITWISE** on hgb and rsb (all bins + rst payload, 4 straight vs 2+rst+2).
- **G3 determinism (BASE): BYTE-IDENTICAL** on hgb and rsb.
- **G2 layout (BASE+r8-cs+T) 1/2 vs 4 GPU: NOT bitwise, round-off level**: hgb max rel dens 2e-13, eint 6e-11, m1_e 2e-10, m1_vet D 2e-10; rsb dens/ener 1.3e-13 (velocity/flux components show rel up to 2 only where |value| ~1e-9).
- **A/B 400 cycles, all arms clean**: 0 NOT POSITIVE, 0 NOT ADMISSIBLE, 0 BE redo, 0 failed stages, m1_pos_flr 0 (hgb; rsb has no log block). dt pinned by the cap in both states (571 s / 4000 s 4-cycle states: no stress events, the RSG fb_local floor issue cannot be tested here).
  s/cycle (1 GH200/arm): hgb BASE 0.597, r8 0.597, r8+T 0.631 (+6 %), r8+pert 0.598; rsb 0.477 / 0.473 / 0.494 (+4 %) / 0.471. Picard max on hgb 22 (BASE) -> 8 (r8-cs).
  Final-state difference vs noise (pert 1e-10): hgb r8 vs BASE dens 6.4e-3, r8+T vs r8 5.2e-3, noise 1.5e-10; rsb 1.5e-4 / 5.0e-5 / 3.0e-11.
- GPU memory per GH200: hgb 14.4 GB (4 GPU), 26.8 (2), 62.6-63.9 (1); rsb 11.7 / 21.1 / 70.4-71.4.
- SU: 3550909 9:20 (1.24 SU) + 3551207 18:41 (2.49 SU) = 3.7 SU of 16.

## Stale-build check (Raven segfault root cause = stale incremental build)
DeltaAI's athena_hesnf_gpu_698b8241329e (md5 136cf6db) was built incrementally: only 15 of 25 src/rad_m1/*.cpp recompiled while rad_m1.hpp changed +473/-19. Clean rebuild in a fresh clone: /work/nvme/bivj/jma20/csdai_build_698b/bin/athena_hesnf_gpu_698b8241329e md5 8fc5e319f7daf7cb5b2b0fbd81e23cf8 (256 objects). **Clean vs old: hgb gate bins byte-identical** (G3a BASE and G1s r8+T); rsb was run on the clean binary only. Use the clean binary from now on anyway.
Also built: vstarfix-1106 3b1d3dc19 hesnf recipe, /work/nvme/bivj/jma20/csdai_build_3b1d/bin/athena_hesnf_gpu_3b1d3dc19 md5 81c33106132acf0d7922febd0e71c6d6 (for BSG).

## BSG holder (USER: queue bsg on deltaai, keep the allocation)
csdai_bsgH 3550945: ghx4 1 node 48 h (<= 192 SU, exits after 60 min idle), PENDING (est. start 11-08 16:38 CST at submission). Queue /work/nvme/bivj/jma20/bsgnv_1107/queue: 05_smoke_v3b (10 cycles; failure -> STOP in long/armB only), 10_long_v3b (claim bsg_long; Raven owns it -> skips), 20_armB_v3b (2 d, vet_gd_tangential=false, claim bsg_step0); all on 3b1d3dc19, all preflighted. Run dirs on /work/hdd/bivj/jma20/bsgnv_1107 (symlinked).

Details: viper /viper/ptmp2/jinma/defaults_1107/cs_dtai/RESULTS.md (copied below).

# cs default-candidate gates on DeltaAI (TASK-2026-11-07-deltaai-cs-defaults-gates), run by viper via the USER's ssh master

Jobs (account bivj-dtai-gh, ghx4-interactive, 1 node = 4 GH200, billing x2 = 8 SU/node-h):
- csdai_hgb 3550909: hgb, OLD binary /u/jma20/ATHENAK/builds/bin/athena_hesnf_gpu_698b8241329e (md5 136cf6db), 9:20, COMPLETED, 1.24 SU
- csdai_both 3551207: rsb then hgb, CLEAN rebuild of 698b8241329e /work/nvme/bivj/jma20/csdai_build_698b/bin/athena_hesnf_gpu_698b8241329e (md5 8fc5e319f7daf7cb5b2b0fbd81e23cf8), 18:41, COMPLETED, 2.49 SU
- total 3.7 SU of the 16 SU budget. (A first csdai_rsb 3550977 on the old binary was cancelled while PENDING.)

States (no rotation: stellar runs; time / cycle of the restart): hgb_hegiant.00001.rst t=571.43 s cycle 4 (md5 596fe2d4), rsb_rsg.00001.rst t=4000 s cycle 4 (md5 4c1c93d0); tables repointed to /work/nvme/bivj/jma20/cs_1107/J (raven_deps_1107.tar, md5 426a2d8e). Both: 24 MeshBlocks (hgb 704x16x16, rsb 548x16x16).
Steps: R1 = +1 cycle BASE from the Caltech rst (writes the extended rad_m1 record); all gates restart from R1. Env: module reset, MPICH_GPU_SUPPORT_ENABLED=1, UCX_TLS=^xpmem, 1 rank per GH200 via CUDA_VISIBLE_DEVICES (gpuwrap.sh), plain srun --cpu-bind=cores.
Preflight (preflight.py = raven pool_1107 copy, identical to the branch copy): all 12 command lines per state OK on the login node and again in-job.
Run trees: /work/nvme/bivj/jma20/cs_1107/runs/{hgb.3550909,hgb_clean.3551207,rsb_clean.3551207}; logs /work/nvme/bivj/jma20/cs_1107/logs.

## Stale-build check
The old binary came from the incremental hesnf_gpu tree: make log make_hesnf_gpu_698b8241329e.log recompiled 67 objects, only 15 of the 25 src/rad_m1/*.cpp, while rad_m1.hpp changed +473/-19 lines since the previous build (ff1f52bc). Clean rebuild in a fresh clone (256 objects): md5 8fc5e319. **Clean vs old, hgb G3a (BASE) and G1s (BASE+r8+T) bins: FULL FILES BYTE-IDENTICAL** (hydro_w, m1, m1_vet); all gate verdicts identical. The old-binary hgb results stand for these keys.

## Gates (all runs: rc 0, no FATAL, no NaN, summary NON-CONVERGED=0)
| gate | hgb (old and clean binary, identical) | rsb (clean) |
|---|---|---|
| G1 restart, BASE+r8+T, 4 straight vs 2+rst+2, 4 GPU | BITWISE: hydro_w, m1, m1_vet bins + rst payload | BITWISE: hydro_u bin + rst payload |
| G2 layout, BASE+r8+T, 1 vs 4 GPU | DIFF (round-off): max rel dens 7.4e-14, eint 9.5e-12, m1_e 2.9e-11, m1_f1 3.8e-11, m1_vet D 5.3e-11; vel2/3, m1_f2/3 rel up to 2.0 only where |value| ~1e-9 (sign of ~0) | DIFF (round-off): dens 1.3e-13, ener 1.3e-13, mom1 2.0e-6 / mom2,3 up to 2.0 on ~0 values |
| G2 layout, 2 vs 4 GPU | DIFF: dens 2.1e-13, eint 5.6e-11, m1_e 1.9e-10, m1_f1 1.1e-10, m1_vet D 2.4e-10 | DIFF: dens 1.2e-13, ener 1.2e-13 |
| G3 determinism, BASE x2, 4 GPU | BYTE-IDENTICAL (bins + rst) | BYTE-IDENTICAL (bin + rst) |

Layout dependence is at round-off (first differences appear in the first 4 cycles at ~1e-16 and grow to <=2e-10 relative in E/eint). Not bitwise, so a layout change mid-production is not reproducible bit-for-bit.

## A/B (400 cycles from R1, 4 arms concurrently, 1 GH200 each, same layout)
| state | arm | s/cycle | Picard mean / max | NOT POS | NOT ADM | BE redo | m1_pos_flr | dt | GPU mem |
|---|---|---|---|---|---|---|---|---|---|
| hgb (old bin) | BASE | 0.597 | 3.287 / 22 | 0 | 0 | 0 | 0 | 142.857 (pinned) | 63.9 GB |
| hgb | +r8 | 0.597 | 3.290 / 8 | 0 | 0 | 0 | 0 | 142.857 | 63.9 GB |
| hgb | +r8+T | 0.631 | 3.286 / 8 | 0 | 0 | 0 | 0 | 142.857 | 63.9 GB |
| hgb | +r8+pert 1e-10 | 0.598 | 3.290 / 8 | 0 | 0 | 0 | 0 | 142.857 | 63.7 GB |
| hgb (clean bin) | BASE / r8 / r8T / pert | 0.614 / 0.615 / 0.648 / 0.615 | same as above (identical) | 0 | 0 | 0 | 0 | 142.857 | 63.9 GB |
| rsb (clean) | BASE | 0.477 | 3.022 / 4 | 0 | 0 | 0 | n/a (no log block) | 1000 (pinned) | 71.4 GB |
| rsb | +r8 | 0.473 | 3.026 / 4 | 0 | 0 | 0 | n/a | 1000 | 71.4 GB |
| rsb | +r8+T | 0.494 | 3.026 / 4 | 0 | 0 | 0 | n/a | 1000 | 71.4 GB |
| rsb | +r8+pert 1e-10 | 0.471 | 3.026 / 4 | 0 | 0 | 0 | n/a | 1000 | 71.4 GB |

T (vet_gd_tangential) costs +6 % (hgb) / +4 % (rsb) per cycle; r8-cs is cost-neutral. Max Picard drops 22 -> 8 with r8-cs on hgb. Concurrent arms: timings comparable between arms, not to a dedicated-node run.
Final-state differences (max rel, last common dump): hgb r8 vs BASE dens 6.4e-3, eint 6.4e-3, m1_e 4.1e-4; r8+T vs r8 dens 5.2e-3, m1_e 4.7e-4; noise (pert vs r8) dens 1.5e-10, m1_e 1.4e-9. rsb r8 vs BASE dens 1.5e-4; r8+T vs r8 5.0e-5; noise 3.0e-11. So both key sets change the solution far above the noise floor (they are scheme changes); whether closer to physics is not decided by this test.
Caveats: dt is pinned by a cap in both small states, so the dt criterion does not discriminate; these states are 4-cycle fresh starts (t = 571 s / 4000 s), so no stress events (zero NOT POS / BE redo in all arms) and the Caltech "fb_local x2.6 floors" on RSG cannot be checked here (rsb has no log output; hgb m1_pos_flr = 0 in all arms).

## Per-run tables (ana.py)
Note: the columns "BE redo lines" and "fallback lines" count log lines that mention the words (the end-of-run summaries), not events; the event counts are in the summaries (their BE redo 0, failed stage attempts 0 in every A/B arm).
### hgb.3550909
| run | ranks | cycles | s/cycle (excl. 1st) | last dt | Picard (last summary) | NOT POS | NOT ADM | BE redo lines | FBDIAG | fallback lines | m1_pos_flr (log) | GPU mem max MiB/GPU | rc | FATAL | NaN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R1 | 4 | 2 | 5.319 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00 | 13309 / 13309 / 13349 / 13349 | 0 | 0 | 0 |
| G1s | 4 | 5 | 1.856 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| G1a | 4 | 3 | 2.716 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 14145 / 14143 / 14217 / 14215 | 0 | 0 | 0 |
| G1b | 4 | 3 | 2.982 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| G2n1 | 1 | 5 | 1.513 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 62609 | 0 | 0 | 0 |
| G2n2 | 2 | 5 | 1.754 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 26781 / 26781 | 0 | 0 | 0 |
| G3a | 4 | 5 | 1.664 | 1.428571e+02 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| G3b | 4 | 5 | 1.660 | 1.428571e+02 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| AB_base | 1 | 401 | 0.597 | 1.428571e+02 | mean=3.287440e+00 max=2.200000e+01 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
| AB_r8 | 1 | 401 | 0.597 | 1.428571e+02 | mean=3.289855e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
| AB_r8T | 1 | 401 | 0.631 | 1.428571e+02 | mean=3.286232e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
| AB_r8pert | 1 | 401 | 0.598 | 1.428571e+02 | mean=3.289855e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
### hgb_clean.3551207
| run | ranks | cycles | s/cycle (excl. 1st) | last dt | Picard (last summary) | NOT POS | NOT ADM | BE redo lines | FBDIAG | fallback lines | m1_pos_flr (log) | GPU mem max MiB/GPU | rc | FATAL | NaN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R1 | 4 | 2 | 4.707 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00 | 12503 / 12503 / 12567 / 12567 | 0 | 0 | 0 |
| G1s | 4 | 5 | 1.754 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| G1a | 4 | 3 | 3.011 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| G1b | 4 | 3 | 3.215 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| G2n1 | 1 | 5 | 1.565 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 62609 | 0 | 0 | 0 |
| G2n2 | 2 | 5 | 1.895 | 1.428571e+02 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00 | 26781 / 26781 | 0 | 0 | 0 |
| G3a | 4 | 5 | 1.900 | 1.428571e+02 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00 | 14339 / 14337 / 14411 / 14409 | 0 | 0 | 0 |
| G3b | 4 | 5 | 1.820 | 1.428571e+02 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00 | 14403 / 14401 / 14475 / 14473 | 0 | 0 | 0 |
| AB_base | 1 | 401 | 0.614 | 1.428571e+02 | mean=3.287440e+00 max=2.200000e+01 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
| AB_r8 | 1 | 401 | 0.615 | 1.428571e+02 | mean=3.289855e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
| AB_r8T | 1 | 401 | 0.648 | 1.428571e+02 | mean=3.286232e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63869 | 0 | 0 | 0 |
| AB_r8pert | 1 | 401 | 0.615 | 1.428571e+02 | mean=3.289855e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | 0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00+0.00000e+00 | 63681 | 0 | 0 | 0 |
### rsb_clean.3551207
| run | ranks | cycles | s/cycle (excl. 1st) | last dt | Picard (last summary) | NOT POS | NOT ADM | BE redo lines | FBDIAG | fallback lines | m1_pos_flr (log) | GPU mem max MiB/GPU | rc | FATAL | NaN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R1 | 4 | 2 | 8.689 | 1.000000e+03 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | - | 11641 / 11641 / 11717 / 11709 | 0 | 0 | 0 |
| G1s | 4 | 5 | 2.688 | 1.000000e+03 | mean=3.375000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | - | 11643 / 11643 / 11719 / 11711 | 0 | 0 | 0 |
| G1a | 4 | 3 | 5.994 | 1.000000e+03 | mean=3.500000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | - | 11643 / 11643 / 11719 / 11711 | 0 | 0 | 0 |
| G1b | 4 | 3 | 5.851 | 1.000000e+03 | mean=3.250000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | - | 11643 / 11643 / 11719 / 11711 | 0 | 0 | 0 |
| G2n1 | 1 | 5 | 1.318 | 1.000000e+03 | mean=3.375000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | - | 70393 | 0 | 0 | 0 |
| G2n2 | 2 | 5 | 1.459 | 1.000000e+03 | mean=3.375000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | - | 21097 / 21075 | 0 | 0 | 0 |
| G3a | 4 | 5 | 1.552 | 1.000000e+03 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | - | 11641 / 11641 / 11717 / 11709 | 0 | 0 | 0 |
| G3b | 4 | 5 | 1.519 | 1.000000e+03 | mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | - | 10901 / 10901 / 10965 / 10965 | 0 | 0 | 0 |
| AB_base | 1 | 401 | 0.477 | 1.000000e+03 | mean=3.022005e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | - | 71377 | 0 | 0 | 0 |
| AB_r8 | 1 | 401 | 0.473 | 1.000000e+03 | mean=3.025672e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | - | 71377 | 0 | 0 | 0 |
| AB_r8T | 1 | 401 | 0.494 | 1.000000e+03 | mean=3.025672e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 7 | - | 71377 | 0 | 0 | 0 |
| AB_r8pert | 1 | 401 | 0.471 | 1.000000e+03 | mean=3.025672e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | 0 | 0 | 1 | 0 | 6 | - | 71427 | 0 | 0 | 0 |
