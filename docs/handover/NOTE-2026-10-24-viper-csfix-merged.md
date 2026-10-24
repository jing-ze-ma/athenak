# NOTE 2026-10-24 (viper -> caltech): csfix-1023 merged into fin51-r-fix (958cac780)

**fork/fin51-r-fix is now 958cac780** (fast-forward from c1681dbe3; the single commit of csfix-1023).
It carries the vet_gd table restore / cs compact-list drop after a failed dt_retry attempt, and the new default
`<time>/dt_retry_t2_be = 0` (was 2).

## Gates (viper; all in /viper/ptmp2/jinma/csfix_1024)
| gate | test | result |
|---|---|---|
| 1 build | CPU none (gcc/14), GPU none_gpu72 (ROCm 7.2), GPU accretor r71 (ROCm 7.1, ry_per_accretor, -ffp-contract=off) | all built, 0 errors |
| 1 style | cpplint (upstream master) on the 4 changed files, c1681dbe3 vs 958cac780 | same 71 messages, all already in c1681dbe3: **no new ones**; no tabs, no new `}}`, `#pragma` lines or trailing whitespace; mode 644 |
| 2 default, no retry deck (CPU) | cbm_r4 with the time2_reject line removed, nlim 6, 1 core | bin data IDENT |
| 2 default, no retry deck (GPU) | D_B1v6 (vphi_acc_1021) with the time2_reject line removed, nx3 128, tlim 4e-4 (29 cycles), apudev 2 GPUs (job 12355681) | 4 hst IDENT, bin data IDENT |
| 3a D_B1v6 prod deck + `dt_retry_t2_be = 2` in the deck | run12 deck acc3d_prod_D_B1v6_extrap2_viper (copy), same keys, 44 cycles (job 12355681) | 4 hst IDENT, bin data IDENT vs c1681dbe3 |
| 3b same deck, new default 0 | (job 12355682) | changes by design, see below |
| 4 restart | 3b deck: full run to t 4e-4 vs restart from its t 1.5e-4 rst (job 12355682) | bin data IDENT; the 18 hst rows after the restart IDENT |
| 5 block-end retry test | | **not done** (skipped on the user's request to make the merge quick). Caltech's A6 RSG test is still the only test of the "vet_gd tables back" path |

Smoke: job 12355532 (all gate binaries and decks, 0 FATAL).

## D_B1v6 production deck: dt_retry_t2_be 0 (new default) vs 2 (old default)
GPU, nx3 128, t 0 -> 4e-4, the same window for both:

| | t2_be = 2 | t2_be = 0 |
|---|---|---|
| cycles to t 4e-4 | 44 | 29 |
| DT_RETRY cuts | 2 (both at cycle 2, dt 1.42e-5 -> 6.37e-6) | 0 |
| backward-Euler redo cycles | 0 | 1 (cycle 2, at the full dt 1.42e-5) |
| hesdirk2 backward-Euler steps | 1 | 2 |
| dt after cycle 2 | 6.37e-6, grows back from cycle 22 (ngood 20) | 1.42e-5 kept |
| dt at the end | 6.94e-6 | 6.94e-6 |
| Picard mean / max | 3.9 / 8 | 6.4 / 150 (the 150 is the BE redo at cycle 2: 11.3 s for that step) |
| NON-CONVERGED | 0 | 0 |
| GPU time (cpu time used) | 41.9 s | 39.3 s |

So on D_B1v6 the new default skips the two cuts and the 20-cycle crawl at half dt. It pays for this with one expensive backward-Euler step, and over this window the cost is about the same. The deck also logs the pre-existing warning "rad_m1/time2_pred_cap never read" (it is in the deck, not related to this change).

## Binaries (958cac780)
- viper accretor r71 (same configuration as the run12 1976bac6 build): `/viper/ptmp2/jinma/csfix_1024/bin/athena_apu_ry3d_958cac78_r71`, md5 a282a9be93358724c2008a83c61f43bc
- viper none_gpu72: `/viper/ptmp2/jinma/builds/bin/athena_none_gpu72_958cac78_csfix`, md5 8f5e9bd10655b18d7d4acfce9c6c31fe
- viper none_cpu: `/viper/ptmp2/jinma/builds/bin/athena_none_cpu_958cac78_csfix`, md5 e3a327e85899f9e2c8846169538cdc5d

## Running productions
- viper run12 (jobs 12327914/16/17/18, still PENDING at merge time) and the Raven HG wedge still read the older 1976bac6 binary.
- A swap of the viper run12 binary to 958cac780 is in progress (USER 10-24). It happens only if 12327914 has not started, and is reported separately.
