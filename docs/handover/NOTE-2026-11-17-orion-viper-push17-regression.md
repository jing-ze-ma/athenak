# NOTE 2026-11-17 Orion -> viper (USER): push 17 CPU regression of b66ff2552 vs dd468ab99

Answers TASK-2026-11-17-viper-orion-push17-regression.md through UPDATE 2 (2455c5c78). This replaces the
preliminary version (d3a90ed90) of this file.

**A: bitwise. Restart: bitwise.** On b66ff2552 with `clean_yres = false` written into every deck, all reg runs are
bitwise to dd468ab99, including multi-rank, the x1line deck, WB and the RT feed decks; the feed blocker is fixed by
7c4805666.
**B (no keys, so clean_yres = local) is larger than the expected ~2e-8, mostly from yres = local:** flz_q changes by
up to 0.9 of max|E|, the RT feed decks by 6.7e-3, and the others by 1e-7 to 3e-6. Scores, conservation and WB levels
are unchanged, and there are no warnings.
10/10 restarts are bitwise in A and in B.
**Opt-in arms** (run on 243146188, which has the same opt-in code): restarts 9/9 bitwise. The forced retry with the
float set is not bitwise on sp3 and thin1, which matches your clean_extrap_share note.

## Binaries (Orion, Release, gcc/13 + openmpi/4.1, -ffp-contract=off, Athena_M1_VGD_NMAX=4096; git archive builds)
| head | md5 |
|---|---|
| b66ff2552 (final candidate) | `90608574d0593dbaa47033327e950cc6` |
| 243146188 | `4fba0a71244eea03ea8fd134a5fc6ccb` |
| e5c54306a | `cf52f1cec16c3441acc78ab5e00782c4` |
| dd468ab99 (push 16) | `824b439cbb0bcdf6d241a781357be5d7` |

- All runs ran in `srun --overlap` steps of hold 230761 (orion002).
- The b66ff2552 arms ran through /orion/ptmp/jinma/clean_1115/tools/pick.sh, at the harness thread counts.
- Six runs lost their first attempt because pick.sh found no free CPU and called taskset with an empty list. They
  were rerun. A guard in pick.sh for an empty CPU list would avoid this.

## A: b66ff2552 with clean_yres = false vs dd468ab99
- **All runs: bin data blocks, tab/hst and every CLEAN line identical.** These are:
  - the cc ladder; thin_p/r/f, dif, lat, pulse0, s88, s81 and rhse;
  - cs uni / closed / thin1 / thin1mb (1 and 2 ranks) / feed1; flz_q, flz_Q and rw x4;
  - RT hcyl, sp3mb, sp3mbmr, thin1, feed, feedx1 (x1line), feedmr (2 ranks), rot1s and csCnr.
- WB fv1/fv/lld/cart: rh.user.hst byte-identical.
- 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL.
- **The feed blocker on e5c54306a / 243146188 is gone.** It was a fresh-run change of 2.8e-3 on the RT feed decks
  (cubed sphere, `clean_face_model = milne`, `force_reference = clean`), and their restart was not bitwise.
  - My bisect put it in the speed2 merge bb9c612c1: 6444c2cbe (mem2) is bitwise, bb9c612c1 changes it. That agrees
    with your cc7b288c5 finding (RHS face fluxes evaluated before CleanCellMilneSetup) and the 7c4805666 fix.
  - On 243146188 the B arm (yres local) happened to restart bitwise on the feed decks; the A arm did not.

## B: b66ff2552 defaults (clean_yres = local) vs dd468ab99
| deck | maxrel final dump | CLEAN lines | metric (B = P16 unless noted) |
|---|---|---|---|
| cc hcyl cbg cb shd nltek1 | 0 | same | scores 1.2061e-02 1.1003 0.45999 1.1805e-02 9.0716e-05 |
| thin_p/r/f, dif, lat, pulse0, s88, s81 | 0 | same | L1 as push 16; pulse0 dE +8.882e-16 |
| cs_uni, cs_closed, cs_thin1, cs_thin1mb (1 and 2 ranks) | 0 | same (thin1mbmr: CLEAN lines differ, data bitwise) | thin1 / thin1mb L1(E) 3.891e-02; closed sum E V rel 0 |
| rhse (hst) | changed | 131/157 | vcs_max 1.85059e-01 (P16 1.85058e-01), it 877 (P16 806) |
| cs_feed1 | 1.6e-6 | 16/19 | dE rel -1.11e-16, dLz 0.0 (P16 -5.84e-17) |
| **flz_q** (sp fv1, chi = 1, rho 1e2, T_gas 1e-4) | **8.98e-1** (m1_e: max\|E - E16\| 0.111 of max 0.124; m1_f 0.8) | 64/76 | dLz 5.87e-16, dE 7.58e-14 (P16 8.40e-14), final Lrz 1.050e-3 (P16 1.086e-3), it 7182 (P16 7909) |
| flz_Q | 2.7e-7 | 73/79 | dLz 1.41e-15, dE 3.19e-14 |
| rw x4 | E 9.8e-11 / 6.2e-11, F1 3.2e-6 / 2.1e-6 (N32 / N64) | changed | odd F1 3.3e-6 / 2.2e-6; shrinks with N |
| RT feed / feedx1 / feedmr | **6.65e-3** (m1_f; m1_e 8.7e-5) | 60/62 | |
| RT sp3mb / sp3mbmr | 2.4e-7 | 46/92, 42/92 | |
| RT csCnr | 8.9e-7 | 56/62 | |
| RT hcyl, thin1, rot1s | 0 | same | |
| WB fv1 / fv / lld / cart | hst differ | | vcs_max printed 8.39985e-08 / 6.83945e-08 / 5.56674e-08 / 1.85442e-07, same as P16 |

- 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL in every B run.
- **flz_q is the one large change.** The radiation field of this cold, radiation-dominated rotating box differs at
  O(1) under yres local, while Lz and E stay conserved to round-off. Its 40-step L_rz changes by 3 %.
  - I have not checked which result is closer to the exact answer.
  - With yres = false (A) it is bitwise to push 16, so the change comes from clean_yres = local alone.
- B vs 243146188 B is bitwise except the RT feed decks (3.1e-3). That is the 7c4805666 fix.

## Opt-in arms (run on 243146188, which has the same opt-in code; maxrel vs the default straight run, 16 cycles)
| deck | fv_lag = extrap | float set (fv_float, pc_factor_float, extrap_share, dt_retry_views = clean) | forced retry (cycle 4, fac 1) + float set, vs float set |
|---|---|---|---|
| sp3 (sp, fv1, force_reference = clean) | 4.0e-5 | 5.1e-8 | 6.4e-6, not bitwise |
| thin1 (cubed sphere, fv1) | 2.7e-12 | 3.9e-9 | 5.7e-16, not bitwise |
| hcyl (Cartesian, clean_rays = fv1) | 6.4e-16 | 1.3e-6 | bitwise |

- All arms ran with rc 0 and 0 CLEAN WARNING; each forced arm shows one `DT_RETRY reason=forced`.
- The forced retry is not bitwise with clean_extrap_share in the set, as you documented. I did not rerun the retry
  without clean_extrap_share.
- **Restarts: bitwise in all 9 opt-in arms**, and for the hcyl fv1 base.
- **clean_yres = local on rhse_cs6** (8 cycles): 0 CLEAN WARNING, 0 DT_RETRY; clean_yres = false is 0 / 0 too.

## Timing
- See NOTE 52611dba2 (unchanged).

## Harness note
- reg/run.sh copies `rw/N*_amp*` per arm, and that glob also matches the earlier arms' copies, so the directories
  grow combinatorially (about 9600 here).
- The fix is `for r in rw/N*_amp*e-3`. The copy in handover_1116_clean/reg has the same line.

## Files (Orion)
- /orion/ptmp/jinma/clean_1115/push16_reg/: build_p17c (b66ff2552), build_p17b, build_p17, bis17/.
- S/reg arms Z17A (A) and Z17B (B); also Y17A/Y17B (243146188) and P17 (e5c54306a); opt-in runs rst/runs/O17_*
  and O18_*.
- S/cmpx_Z17A_DEFF.txt, S/cmpx_Z17B_DEFF.txt, S/cmpx_Z17B_Y17B.txt.
- All .rst dumps were deleted.

-- Orion
