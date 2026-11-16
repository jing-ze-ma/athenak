# NOTE 2026-11-16 Orion -> viper (USER): CPU regression of push 16 = 0aee9a208 vs push 15 = dfe21a746

Answers TASK-2026-11-16-viper-orion-push16-regression.md, including UPDATE (a97e368cb) and UPDATE 2 (0aee9a208).

**A: bitwise.** On 0aee9a208 and on a97e368cb, every bin data block, tab and hst is identical to push 15 on every deck,
including multi-rank, x1line and restart. The one log-only difference is in flz_q, see below.
No harness deck has etotgrav or MHD, so there is no separate list.
**B (0aee9a208, no keys): table below.** 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL in every run.
**Restart: bitwise.** 7/7 RT in every arm (A and B, both heads).
**a97e368cb vs 51e2ce016: every file is bitwise**, as expected.
**rhse_cs6 / cs24 with clean_yres = true: the CPU warns too** (1 CLEAN WARNING + 1 DT_RETRY each, cycle 0).
**CPU timing: not yet.** The standalone timing job is still queued; a follow-up NOTE will bring it.

## Setup
- Binaries: Release, gcc/13 + openmpi/4.1, MPI+OpenMP, `-ffp-contract=off`, `Athena_M1_VGD_NMAX=4096`. md5s:
  - 0aee9a208 (push-16 head): `1abe7a4a4aa64b22674e3c94c9cf80ce`
  - a97e368cb: `07fe959ecfc5352b9948e3b9501d0e62`
  - 51e2ce016: `59aeb3e8b04bb1925f96ded80e524ee4`
  - push 15 dfe21a746: `0f714f25ea44b2865bed11392edacaa8`
- Harness: reg/ from handover_1116_clean (29 runs + wb.sh 4 + RT 4), plus 4 added runs:
  - **x1line deck:** the RT feed deck (cubed sphere, fv1, milne) with `clean_pc = x1line` in the deck, run as RT
    `feedx1` (16 cycles vs 8 + restart + 8) in every arm;
  - **multi-rank:** cs thin1mb (2 radial blocks) under `mpirun -np 2`, and RT sp3mb and RT feed under `mpirun -np 2`.
- Arms:
  - PREV = push 15;
  - OFFA / OFFC = arm A on a97e368cb / 0aee9a208;
  - DEFB = a97e368cb defaults (clean_yres on);
  - DEF0 = 51e2ce016 defaults;
  - DEFC = 0aee9a208 defaults (B).
- How arm A passes the keys: an athena command-line override of a key that is not in the input is FATAL. So a wrapper
  writes the five keys (`clean_gas_exact = false`, `clean_yres = false`, `clean_pc_factor = false`,
  `clean_fv_joint = false`, `clean_rot_ring = off`) under `<rad_m1>` of a copy of every `-i` deck. It also passes
  them on the command line. Restarts take the keys from the restart file.
  - Raw .bin bytes therefore differ through the input text in the header. The compare is per variable on the bin
    data blocks, plus raw tab/hst and every `CLEAN` log line.
- Jobs:
  - 230765 (p.exclusive orion097): PREV, OFFA, DEFB and DEF0, concurrently, 333 s.
  - 0aee9a208 arms DEFC and OFFC, and the cs6/cs24 check: a 56-core `srun --overlap` step 230761.4 inside our
    p.shared allocation on orion002, 312 s.
  - Timing: job 230771, queued.

## A vs push 15
- **OFFC (0aee9a208) and OFFA (a97e368cb):** all 33 run dirs are bitwise in data, as are WB fv1/fv/lld/cart (hst and
  CLEAN lines) and RT 7/7. The ana.sh scores, iterations and conservation lines are identical. OFFC vs OFFA is
  byte-identical everywhere.
- **The only difference is in the log of flz_q (fv1, chi = 1): 5 of 76 CLEAN lines.** The data stays bitwise.
  - Example, cycle 1: push 15 prints `pic 12 6.79e+285` (garbage), A prints `pic 12 3.95e-11`.
  - In cycles 3 and 32, `it` is 1 vs 0 and `solves` is 54/52 vs 56.
  - This looks like the always-on relative Newton change, 0a3df762d: the old printout divided by about |E_old| = 0.
    I did not bisect it.

## B (0aee9a208 defaults) vs push 15
- maxrel: final dump, max over variables of max|B - P15| / max|P15| of the variable group.
- CLEAN: number of `CLEAN` lines that differ.
- B = P15 unless noted.

| deck | maxrel final dump | CLEAN lines | metric |
|---|---|---|---|
| cc hcyl cbg cb shd nltek1 | 0 (bitwise) | same | scores 1.2061e-02 1.1003 0.45999 1.1805e-02 9.0716e-05 |
| thin_p / thin_r / thin_f | 0 | same | L1(E) 1.289e-01 / 1.179e-01 / 2.942e-02 |
| dif, lat | 0 | same | L1(E) 1.536e-02, 2.915e-01 |
| pulse0 | 0 | same | dE +8.882e-16 |
| s88, s81 | 0 | same | steady 2.578e-07 |
| rhse (hst only) | 1.3e-17 | 107/157 | vcs_max 1.85058e-01, it 806 |
| cs_uni | 0 | same | exact |
| cs_closed | 0 | same | sum E V rel 0, \|dP\|/sum\|F\|V 1.58e-01 |
| cs_thin1 | 7.6e-16 | 12/199 | **thin1 L1(E) 3.891e-02**, it 118 (P15 119) |
| cs_thin1mb | 1.4e-15 | 13/199 | **thin1mb L1(E) 3.891e-02**, it 122 (P15 124) |
| cs_thin1mb, 2 ranks | 0 | same | L1(E) 3.891e-02 |
| cs_feed1 | 1.1e-13 | 16/19 | dE rel -1.11e-16, dLz -5.84e-17 |
| flz_q | **2.0e-5** | 36/76 | dLz 5.87e-16, dE 8.40e-14 (P15 6.72e-14), it 7909 (P15 7864) |
| flz_Q | 5.5e-14 | 66/79 | dLz 1.25e-15 (P15 1.72e-15), dE 2.56e-14 (P15 3.83e-14) |
| rw x4 | E 1.2e-10 / 6.4e-11, F1 3.9e-6 / 4.2e-6 (N32 / N64) | changed | odd part F1 3.7e-6 / 4.3e-6 |
| RT feed / feedmr / **feedx1** (16 cycles) | 1.05e-13 / 1.05e-13 / 1.15e-13 | same | |
| RT sp3mb / sp3mbmr | 6.3e-10 / 6.3e-10 | 42/92, 38/92 | |
| RT thin1 / hcyl | 2.4e-12 / 0 | 8/62 / same | |
| **WB vcs_max** fv1 / fv / lld / cart | | | 8.39985e-08 / 6.83945e-08 / 5.56674e-08 / 1.85442e-07; P15 8.39985e-08 / 6.87442e-08 / 5.57512e-08 / 1.85442e-07 |

- **No warnings.** CLEAN WARNING 0, DT_RETRY 0 and FATAL 0 in every B run, including WB, RT and multi-rank. The same
  holds for A, P15, DEFB and DEF0.
- **rw is not an N64 shrink.** The rw changes are at the 1e-6 level of the flux and do not shrink from N32 to N64
  (F1 3.9e-6 -> 4.2e-6). With clean_yres on (a97e368cb) they were 7.8e-6 -> 5.7e-6. cmpx's raw tab maxrel for rw
  (up to about 0.8) comes from a column at round-off level (about 1e-11 of the flux), not from the physics.
- **Decks that change** are those with gas coupling, whether from gas_exact or a Newton-path change: flz, feed, rhse,
  rw, sp3mb, thin1. The v = 0 radiation-only decks and the cc ladder stay bitwise.

## a97e368cb (clean_yres on) vs 0aee9a208
- Bitwise on: cc, the thin/dif/lat/pulse0/s88/s81 decks, cs_uni/closed/thin1/thin1mb (1 and 2 ranks), RT hcyl and
  RT thin1. These are all decks where clean_yres does not act.
- These change:

| deck | maxrel |
|---|---|
| flz_q | 4.5e-8 |
| flz_Q | 1.3e-14 |
| cs_feed1 | 5.7e-13 |
| RT feed / feedmr / feedx1 | 7e-15 / 7e-15 / 5.7e-14 |
| RT sp3mb, sp3mbmr | 2.0e-9 |
| rhse hst | 1.7e-17 |
| rw N64 | 2.4e-5 |

- With clean_yres on, B vs P15 had flz_q 2.0e-5, sp3mb 2.0e-9, cs_feed1 5.7e-13, and the same zero
  warnings/retries on all harness decks.

## clean_yres on the GPU failing decks, on CPU (0aee9a208, 8 cycles, 4 threads)
| deck | clean_yres = true | false |
|---|---|---|
| g4_rhse_cs6 | 1 CLEAN WARNING + 1 DT_RETRY | 0 / 0 |
| g4_rhse_cs24 | 1 CLEAN WARNING + 1 DT_RETRY | 0 / 0 |

- The true case fails at cycle 0, t = 0, in the first backward-Euler solve, with the same line in both decks:
  `CLEAN WARNING cycle 0: backward-Euler solve not converged (dt_retry, 1 so far)` and
  `DT_RETRY ... dt 1.530931e-02 -> 7.654655e-03 reason=NON-CONVERGED`.
- So it is not GPU-specific, and opt-in is the right default.

## Restart
- `hydro_u same m1 same` for hcyl, sp3mb, sp3mbmr (2 ranks), thin1, feed, feedx1 (x1line) and feedmr (2 ranks).
- This holds in OFFA, OFFC, DEFB, DEFC, DEF0 and P15. The predictor-increment WARNING text there is stale, as before.
- No harness deck rotates the vet_gd FV rays, so the live ring is never switched on. I did not test the case of a
  push-15 restart file in a rotating deck.
- All .rst dumps were deleted afterwards (378 files).

## CPU time
- The standalone timing (job 230771) is still queued: 0aee9a208, then push 15, then push 14, one run at a time on
  one node, on flz_q, flz_Q, cs_closed, cs_feed1, cs_thin1, thin_f, RT feed 16 cycles and feedx1 16 cycles.
- It is not run on a shared node, because timing there is not valid. Values will follow in a separate NOTE.
- The cpu times in the regression logs are not usable as timing, because the arms ran concurrently.

## Files (Orion)
- /orion/ptmp/jinma/clean_1115/push16_reg/: build_p16c / build_p16b / build_p16 / build_p15.
- S/reg (arms PREV, OFFA, OFFC, DEFB, DEFC, DEF0).
- S/cmpx_*.txt, S/ana_*.txt, S/tools (cmpx.py, rwcmp.py, tim.sh), S/athenaA and S/athenaC (A wrappers).
- yres/ (cs6/cs24 on/off).
- Nothing was pushed to clean-cc or clean-cc-push16.

-- Orion
