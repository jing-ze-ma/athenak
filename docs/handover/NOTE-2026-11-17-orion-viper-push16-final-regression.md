# NOTE 2026-11-17 Orion -> viper (USER): final CPU regression of push 16 = 0678599eb vs push 15 = dfe21a746

Answers UPDATE 3 (2911647f6) of TASK-2026-11-16-viper-orion-push16-regression.md. It builds on
NOTE-2026-11-16-orion-viper-push16-regression.md (0aee9a208).

**A: not bitwise. First deck/file: RT `feed`, `cs.m1.00001.bin`.** Only the five runs with `force_reference = clean`
change: RT feed / feedx1 / feedmr, RT sp3mb / sp3mbmr, and the 4 WB runs. These are the 4ac0c0e49 (dt_retry 0002)
reference change, the same in A and B, see below. All other decks are bitwise in data, including multi-rank and
x1line. No harness deck has etotgrav or MHD.

**B: table below.** 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL in every run.

**Restart: bitwise.** 7/7 RT in A and in B.

**CPU:** push 16 is 1.04-1.16x per cycle slower than push 15 and 0.65-0.94x of push 14 (standalone, p.exclusive).

## Setup
- Binaries: Release, gcc/13 + openmpi/4.1, MPI+OpenMP, `-ffp-contract=off`, `Athena_M1_VGD_NMAX=4096`. md5s:
  - 0678599eb: `c69a9fddcaba6b5a7174116f030dbe87`
  - 0aee9a208: `1abe7a4a4aa64b22674e3c94c9cf80ce`
  - push 15 dfe21a746: `0f714f25ea44b2865bed11392edacaa8`
  - push 14: `e0899b8cb1789f4a5275c18492ee70ca`
- Harness: the same as yesterday's NOTE, 33 runs + WB 4 + RT 7. That is reg/ 29 plus:
  - the x1line RT deck `feedx1` (`clean_pc = x1line` in the deck);
  - cs thin1mb with `mpirun -np 2`;
  - RT sp3mb and RT feed with `mpirun -np 2`.
- Arms, run concurrently:
  - OFFD = A: a wrapper writes `clean_gas_exact = false`, `clean_pc_factor = false` and `clean_rot_ring = off` under
    `<rad_m1>` of a copy of each `-i` deck and also passes them on the command line. A command-line key not in the
    deck is FATAL. Restarts take the keys from the restart file.
  - DEFD = B, no keys.
  - Compared against PREV (push 15) and DEFC (0aee9a208), from the same harness copy.
- The compare is per variable on the bin data blocks. Raw .bin bytes differ through the input text in the header,
  so they are not used. Raw tab/hst and every `CLEAN` log line are also compared.
- maxrel = max over variables of max|X - Y| / max|Y| of the variable group; vector components share the group max.
- Jobs:
  - DEFD + OFFD: `srun --overlap -c 56` step in our p.shared allocation 230761 (orion002), 481 s.
  - Timing: 230771 (orion083) and 230834 (orion077), both p.exclusive, standalone, one run at a time.

## A (OFFD) vs push 15
- **Bitwise data in 28 of the 33 runs:** cc ladder, thin_p/r/f, dif, lat, pulse0, s88, s81, rhse, cs uni / closed /
  thin1 / thin1mb / thin1mb on 2 ranks / feed1, flz_q, flz_Q, rw x4, RT hcyl, RT thin1.
- **flz_q:** 5 of 76 CLEAN lines differ while the data stays bitwise. This is the same log-only diagnostic as
  yesterday: push 15 prints a garbage pic value (6.79e+285), from 0a3df762d.
- **Not bitwise: every run with `force_reference = clean`.**

| run | maxrel (final dump) | CLEAN lines |
|---|---|---|
| RT feed, feedx1, feedmr (cubed sphere, milne) | 8.37e-3 (m1_f; m1_e 1.6e-3, hydro_u 3.8e-3) | 60/62 |
| RT sp3mb, sp3mbmr (sp, rh_seed 1e-3) | 1.24e-5 (m1 5.4e-6) | 42/92, 44/92 |
| WB fv1 / fv / lld / cart | rh.user.hst differs | vcs_max below |

  - The diff is already at full size at the first dump (cycle 8) and the same at cycle 16. So it is an offset set
    at the start, not a growth: the reference force now comes from the converged, accepted first coupled solve
    (4ac0c0e49).
  - The CLEAN lines differ from cycle 0 on, for example `pic 12 2.88e-04 solves 19` vs `2.75e-04 solves 18`.
  - OFFD vs OFFC (A on 0aee9a208) differs in exactly the same runs with the same numbers. OFFC was bitwise to
    push 15 on them.
  - The RT feed size, 8e-3, is above the 1.6e-5 you quote for rhse_sp; sp3mb is 1.2e-5. Please confirm 8e-3 is
    expected for the cubed-sphere feed deck.

## B (DEFD = 0678599eb defaults) vs push 15
| deck | maxrel final dump | CLEAN lines | metric (B = P15 unless noted) |
|---|---|---|---|
| cc hcyl cbg cb shd nltek1 | 0 (bitwise) | same | scores 1.2061e-02 1.1003 0.45999 1.1805e-02 9.0716e-05 |
| thin_p / thin_r / thin_f | 0 | same | L1(E) 1.289e-01 / 1.179e-01 / 2.942e-02 |
| dif, lat, pulse0 | 0 | same | L1(E) 1.536e-02, 2.915e-01; dE +8.882e-16 |
| s88, s81 | 0 | same | steady 2.578e-07 |
| rhse (hst only) | 1.3e-17 | 107/157 | vcs_max 1.85058e-01, it 806 |
| cs_uni, cs_closed | 0 | same | closed sum E V rel 0, \|dP\|/sum\|F\|V 1.58e-01 |
| cs_thin1, cs_thin1mb (1 and 2 ranks) | **0 (bitwise again: fv_joint off)** | same | **thin1 / thin1mb L1(E) 3.891e-02** |
| cs_feed1 | 1.2e-13 | 13/19 | dE rel 0.0 (P15 -1.11e-16), dLz -5.84e-17 |
| flz_q | 2.0e-5 | 36/76 | dLz 5.87e-16, dE 8.40e-14 (P15 6.72e-14), it 7909 (P15 7864) |
| flz_Q | 5.5e-14 | 66/79 | dLz 1.25e-15 (P15 1.72e-15), dE 2.56e-14 (P15 3.83e-14) |
| rw x4 | E 1.2e-10 / 6.4e-11, F1 3.9e-6 / 4.2e-6 (N32 / N64) | changed | odd part F1 3.7e-6 / 4.3e-6 |
| RT feed / feedx1 / feedmr | 8.37e-3 | 60/62 | force_reference = clean, see A |
| RT sp3mb / sp3mbmr | 1.24e-5 | 48/92, 44/92 | force_reference = clean |
| RT thin1, RT hcyl | 0 | same | |
| **WB vcs_max** fv1 / fv / lld / cart | hst differs | | B 8.39985e-08 / 6.83945e-08 / 5.56674e-08 / 1.85442e-07; P15 8.39985e-08 / 6.87442e-08 / 5.57512e-08 / 1.85442e-07; A 8.39985e-08 / 6.88374e-08 / 5.57380e-08 / 1.85468e-07 |

- CLEAN WARNING 0, DT_RETRY 0 and FATAL 0 in every run of A and B, including WB, RT and multi-rank.

## B (0678599eb) vs 0aee9a208 (DEFC): which decks change and why
- **fv_joint back to opt-in:** cs_thin1 (7.6e-16) and cs_thin1mb (1.4e-15) change, and so does RT thin1 (2.4e-12).
  They are now bitwise to push 15 again.
- **force_reference = clean (4ac0c0e49, dt_retry 0002):** RT feed / feedx1 / feedmr 8.37e-3, RT sp3mb / sp3mbmr
  1.24e-5, and the WB hst.
  - The WB rh.user.hst files differ, but vcs_max is the same to the printed digits in B.
  - In A, vcs_max moves: fv 6.874 -> 6.884e-8, lld 5.575 -> 5.574e-8, cart 1.85442 -> 1.85468e-7.
- **cs_feed1:** 4.5e-15. This is the only other change; it is round-off.
- **Unchanged (bitwise to DEFC):** all others, including the steady-start decks s88, s81 and rhse (`clean_init =
  steady`). So the steady-start predict fix 5e33f056a and the init audit (4cacdae21, 8db6f6141) leave every harness
  deck bitwise.
  - Our steady-start decks use the sph closure (s88/s81) or one steady start with force_reference = none (rhse).
    There is no steady + predict deck with force_reference = clean in the harness, so the 1.6e-5 rhse_sp case is not
    in it. The nearest is RT sp3mb at 1.24e-5.
- **d320347b3** (dt_retry checkpoint restores cln_fvi / cln_fref_set): no run had a DT_RETRY, so it is not exercised.

## Restart
- `hydro_u same m1 same` for hcyl, sp3mb, sp3mbmr (2 ranks), thin1, feed, feedx1 (x1line) and feedmr (2 ranks), in A
  and in B.
- No harness deck rotates the vet_gd FV rays, so the ring default (live) is never switched on.
- All .rst dumps were deleted.

## CPU time per cycle (standalone, p.exclusive, one run at a time)
- Each value is athena "cpu time used" / cycles, in s.
- Push 16 and push 15 ran on the same node (orion077, job 230834).
- Push 14 is from job 230771 (orion083). In that job push 15 measured within 0-6 % of the orion077 values on every
  deck except cs_closed, 1.545 vs 1.394 (11 %).

| deck | threads | cycles | P14 | P15 | P16 0678599eb | P16 / P15 | P16 / P14 |
|---|---|---|---|---|---|---|---|
| flz_q (sp, fv1) | 8 | 40 | 2.359 | 1.680 | 1.878 | 1.12 | 0.80 |
| flz_Q | 8 | 40 | 0.362 | 0.264 | 0.277 | 1.05 | 0.77 |
| cs_closed | 6 | 4 | 1.910 | 1.394 | 1.572 | 1.13 | 0.82 |
| cs_feed1 | 6 | 10 | 1.309 | 1.026 | 1.185 | 1.15 | 0.91 |
| RT feed deck, 16 cycles | 8 | 16 | 1.734 | 1.301 | 1.466 | 1.13 | 0.85 |
| cs_thin1 | 6 | 100 | 0.660 | 0.558 | 0.581 | 1.04 | 0.88 |
| thin_f (sp full sphere) | 8 | 150 | 0.538 | 0.485 | 0.508 | 1.05 | 0.94 |
| x1line feed deck, 16 cycles | 8 | 16 | 1.764 | 1.045 | 1.147 | 1.10 | 0.65 |

- **Push 15 dfe21a746 already recovered the push-15 slowdown.** Last time, push 15 e3f4881dc was 1.2-3.2x slower than
  push 14. The host line solve with one thread per line makes dfe21a746 0.6-0.9x of push 14.
- **Push 16 costs 4-16 % per cycle on CPU over dfe21a746 and stays faster than push 14 on every deck.** The largest
  costs are on gas-coupled decks: cs_feed1 1.15x, cs_closed 1.13x, RT feed 1.13x, flz_q 1.12x. That points at
  clean_gas_exact (default on).
- For reference, 0aee9a208 (fv_joint on) on orion083: cs_thin1 0.634, so 1.17x push 15. With fv_joint off it is now
  1.04x.

## Files (Orion)
- /orion/ptmp/jinma/clean_1115/push16_reg/: build_p16d (0678599eb), S/reg (arms PREV, OFFD, DEFD, DEFC, ...).
- S/cmpx_{OFFD_PREV,DEFD_PREV,DEFD_DEFC,OFFD_OFFC}.txt, S/ana_DEFD.txt, S/athenaD (A wrapper), S/tools.
- Nothing was pushed to clean-cc or clean-cc-push16.

-- Orion
