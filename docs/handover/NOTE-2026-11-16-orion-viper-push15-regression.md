# NOTE 2026-11-16 Orion -> viper (USER): CPU regression of push 15 = e3f4881dc vs push 14 = d6b6ed4e2

Answers TASK-2026-11-16-viper-orion-push15-regression.md (7f4c9229d).

**A: not bitwise.** First deck/file: `runs/pulse0`, `sp.m1.00001.bin`. This is a round-off difference that
starts at 37e0e566e; details below. **B: table below.** 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL. **Restart: bitwise**,
4/4 RT for A and 4/4 for B. **Also: CPU speed regression**, 1.2-3.2x per cycle on the clean decks. It comes from
f7273d34d and the memory branch; see the CPU section.

## Setup
- Binaries: Release, gcc/13 + openmpi/4.1, MPI+OpenMP, `-ffp-contract=off`, `Athena_M1_VGD_NMAX=4096`. All runs
  are OpenMP only.
  - push 15 e3f4881dc: md5 `c88ad38674c36c9d2d9bda582090b5f5`
  - push 14 d6b6ed4e2: md5 `e0899b8cb1789f4a5275c18492ee70ca` (the same binary as the push-14 regression)
- Harness: reg/ from handover_1116_clean (29 runs + wb.sh 4 runs + RT 4). Each arm ran on its own p.exclusive
  node, srun -c 112.
- Jobs:
  - 230519: B + push 14, plus a first A attempt that failed (see below);
  - 230521: A;
  - 230522: standalone timing + A determinism;
  - 230526 / 230558: commit bisect (pulse0, wb fv, timing).
- How A passed the key: an athena CLI override of a key that is not in the input file is FATAL:
  `Parameter 'clean_fv_tol' in block 'rad_m1' on command line not found`. So in A a wrapper wrote
  `clean_fv_tol = 1.0e-13` under `<rad_m1>` of every `-i` deck and also appended `rad_m1/clean_fv_tol=1.0e-13`
  on the command line. Restarts pick the key up from the restart file. This covers ladder, gtg, fb, cs, rw, RT
  and wb.
  - Consequence for the comparison: the bin header carries the input text, so the raw .bin bytes differ in every A
    run. For A I therefore compared every variable of every bin data block, raw tab/hst bytes, and every `CLEAN`
    log line.

## A (clean_fv_tol = 1e-13) vs push 14
- **Bitwise on 32 of the 33 run directories.** These are all 29 harness runs except pulse0, plus RT 4/4. In them,
  every bin data block, tab, hst and CLEAN line is identical. Every ana.sh line is identical too (scores, it, warn,
  conservation).
- **runs/pulse0** (gtg pulse0: sp, gd, default fv1 rays, reflecting x1max) is not bitwise.
  - 19 of 20 bins differ, from `sp.m1.00001.bin` onward.
  - The final dump differs by max|diff|/max|group| = 4.6e-14 (m1_f3; m1_e 2e-18).
  - The CLEAN lines and the energy line (+8.882e-16) are identical.
- **wb fv** (rhse_sp, clean_rays=fv) is not bitwise.
  - rh.user.hst and rh.hydro.hst differ from row 4.
  - 32 of 534 CLEAN lines differ, and only in the pic residual value. The sum/pic counts are identical.
  - vcs_max is the same to the printed digits: 6.83945e-08 for both.
  - wb fv1, lld and cart are bitwise.
- **It is deterministic, not thread order.**
  - Two repeat A runs are bitwise to each other and to the regression A run, for both pulse0 and wb fv.
  - With OMP_NUM_THREADS=1, A still differs from push 14 on pulse0.
  - push 14 is reproducible run to run: the pulse0 bins are bitwise to the earlier push-14 regression run.
- **Bisect.** Plain binaries, no key; before 401e274bf the FV-ray stop is the fixed 1e-13. Tested on pulse0 and
  wb fv:
  - 4690abb2d is bitwise to push 14 on both. So the A100 value-copy diff is CPU-bitwise, as before.
  - **37e0e566e is the first non-bitwise commit** ("FV rays: device cell order, wavefront levels, async level
    reductions, r face fluxes of the column sums at once").
  - 37e0e566e, 8701567f5, e472c5e54, f7273d34d and 97ff0a449 are all identical to A on pulse0 and wb fv. So nothing
    later adds a further change on these decks.
  - Likely cause: a different CPU summation or reduction order in the new FV-ray sweep. Your GPU run was bitwise,
    so the reordering seems to show up on the host path only.

## B (new default 1e-10) vs push 14
- maxrel: final dump, max over variables of max|B - P14| / max|P14| of the variable group (vector components share
  the group max).
- CLEAN: number of `CLEAN` lines that differ in any way. "counts" is how many of them differ in sum/pic counts and
  not only in the residual value.

| deck | maxrel final dump | CLEAN lines (counts) | it B / P14 | metric (B = P14 unless noted) |
|---|---|---|---|---|
| cc hcyl cbg cb shd nltek1 | 0 (bitwise) | same | same | scores 1.2061e-02 1.1003 0.45999 1.1805e-02 9.0716e-05 |
| thin_p | 2.6e-16 | 3/299 (2) | 101 / 100 | L1(E) 1.289e-01 L1(F) 1.464e-01 |
| thin_r | 3.7e-16 | 4/299 (3) | 101 / 99 | L1(E) 1.179e-01 L1(F) 1.418e-01 |
| thin_f | 1.3e-15 | 9/299 (1) | 120 / 120 | L1(E) 2.942e-02 L1(F) 6.445e-02 |
| dif | 2.4e-15 | same | 963 | L1(E) 1.536e-02 |
| lat | 4.1e-16 | same | 963 | L1(E) 2.915e-01 |
| pulse0 | 2.3e-11 | same | 58 | dE max +8.882e-16 |
| s88, s81 | 0 (bitwise) | same | same | steady 2.578e-07 |
| rhse (no dumps) | user.hst bitwise; hydro.hst only 2-mom (round-off, 8e-18 abs) | 100/157 (0) | 806 | vcs_max 1.85058e-01 |
| cs_uni, cs_closed | 0 (bitwise) | same | same | closed sum E V rel 0; abs dP/sum abs(F)V 1.58e-01 |
| cs_thin1 | 7.6e-16 | 12/199 (3) | 119 / 117 | **thin1 L1(E) 3.891e-02** (ct.py score thin) |
| cs_thin1mb | 8.8e-16 | 13/199 (2) | 124 / 123 | **thin1mb L1(E) 3.891e-02** |
| cs_feed1 | 2.5e-15 | 16/19 (0) | 257 | dE rel -1.11e-16, dLz total -5.84e-17 |
| flz_q (fv1, chi=1) | **9.0e-6** | 46/76 (42) | 7864 / 7910 | dLz 5.87e-16 (P14 3.91e-16), dE 6.72e-14 |
| flz_Q | 0 (bitwise) | same | 711 | dLz 1.72e-15, dE 3.83e-14 |
| rw x4 | 0 (bitwise) | same | same | |
| RT hcyl / feed / sp3mb / thin1 (16-cycle final) | 0 / 2.7e-15 / 1.0e-13 / 2.4e-12 | 0 / 0 / 23 of 46 (0) / 4 of 31 (1) | | |
| **WB vcs_max** fv1 / fv / lld / cart | | counts differ 30 / 12 / 37 / 0 | | 8.37529e-08 / 6.87442e-08 / 5.57512e-08 / 1.85442e-07, P14 8.39985e-08 / 6.83945e-08 / 5.56674e-08 / 1.85442e-07 (cart bitwise) |

- Our harness has thin1 and thin1mb (2 radial blocks), not thin1l. Both L1(E) are unchanged, and so are A and P14.
- The largest change is flz_q, the fv1 rotating box at 9e-6. Conservation there stays at round-off.
- CLEAN WARNING 0, DT_RETRY 0, FATAL 0 in every B run, including wb and RT. The same holds for A and P14.

## Restart
- `hydro_u same m1 same` for all four RT decks (hcyl, sp3mb, feed, thin1), in A and in B. P14 also 4/4.
- All .rst dumps were deleted afterwards.

## CPU time per cycle (push 15 is slower on CPU)
- Standalone runs, one at a time on an idle node (job 230522). Each value is athena's "cpu time used" divided by
  cycles.

| deck | threads | cycles | P14 s/cycle | A | B | B / P14 |
|---|---|---|---|---|---|---|
| flz_q (sp, fv1) | 8 | 40 | 2.357 | 6.646 | 7.446 | **3.16x** |
| flz_Q (sp, fv1) | 8 | 40 | 0.354 | 0.737 | 0.746 | 2.11x |
| cs_closed | 6 | 4 | 1.849 | 4.975 | 4.633 | 2.51x |
| cs_feed1 | 6 | 10 | 1.330 | 2.853 | 2.668 | 2.01x |
| RT feed deck, 16 cycles (sp, milne) | 8 | 16 | 2.480 | 4.437 | 4.841 | 1.95x |
| cs_thin1 | 6 | 100 | 0.668 | 0.889 | 0.790 | 1.18x |
| thin_f (sp full sphere) | 8 | 150 | 0.557 | 0.618 | 0.600 | 1.08x |

- In the 30-way harness run the arm wall time was also longer: B 310 s, A 285 s, P14 107 s. flz_q alone was 310 s
  vs 107 s.
- Bisect timing: 6 binaries at once on one node, 8+6 threads each, so the load is mild. Values are cpu time in s.

| commit | flz_q | cs_closed |
|---|---|---|
| 4690abb2d | 86.0 | 9.28 |
| 37e0e566e | 80.8 | 9.33 |
| 8701567f5 | 83.7 | 9.10 |
| e472c5e54 | 83.1 | 9.39 |
| f7273d34d | **208.8** | **15.1** |
| 97ff0a449 (memory branch) | **286.8** | **17.9** |

- So the CPU slowdown comes from **f7273d34d** (operator face fluxes once per face, line-solve team size by line
  length), about 2.5x on flz_q and 1.6x on cs_closed. The memory branch 98c731e30..97ff0a449 adds about 1.4x / 1.2x
  more. The cs decks do not use the FV rays, which points at the cell operator / line-solve team sizing on the
  OpenMP host.
- The 1e-10 tolerance itself gives no CPU gain on these decks: B is about the same as A.
- Your q-deck 2.1x speed-up on MI300A is a GPU result. On CPU the same changes cost about 2-3x.

## Files (Orion)
- /orion/ptmp/jinma/clean_1115/push15_reg/: build_new, reg/ (arms TOLA = A, DEFB = B, PREV = P14), wbA/wbB/wbP,
  tools/cmpx.py (data-block compare + maxrel), athenaA (A wrapper), bis/ (bisect builds + worktrees).
- Nothing was pushed to clean-cc.

-- Orion
