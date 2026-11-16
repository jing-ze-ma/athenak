# NOTE 2026-11-17 Orion -> viper (USER): race3 a9f481263 CPU check

Answers request 2 of NOTE-2026-11-17-viper-orion-bugfixes-branches.md.

**Result: a9f481263 is bitwise to dd468ab99 on the whole reg harness and on all 10 restarts. The extra hydro and MHD
nscalars decks are bitwise too. The one exception is mhd_compress_shwave: there dd468ab99 is itself not reproducible
with OpenMP threads on CPU, and a9f481263 fixes it.** So E2 is not CPU-neutral for shearing-box runs with more than
one thread. Every other fix is CPU-neutral.

## Binaries (Orion, Release, gcc/13 + openmpi/4.1, MPI+OpenMP, -ffp-contract=off, Athena_M1_VGD_NMAX=4096)
- a9f481263 (clean-cc-race3-viper): md5 `9a15ab9dd2793157a76155aef477f13e`, built from `git archive`.
- dd468ab99 (push 16): md5 `824b439cbb0bcdf6d241a781357be5d7`, the same binary and harness results as
  NOTE-2026-11-17-orion-viper-push16-dd468ab99.
- Everything ran in `srun --overlap -c 112 --cpu-bind=none` steps of our p.shared allocation 230761 (orion002).
  The harness arm took 221 s and the whole step 857 s.

## reg harness (30 runs + WB 4 + 10 restarts), a9f481263 vs dd468ab99, push-16 defaults
- **All run directories: bin data blocks, tab/hst and every CLEAN line identical.**
  - Decks: cc ladder; thin_p/r/f, dif, lat, pulse0, s88, s81, rhse; cs uni/closed/thin1/thin1mb (1 and 2
    ranks)/feed1; flz_q, flz_Q, rw x4.
  - Also the straight runs of the RT decks, rot1s and csCnr.
  - rhse hst files are byte-identical to dd468ab99. The duplicate copy in rhse/DEFF's hst comes from the earlier
    stopped launch.
- **WB fv1/fv/lld/cart: rh.user.hst byte-identical.** vcs_max 8.39985e-08 / 6.83945e-08 / 5.56674e-08 /
  1.85442e-07.
- **Restart: `hydro_u same m1 same` on all 10.** These are hcyl, sp3mb, sp3mbmr (2 ranks), thin1, feed,
  feedx1 (x1line), feedmr (2 ranks), rot1s (ring live, 20 vs 10 + rst + 10) and csCnr (16 vs 8 + rst + 8).
- 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL.
- The sp hydro decks of the harness (rhse, wb fv1/fv/lld, thin*, sp3mb, flz) run with fofc off (default) and no
  sp_face_avg, which is the R4 exposure. They are all bitwise.

## Extra decks for the fixes (a9f481263 vs dd468ab99, 8 threads)
| deck | what it exercises | result |
|---|---|---|
| adv3_hyd: tests/advect_hyd made 3-D (64x32x32, MeshBlocks 32x16x16), nscalars = 2, flow_dir = 2, dynamic, hlle, 40 cycles | R1/R2, hydro x2/x3 scalar flux sweeps | all bins byte-identical |
| adv3_mhd: the same with `<mhd>`, nscalars = 2 | R1, MHD x2/x3 scalar flux sweeps | all bins byte-identical |
| mhd_compress_shwave (inputs/shearing_box, 32^3, 100 cycles, bin output) | E2, shearing box | **differs**, see below |

## mhd_compress_shwave: dd468ab99 is not reproducible on CPU, a9f481263 is
Run-to-run check on the cycle-100 dump (`mhd_w_bcc.00002.bin`):

| binary | 1 thread (2 runs) | 8 threads (3 runs) |
|---|---|---|
| a9f481263 | identical | identical, and the same as 1 thread |
| dd468ab99 | identical, and **the same as a9f481263** | **2 of 3 runs differ from each other** |

- In the 8-thread runs that differ, the differences are 8e-3 and 1e-1 relative to the variable's own max. The 1e-1
  is on bcc3, whose max is 7e-10; the absolute differences are about 3e-10 to 8e-10.
- So the one-past-the-end write in `par_for_inner(0, nj)` also races between OpenMP threads on the host.
- With E2, every 8-thread run equals the 1-thread result, which the old code also gave serially.
- **E2 is therefore a CPU bug fix as well,** for shearing-box runs with OpenMP threads. The CPU-neutral claim holds
  for R1-R4.

## Files (Orion)
- /orion/ptmp/jinma/clean_1115/push16_reg/: build_r3 (from src_r3).
- S/reg arm R3; S/cmpx_R3_DEFF.txt.
- S/x3: extra decks adv3_hyd, adv3_mhd and shwave; runs run_R3, run_DF and rep_*.
- All .rst dumps were deleted.

-- Orion
