# NOTE 2026-11-17 Orion -> viper (USER): dd468ab99 CPU check, and the rotating mid-run restart deck in reg/

This answers UPDATE 4 (00164a70e) of TASK-2026-11-16-viper-orion-push16-regression.md. It also answers item 2 of
NOTE-2026-11-17-viper-orion-bugs-races.md: reg/ now has a rotating deck with a restart in the middle of the run.

**Fresh runs: bitwise.** dd468ab99 is identical to 0678599eb on the full reg harness (data, tab/hst, every CLEAN
line, WB). This includes multi-rank, x1line, the new rot1s deck and csCnr.
**Restart: bitwise on dd468ab99.** 10/10 RT: the 8 RT decks, rot1s and csCnr.
**0678599eb is not bitwise on restart** for rot1s (ring = live) and for csCnr. push 15, ring = off and 0a03129a8
are bitwise.
0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL in every run.

## Binaries
- All builds: Orion, Release, gcc/13 + openmpi/4.1, MPI+OpenMP, `-ffp-contract=off`, `Athena_M1_VGD_NMAX=4096`.
  All runs are on Orion; nothing ran on Raven.
- md5s:
  - dd468ab99: `824b439cbb0bcdf6d241a781357be5d7` (built from `git archive dd468ab99`; ptmp was too slow for a
    worktree)
  - 0a03129a8: `1b19775763fa0b92290ce91fd5dd7d8f`
  - 0678599eb: `c69a9fddcaba6b5a7174116f030dbe87` (the binary from the final regression)
  - push 15 dfe21a746: `0f714f25ea44b2865bed11392edacaa8`
- The dd468ab99 arm ran as an `srun --overlap -c 112` step in our p.shared allocation 230761 (orion002): 203 s for
  the harness, WB, 10 RT and csCnr on 0678599eb.
  - `-c 224` on 230761 and `-c 128` on 230775 are refused ("More processors requested than permitted";
    CPUs/Task 112 and 64).

## New harness deck: rst/rot1s.athinput + rst/rtn.sh (committed to handover_1116_clean/reg)
- **rot1s** is rst/thin1.athinput (cubed sphere 64x8x8, gd fv1 rays, hesdirk2) with these rotation keys:
  `vet_gd_rotate_dt = 2.25` (a new set every 4 cycles), `vet_rotate_avg = 4`, `vet_rotate_advance = true` and
  `vet_gd_rotate_xfade = true`.
  - With these, ring_ok holds and push 16 takes `clean_rot_ring = live` by default.
  - Bins are written at cycle 20 and a restart at cycle 10.
- **rtn.sh NAME BIN DECK THREADS NTOT NRST** runs NTOT cycles straight vs NRST + restart + the rest, then
  cmpbin.py on the final bins. It is a generalised rt.sh.
- cmds.txt gets `./rtn.sh ${ARM}_rot1s $BIN $REG/rst/rot1s.athinput 6 20 10`. The README gate table and the run
  count (30) are updated.

## Rotating restart results (Orion CPU, 6 threads)

| binary | ring | 20 straight vs 10 + rst + 10 | 100 straight vs 52 + rst + 48 |
|---|---|---|---|
| push 15 dfe21a746 | (none) | bitwise | bitwise |
| 0678599eb | live (default) | **not bitwise**: m1_e 9.1e-6, m1_f 1.4e-5 (group max) | **not bitwise**: m1_e 5.7e-16, m1_f 4.3e-16 (f2/f3 own max 1.7e-14) |
| 0678599eb | off | bitwise | bitwise |
| 0a03129a8 | live | bitwise | bitwise |
| 0a03129a8 | off | bitwise | bitwise |
| dd468ab99 (harness arm) | live | bitwise | |

- **The ring is really live with defaults.** On fresh runs, 0678599eb with defaults vs ring = off differs by
  m1_e 4.6e-5 at cycle 100 and 1.3e-4 at cycle 20.
- **The fix does not change fresh runs.**
  - 0a03129a8 straight is bitwise to 0678599eb straight, at both lengths.
  - 0678599eb ring = off is bitwise to push 15.
- **The 0678599eb restart error depends on where the restart falls in the rotation.** A restart at 10 gives about
  1e-5, a restart at 52 gives about 1e-14. The latter is your A100 2.5e-14 case.
  - No crash on CPU.
  - The hydro_u dumps are bitwise in all cases.
- The predictor-increment WARNING after the restart is the stale text seen in every run.

## dd468ab99 vs 0678599eb, fresh runs (full harness, defaults)
- Every run dir has bitwise bin data blocks, tab/hst and identical CLEAN lines. That covers:
  - the cc ladder; thin_p/r/f, dif, lat, pulse0, s88, s81 and rhse;
  - cs uni / closed / thin1 / thin1mb (1 and 2 ranks) / feed1; flz_q, flz_Q and rw x4;
  - the straight runs of all RT decks, including rot1s and csCnr.
- The WB rh.user.hst files are byte-identical for fv1/fv/lld/cart. vcs_max is 8.39985e-08 / 6.83945e-08 /
  5.56674e-08 / 1.85442e-07.
- One harness note: rhse/DEFF's hst holds two identical copies of the run. An interrupted first launch appended to
  it; each copy is byte-identical to 0678599eb.

## Restarts on dd468ab99 (all `hydro_u same m1 same`)
- The 8 RT decks: hcyl, sp3mb, sp3mbmr (2 ranks), thin1, feed, feedx1 (x1line) and feedmr (2 ranks).
- rot1s (ring live).
- **csCnr** (/orion/ptmp/jinma/clean_1115/audit2/rst/csCnr.in: cubed sphere, cell layout, no clean_rays named, 16
  vs 8 + restart + 8): bitwise.
  - On 0678599eb csCnr is **not bitwise**: hydro_u 7.6e-4 abs, m1 8.5e-4 abs, grouped maxrel 8.2e-4.
  - The restart logs `### WARNING: restart file has no <rad_m1> hesdirk2 slope; the first step is backward
    Euler`, the same channel-count misread as the ring.
  - c40170eab / 665ae99de fix it.
- csSnr (staggered) was not run: dd468ab99 has the cell layout only.

## Files (Orion)
- /orion/ptmp/jinma/clean_1115/push16_reg/: build_p16f (dd468ab99, from src_p16f) and build_p16e (0a03129a8).
- S/reg: arm DEFF; rst/runs/ROT_* (100/52) and ROTs_* (20/10); DEFF_csCnr and DEFD_csCnr.
- S/cmpx_DEFF_DEFD.txt.
- All .rst dumps were deleted.

-- Orion
