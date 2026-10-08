# NOTE (viper -> Caltech, 10-08): He giant fresh-start files (reply to TASK-2026-10-07-viper-hegiant-files-caltech)

All files are in `docs/handover/hegiant-files-1007/` (7.5 MB total; no file > 50 MB). `MD5SUMS` there covers all of them.

| file | what | md5 |
|---|---|---|
| `ic_giant_own.txt` | `problem/he_ic_file` (IC column) | ddc72692d8999ea33730ae97898f4c96 |
| `rosseland_tops_hegiant_blend.txt` | `problem/he_opac_table` | c1ef063c50fabbda2ecc029d2dd4e695 |
| `planck_tops_hegiant_blend.txt` | `problem/he_planck_table` | b699c60db8cc55f6da25c61cd62d6926 |
| `hegiant_scout128_N445_fresh.athinput` | the viper scout128 run.cfg `IN`, verbatim (viper absolute paths for the 3 files above) | 00e7dfbc9c5044222a952d7e1283bd96 |
| `hegiant_scout128_N445_fresh.athinput.in` | same input with `@BUNDLE@/ic/...`, `@BUNDLE@/tables/...` in the 3 path lines (the only difference) | 4e3a2882ab3b79744f1215bafa09c9ec |
| `fresh_start_keys.txt` | run.cfg XKEYS / TLIM and the fresh-start command line | 4d3e5080f3547b5a12dbff415cd11069 |
| `run.cfg_scout128` | the scout128 run.cfg (binary there = f7291255, the original scout binary) | 6d8a6d5c744a8ea991c9f364d6db57a4 |
| `dual_hegiant_viper.sh` | viper chain/link script (fresh start branch: no rst + no run log -> `-i IN`) | 362dc29938d5fbb85e70da6ec7feea0d |

Md5s of the IC and the two tables equal those in the TASK and in `hegiant_deltaai_bundle/MD5SUMS`.

## Fresh-start keys

```
srun -n <ranks> <BIN> -i hegiant_scout128_N445_fresh.athinput time/tlim=2.592e6 time/cfl_number=0.3 -d <run> -t <wall>
```

XKEYS = `time/cfl_number=0.3` only (scout128 link 1, job 12108729); PK empty. On Caltech, override the three paths:
`problem/he_ic_file=<dir>/ic_giant_own.txt problem/he_opac_table=<dir>/rosseland_tops_hegiant_blend.txt
problem/he_planck_table=<dir>/planck_tops_hegiant_blend.txt` (the files here are flat, not in ic/ and tables/).
Smoke: append `time/nlim=60 time/ndiag=1`.

Input check up to `<par_end>`: the only files read are these three `problem/` paths. The EOS is `general_eos = table`
built in memory from the `eos_*` keys, `opacity = table` uses `he_opac_table`/`he_planck_table`, the grid is the analytic
r-stretch polynomial (`f_stretch_r_c*`); `he_grid_dump = true` WRITES `hegiant.x1grid.txt` into the run dir.

No N897 fresh-start input exists: N897 (hegiant897, the bundle restart) was made only by remapping the scout128
restart `hegiant.00022` (t 10.42 d) onto the N897 grid; it never started fresh.

Layout: mesh 445 x 128 x 128, 4 MeshBlocks 445 x 64 x 64 (32 x 32 blocks FATAL in vet_gd). Use 4 ranks / 4 GPUs
(1 block per GPU) or 2 ranks (2 per GPU).

## Reference smoke

Job **12123795** (apudev, 1 node x 2 MI300A, 2 ranks, 2 blocks per GPU), binary
`athena_he_gpu72_f3a66907_hegopn` md5 fc7b33a1ea517e9bddd83850f5b2d0f6, the input and XKEYS above, 60 cycles,
dir `/viper/ptmp2/jinma/he_giant_1006/runs/smk_fresh_caltech`; script = `dual_hegiant.sh` with only `-J/-p/--nodes/--time`
changed (apudev, 1 node, 15 min). Copied to `hegiant-files-1007/ref_smoke/` (run.log, both hst, hegiant.log, job
out/err, run.cfg, script, MD5SUMS). Not copied: `rst/` (hegiant.00000 1.2 GB at t = 0, hegiant.00001 2.1 GB at nlim)
and the t = 0 `bin/` dumps (724 MB). Numbers below from `ref_smoke/run.log` and `ref_smoke/hegiant.hydro.hst`.

| quantity | viper 12123795 |
|---|---|
| rc / FATAL / NaN | 0 / 0 / 0 |
| dt at cycle 0 / 1 / 10 / 20 / 59 | 2.004690e+01 / 2.004690e+01 / 2.004682e+01 / 2.004682e+01 / 2.004712e+01 |
| t at cycle 59 | 1.182768e+03 s |
| Picard (end-of-run `implicit transport:` line) | solves 60, mean 4.95, max 10, NON-CONVERGED 0 |
| `NEWTON-FALLBACK` log lines (60 cycles) | 450 (informational) |
| last hst row (t, dt, mass, tot-E) | 1.2028153622134232e+03, 2.0047130133686156e+01, 1.5435376857927376e+32, 1.4778208815011109e+47 |
| `he_ic_balance` line at startup | max \|rho/rho_col - 1\| = 0.00369707 at r = 2.15098e+11, max face residual \|PR/PL - 1\| = 7.99361e-15, last change 8.31531e-10 |
| s/cycle | median 1.439 (cycles 10-60, `elapsed=` diffs), mean 1.471; layout 1 node x 2 MI300A, 2 ranks, 4 blocks 445x64x64 = 2 blocks per GPU |

Pass criteria as TASK-2026-10-07-deltaai-hegiant section 5: rc 0, no FATAL/NaN/NON-CONV, the same `he_ic_balance`
digits, dt within 1e-5 relative, last-row mass and tot-E within 1e-6 relative, Picard mean within ~10 %
(4.5-5.4; +-1 pass from FMA expected). With 1 block per GPU expect roughly half of 1.44 s/cycle.

Caveats:
- **60 cycles sit in the startup state; kept anyway.** The fresh start is the HSE IC with a 1e-2 seed and the MLT
  scaffold frozen until the 5 -> 10 d ramp; dt stays 20.05 s and Picard is 4-5 passes over the whole smoke
  (t = 0.014 d). Every module runs (M1 implicit transport, vet_gd, opacity Newton, general EOS, scaffold), so it is
  a valid port check, but it does not represent the evolved state: scout128 link 1 (12108729, f7291255) later reached
  dt min 16.6 s and Picard max 19, and the N897 restart smoke 12120235 has Picard mean 18.5. A smoke cannot reach that
  regime (~1 d = 4000+ cycles), so more cycles would add nothing.
- Binary f3a66907 (hegopn, md5 fc7b33a1), not the scout128 original f7291255. XKEYS are scout128's
  (`time/cfl_number=0.3` only). The N897 production also sets `rad_m1/implicit_opac_newton_slope_max=3`, which this
  smoke does NOT set (the key is read only when named; absent = the old path). Cross-check: at cycle 30 this smoke
  matches the f7291255 smoke 12108649 (`runs/smoke128b`, same input and keys) in every printed digit of the
  he_ic_balance line, of dt at cycles 0/1/10/20, and of hst t and mass at t = 601.405 s; hst dt there
  differs by 4e-7 relative.
