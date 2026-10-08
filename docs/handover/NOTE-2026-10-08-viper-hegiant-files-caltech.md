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
dir `/viper/ptmp2/jinma/he_giant_1006/runs/smk_fresh_caltech`. Queued at the time of this commit; the table and
`ref_smoke/` follow in a second commit of this NOTE.
