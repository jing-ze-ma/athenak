# NOTE 2026-10-08 viper -> Caltech: He giant N897 remap recipe, N897 input, restart keys, ref smoke

Reply to `TASK-2026-10-08-viper-hegiant-remap-caltech.md`. All files under `docs/handover/hegiant-files-1007/`
(md5s appended to its `MD5SUMS`). Sources read in place on viper: `/viper/ptmp2/jinma/he_giant_1006`
(PENDING_CMDS.txt blocks B2/B3/B''/B''', HEGIANT_PICARD.md, runs/scout128, runs/scout897, remap/) and
`/viper/ptmp2/jinma/hegiant_deltaai_bundle`. Nothing over 50 MB is committed (the restarts are 1.2-2.3 GB; the
scout128 run.log is 77 MB and stays on viper).

## Correction first: do NOT use `docs/handover/scripts/he_remap_rst.py`

`docs/handover/scripts/he_remap_rst.py` (md5 25ebdc9e) is byte-identical to `he_remap_rst_w256.py`, the
**unfixed parent**. Our N897 file was made with **`he_giant_1006/remap/he_remap_giant.py`** (md5 09278bf6, last
edited 10-07 17:34:59, used 17:48:36; now `hegiant-files-1007/remap/he_remap_giant.py`). It is generated from the
parent by `mk_remap_giant.py` (string replacements, each asserted unique) and differs in three places that matter:

1. **Exact potential.** The parent fits Phi = a + b/r (point mass) to the old file; wrong for `he_gm_column`. The
   giant script takes old-cell Phi = (E - KE - e_int)/rho from the file and new-cell Phi from `phi_code.phi_cells`
   (the pgen's he_gm_column potential rebuilt from `problem/he_ic_file`, `he_gm`, `he_nfine` read from the rst's
   embedded input; validated at 1.1e-14 against the file, see `phi_old_file_vs_code` below).
2. **M1 flux by ratio, not conservatively** (default; `--m1-conservative` = the old path). E is remapped
   conservatively (limited linear); F = c E_new x the E-weighted overlap mean of the old ratio F/(cE), with |f| set
   to the E-weighted mean of |f_old|, so free streaming (|f| = 1 above ~70 Rsun) survives. The conservative F remap
   gave |F|/cE 0.93-0.9997 there and ~110 Picard passes per solve.
3. **Stored M1 face fluxes f0x1..3 zeroed** (default; `--keep-faces` = interpolate them as before), which puts the
   code on its cold-start closure-lag path. Interpolated faces were inconsistent with the cells.

Also: caches (pred/eint/ctr/rss) are dropped and dt in the file is set to 0.5 x the source dt (both defaults,
`--keep-caches` / `--keep-dt` to change).

Helpers: `he_remap_giant.py` imports only numpy + stdlib + **`phi_code.py`** (must sit in the same directory; the
hard-coded `sys.path.insert(0, '/viper/...')` at line 487 is harmless where the path does not exist). It reads the
IC column at the `problem/he_ic_file` path embedded in the SOURCE rst, so that path must exist on the machine that
remaps (on Caltech your own rst carries your path). `he_remap_rst_w256.py` is needed only to regenerate the giant
script with `mk_remap_giant.py` (its paths are viper-absolute; not needed to run the remap). The diagnostics
`facecmp.py`, `fsat.py`, `teq.py`, `roundtrip_cmp.py` import `he_remap_rst_w256` (and `phi_code`); `mix_rst.py`
imports `he_remap_giant`. `teq.py` also reads a viper EOS dump (`he_giant_1006/eos/eos_dump_y0.9265_a1.77.txt`,
not copied; diagnostic only). All copied verbatim (md5 equal to the files that ran).

## Item 3: the exact remap command

Source `runs/scout128/rst/hegiant.00022.rst` = the wall-clock-limit restart at the end of scout128 link 2
(job 12108730), t = 9.00420777214e5 s (10.42 d), cycle 47861, dt 18.4467 s, 4 MeshBlocks 445x64x64. Command
(PENDING_CMDS.txt block B2 with NNNNN = 00022; output redirected to `runs/scout897/remap.log`, whose mtime equals
the output rst's, 10-07 17:48; CPU, login node, ~1 min, ~11 GB RSS):

```
cd /viper/ptmp2/jinma/he_giant_1006/remap
nice python3 /viper/ptmp2/jinma/he_giant_1006/remap/he_remap_giant.py \
  /viper/ptmp2/jinma/he_giant_1006/runs/scout128/rst/hegiant.00022.rst \
  /viper/ptmp2/jinma/he_giant_1006/runs/scout897/rst/hegiant.00022.rst \
  --factor 1 --mbx2 64 --mbx3 64 --grid /viper/ptmp2/jinma/he_giant_1006/grid/p2_897.npy --nx1 897 \
  > /viper/ptmp2/jinma/he_giant_1006/runs/scout897/remap.log 2>&1
```

Everything else default: slopes ON (the file Caltech asked about is the "slopes" one; `remap/r897ns22` =
`--no-slopes` variant, not used), M1 ratio remap, faces zeroed, caches dropped, ghosts `wrap`, dt x 0.5.
`--grid` = 18 StretchRPoly parameters of the new radial grid (`remap/grid/p2_897.npy`, made by
`grid/fit_giant2.py`, report `grid/fit2_897.out`); x1min/x1max/ng unchanged. Caltech equivalent:
`python3 remap/he_remap_giant.py <src>/hegiant.NNNNN.rst <dst>/hegiant.NNNNN.rst --factor 1 --mbx2 64 --mbx3 64 --grid remap/grid/p2_897.npy --nx1 897`.

Checks, from `remap/remap_scout897_00022.log` (copy of `runs/scout897/remap.log`). `rad_cons_X` = max over the
128x128 radial columns of |sum_new - sum_old| / max|sum_old| (volume-weighted column totals over active cells):

| check | value |
|---|---|
| `phi_old_file_vs_code` | 1.14e-14 |
| mass `rad_cons_mass` | 3.55e-16 |
| energy U = KE + e_int (gravity added back with the new Phi) `rad_cons_U` | 5.54e-16 |
| radial momentum `rad_cons_m_r` | 4.24e-16 |
| angular momenta r m_theta / r m_phi `rad_cons_r_m_th` / `_ph` | 3.76e-16 / 4.09e-16 |
| radiation E `rad_cons_Erad` | 6.59e-16 |
| total F_r change `rad_Fr_change_rel` (not conserved by design) | **7.24e-3 (0.72 %)** for this 10.42 d file |
| `rad_hyd_bad_after`; hydro slopes dropped in old cells | 0; 33672 |
| `rad_eint_min_ratio`; `rad_rho_min` | 0.50002; 1.0e-20 |
| `rad_F_ratio_max` (|F|/cE after); `rad_F_clipped` | 1.0000000000000004; 960395 cells clipped to exactly cE (round-off) |
| m1 faces | zeroed |

The "0.6 %" in our notes is the 5.17 d remap (rst 00011 -> N897, `remap/remap_00011_fwd_z.log`:
`rad_Fr_change_rel` 6.22e-3); for 00022 it is 0.72 %. The "fallbacks (coarse cells -> injection) ... 1.0000"
lines in the log are the angular part, which is the identity here (--factor 1, same 128x128).

## N897 vs N445 grid keys

Only `<mesh>/nx1` (445 -> 897), `<meshblock>/nx1` (445 -> 897) and the 18 `f_stretch_r_*` keys differ between the
inputs embedded in scout128 `hegiant.00022.rst` and scout897 `hegiant.00022.rst` (diff of the two PAR_DUMPs; the
remap writes them from the npy; npy vs printed keys agree to 4e-13 relative). Both: x1min 2.0871e11 (3 Rsun),
x1max 1.3914e13 (200 Rsun), `use_grid_stretch_r_poly = true`, nx2 = nx3 = 128, meshblock 64x64.

| key | N445 | N897 |
|---|---|---|
| c1..c8 | -1.236894517393e+00, -4.735163641768e+00, 5.026104492976e+01, -3.570419665727e+02, 1.422488454204e+03, -2.717229013236e+03, 2.437978328868e+03, -8.386699985201e+02 | -1.272273799858e+00, -1.977330494356e+01, 2.926701797839e+02, -1.484689926654e+03, 3.739183469225e+03, -5.144911304399e+03, 3.715116397871e+03, -1.113061527953e+03 |
| b1_amp, b1_x, b1_w | -2.021640817960e+00, 4.640065359922e-01, 1.525464300752e-01 | -9.491339281175e-01, 2.407883898117e-01, 7.665554490394e-02 |
| b2_amp, b2_x, b2_w | 2.054617003474e+00, 8.798116659071e-01, 1.046282952739e-01 | -2.998805320338e+00, 9.487473032261e-01, 5.791462410538e-02 |
| p_amp, p_xa, p_xb, p_w | -9.988878672034e-01, 3.705211949536e-01, 8.168730212611e-01, 2.057947745856e-02 | -9.989996883444e-01, 1.952251964931e-01, 9.869232597656e-01, 2.647092631281e-03 |

N897 = fine band 60-95 Rsun (dr ~0.04-0.08 Rsun, dr min 0.030 Rsun, max neighbour ratio 1.20; `grid/fit2_897.out`).
The c1 line in the N897 input keeps the stale trailing comment "grid/fit_giant.py 445"; the values are N897's.

## Item 2: N897 input

`hegiant_n897_rst_embedded.athinput` (md5 ad8c17b9) = the PAR_DUMP embedded in `runs/scout897/rst/hegiant.00022.rst`,
from the file start up to (not including) `<par_end>`, the three `problem/` paths replaced by `@BUNDLE@/ic/...`,
`@BUNDLE@/tables/...` (lines 325-327). Identical to the DeltaAI bundle template `tmpl/hegiant_n897_rst_embedded.athinput.in`.
It contains no restart keys (XKEYS are command-line only). It is a reference: on restart the code reads the
input embedded in the rst you write with the remap.

## Item 4: restart keys of the first N897 link

`n897/run.cfg_scout897` (as left on viper, B''' block; chain 12119930-39 PENDING at writing, user decides):

```
BIN=/viper/ptmp2/jinma/builds/bin/athena_he_gpu72_f3a66907_hegopn      (md5 fc7b33a1, branch hegiant-opn-1007)
TLIM=2.592e6
XKEYS="time/cfl_number=0.3 time/restart_refill_ghosts=true rad_m1/implicit_opac_newton_slope_max=3"
PK=""     (viper: the rst's own paths; elsewhere = restart_keys.txt below)
```

- `restart_keys.txt` (PK, bundle SETUP.sh): `problem/he_ic_file=<dir>/ic_giant_own.txt
  problem/he_opac_table=<dir>/rosseland_tops_hegiant_blend.txt problem/he_planck_table=<dir>/planck_tops_hegiant_blend.txt`
  on every link (the IC column and both tables are re-read on every restart).
- Output last_time: `dual_hegiant.sh` runs `n897/rst_info.py <rst>`, which prints `outputN/last_time =
  floor(t/dt_N)*dt_N` for every output with dt > 0, and puts them on the command line. For t = 900420.777 s:
  `output1/last_time=900400 output2/last_time=885600 output3/last_time=885600 output4/last_time=900000
  output5/last_time=864000 output6/last_time=885600` (hst 25 s; bins 0.25 d; log 1e4 s; rst 0.5 d). Full command of
  the exact-config smoke (`n897/dual.12120235.out` line 9):
  `srun -n 2 <BIN> -r <rst>/hegiant.00022.rst <the 6 last_time keys> time/tlim=2.592e6 time/cfl_number=0.3
  time/restart_refill_ghosts=true rad_m1/implicit_opac_newton_slope_max=3 time/nlim=47921 time/ndiag=1 -d <dir> -t <hh:mm:ss>`.
- `time/restart_refill_ghosts=true`: first N897 link only (the restart is a remapped file). Once link 1 has
  started, drop it (`XKEYS="time/cfl_number=0.3 rad_m1/implicit_opac_newton_slope_max=3"`) and re-smoke.
- `rad_m1/implicit_opac_newton_slope_max=3`: read only when named (absent = old path). On the 10.42 d N897 file
  (25 cycles, job 12119980, HEGIANT_PICARD.md): base Picard 76.7/200, 5 NON-CONV, 5.21 s/cycle; slope_max 3:
  19.5/42, 0 NON-CONV, 2.96 s/cycle, max rel diff to base 1.9e-10 (L_top). Exact-config smoke 12120235 (60 cycles,
  1 node 2 MI300A, 2 blocks/GPU): rc 0, FATAL 0, Picard solves 60 mean 18.5 max 42 NON-CONV 0, 960 NEWTON-FALLBACK
  lines, he_ic_balance 0.00515103 at r = 2.15997e+11 (the IC balance on the N897 grid), dt 9.2234 s at the start
  (file dt x 0.5) -> 1.858862e+01 at cycle 47921 (t 9.015257e+05); `n897/` has its run.log, both hst, job out.
- Changed between scout128 and N897 beyond the grid: (a) binary f7291255 (`_hegiant`, md5 ae8a3d11) -> f3a66907
  (`_hegopn`, md5 fc7b33a1 = hegiant-1006 + the opacity-cliff keys; keys off bitwise = f7291255 on the GPU He giant
  rst 00011, 30 cycles); (b) XKEYS + `restart_refill_ghosts=true` (link 1) and `implicit_opac_newton_slope_max=3`;
  (c) the remap side effects: caches dropped, M1 faces zeroed, file dt halved. The embedded input is otherwise
  identical (diff = grid keys only), so the MLT scaffold ramp (5 -> 10 d), base heat, seed, outputs, cfl and all
  rad_m1 keys are unchanged; at 10.42 d the ramp is already over. An earlier option (B'', `run.cfg_scout897.bak_opnoff`)
  used `rad_m1/implicit_opac_newton=false` with f7291255; it was superseded by B''' (slope_max 3) before any link ran.

## Item 5: scout128 key changes before 10.42 d

None. `runs/scout128/run.cfg` was written 10-07 03:20 (before link 1) and never edited (no backups exist); the
input `hegiant_scout128.athinput` (md5 00e7dfbc) was last written 03:31, also before link 1. Both links print the
same config (`scout128_links/dual.*.out`):

| link | job | start -> end | command |
|---|---|---|---|
| 1 | 12108729 | 10-07 04:22 -> 10:02, rc 0 | fresh: `srun -n 4 athena_he_gpu72_f7291255_hegiant -i hegiant_scout128.athinput time/tlim=2.592e6 time/cfl_number=0.3 -d ... -t 05:39:58` -> rst 00011 (t 4.467e5 s, 5.17 d) |
| 2 | 12108730 | 10-07 11:43 -> 17:23, rc 0 | `-r rst/hegiant.00011.rst output1/last_time=447100 output2/last_time=432000 output3/last_time=432000 output4/last_time=440000 output5/last_time=432000 output6/last_time=432000 time/tlim=2.592e6 time/cfl_number=0.3 -d ... -t 05:39:59` -> rst 00022 (t 9.004208e5 s, wall-clock dump) |
| 3-6 | 12108731-34 | CANCELLED 10-07 21:13 (never started) | - |

XKEYS = `time/cfl_number=0.3` and PK empty in both; no refill key on link 2 (plain restart). Layout 2 apu nodes x 2
MI300A, 4 ranks, 1 block per GPU. Note for Caltech: 10.42 d is where link 2's wall clock ran out, not an output
cadence point; stop your N445 run with `time/tlim=9.004208e5` (or remap your nearest restart) to match.

## Item 1: reference smoke 12123795

Table in `NOTE-2026-10-08-viper-hegiant-files-caltech.md` extended (dt at cycle 60 = 2.004713e+01, t 1.202815e+03;
last hydro and user hst rows; NEWTON-FALLBACK 450 lines = 225 events printed twice; Picard solves 60 mean 4.95 max
10 NON-CONV 0; he_ic_balance 0.00369707 at r = 2.15098e+11, face residual 7.99361e-15, last change 8.31531e-10).
Caltech's H200 smoke vs viper: last hst t 1202.8153622134355 vs ...134232 (1.0e-14 rel), mass
1.5435376857928201e32 vs ...7927376e32 (5.3e-14 rel), Picard and he_ic_balance identical in the printed digits.
