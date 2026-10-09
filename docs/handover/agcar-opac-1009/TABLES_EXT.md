# TABLES_EXT.md (10-09): extended opacity tables for the AG Car runs (X 0.36, Z 0.02, GS98)

Everything here ran on the login node with nice, using at most 4 cores. There was no sbatch, no build, no git and no code edit.
Paths are relative to `/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/` unless absolute. No restart was used: the 1-D
runs start fresh from the general-EOS ICs `geos/icA`, `geos/icB`, so there is no restart rotation to state.

## Deliverables (md5 in `MD5SUMS`)
| item | file | md5 |
|---|---|---|
| Rosseland (incl. electron scattering) | `rosseland_ext_gs98_x0.36_z0.02.txt` | 9720fb29 |
| Planck (absorption only) | `planck_ext_gs98_x0.36_z0.02.txt` | a88f5373 |
| **ext2 (recommended) Rosseland, Ferguson pair** | `rosseland_ext2_gs98_x0.36_z0.02.txt` | b28e97c5 |
| **ext2 Planck** | `planck_ext2_gs98_x0.36_z0.02.txt` | dc482452 |
| builder / Saha e-scattering / checks | `scripts/build_ext.py`, `scripts/saha_es.py`, `scripts/table_check_ext.py`, `scripts/check_overlap.py`, `scripts/column_compare.py`, `scripts/ana_run1d.py`, `scripts/fetch_aesopus.py` | MD5SUMS |
| raw data | `raw/tops_lowrho/tops_lo.dat` (TOPS, rho 1e-21..1e-14, f29ab5ef), `raw/aes21_gas/03-GS98.zip` (AESOPUS 2.1, 975317be), `raw/aesopus_form/x036_partial.dat` (web X 0.36, 6 rows) | |
| check outputs | `checks/table_check_ext.txt`, `checks/overlap.txt`, `checks/kappa_columns.{txt,png}`, `checks/run1d_{A,B}.{txt,png}` | |
| 1-D runs | `run1d/{A,B}_{old,new}` (inputs `run1d/in*.athinput`, runner `run1d/run.sh`) | |

The 1-D runs read an earlier build whose data values are identical (`cmp` of all non-comment lines). That build
had md5s 2dcae9ec / d507dd50; only two header lines were corrected afterwards.

**Format and grid.** The format is unchanged: comment lines, then `# nT nD lTmin dlT lDmin dlD`, then T-slowest
log10 kappa. The grid is `217 421 2.6 0.025 -21 0.05`, i.e. log T 2.6-8.0 (unchanged) and log rho **-21..0**
(was -14..0). The density nodes above -14 sit at the old positions.

**Reader.** `HsReadOpacityTable` (he_star_m1.cpp l. 218) takes any uniform grid. `RosselandTable` bisects. No
code change is needed. With the table starting at -21, `he_opac_logd_min = -21` becomes a no-op: the A_new log
has no "extended from" line, while the old-table run prints two.

## 1. Low-T source (molecules, no grains)
- **Rosseland: AESOPUS 2.1, gas only** ("gasbroad": molecules, no grains). Source: the precomputed GS98 set
  at https://stev.oapd.inaf.it/aesopus_2.1/tables/SOLAR/GAS/03-GS98.zip. I used Zref 0.02 at X = 0.35 and
  0.5, linear in X to 0.36 (w 0.0667). Coverage: log T 2.0-4.5, log R -8..+6 (R = rho/T6^3; per T row
  log rho = log R + 3 log T - 18).
  - Access works. The on-demand web form (scripted, `scripts/fetch_aesopus.py`) computes the exact X 0.36,
    but only about 1 T row per 4-5 min. Six rows (log T 4.40-4.50) finished and differ from the X
    interpolation by at most **0.0042 dex**. The interpolated set is therefore used.
  - The mixture is the same as TOPS: GS98 scaled solar, Fe/O by number 0.04677 in both.
- **Planck: Ferguson et al. 2005** (`g98.pl.35.02` and `g98.pl.5.02`, linear in X to 0.36), log T 2.7-4.5,
  log R -8..1. AESOPUS gives no Planck means.
  - **Grain choice:** the Rosseland mean is grain-free. Ferguson's Planck tables include grains, and there is no
    grain-free Ferguson set. Grains condense only below ~1500-2000 K (log T < 3.3), well below tfloor 3000 K
    (A) and 5000 K (B), so they never enter a run.
- **Cross-check at the TOPS temperature nodes** (no T interpolation), rho <= 1e-8, log R >= -8. AESOPUS-Ferguson
  agree to a median of 0.02-0.06 dex and a max of 0.13.

## 2. Low density
- **TOPS accepts rho < 1e-14 but its numbers are clamped.** Querying to 1e-21 returns identical values below a
  T-dependent floor:

  | log T | 3.76 | 4.16 | 4.54 | 4.84 | 5.06 | 5.46 | 6.06 | 6.54 | 6.91 |
  |---|---|---|---|---|---|---|---|---|---|
  | floor log rho | -15.6 | -15.2 | -14.6 | -14.2 | -13.8 | -13.2 | -12.4 | -11.6 | -11.2 |

  The real data from the floor up to -14 are now used for log T < ~5. Above ~1e5 K the old table already
  held clamped values between -14 and the floor. They are electron-scattering dominated (kR 0.27). I kept them
  bitwise.
- AESOPUS and Ferguson stop at log R -8, i.e. rho 2.7e-16 at 3000 K, 1.3e-15 at 5000 K and 7e-15 at 9000 K.
  A's atmosphere (1e-16 at the IC, the 1e-20 floor when evolved) lies below every source.
- **Physical extension below each source's floor rho_e(T)**, per source T row:
  - **Rosseland:** kR = kes(T, rho) + (kR_e - kes(T, rho_e)) (rho/rho_e)^s.
    - kes is LTE Saha electron scattering (`saha_es.py`). It covers H, He I/II and the first ionisation of
      the 19 TOPS metals, which supply n_e where H is neutral.
    - Check against TOPS: at 5800 K and 1e-16 Saha gives 0.142 vs TOPS 0.150; at full ionisation 0.2677.
  - **Planck:** kP = kP_e (rho/rho_e)^s.
  - **Slope s:** the local log-log slope (of the absorption part, for Rosseland) over the 0.6 dex above the
    edge, clamped to [0, 1]. The bounds are physical: 0 = bound-bound lines of the dominant ion stage (constant
    per gram), 1 = two-body continuum (ff, bf, H-).
    - s is boxcar-averaged over ±0.1 dex in log T. The raw per-row slopes made 1.2 dex row-to-row steps at
      log rho -21; after averaging the maximum is 0.57 dex at log T 3.95-4.0, rho <= 1e-8 (`checks/table_check_ext.txt`).
    - Typical TOPS slopes are sR 0-0.34 and sP 0-0.49.
  - This is a model, not data. It is labelled as such in the header.
- es limit (`checks/table_check_ext.txt`): kR(log T 7, rho 1e-10) = 0.265 (old 0.265); (7, 1e-20) 0.265;
  (6, 1e-20) 0.2706; (4.0, 1e-20) 0.206; (3.8, 1e-20) 0.147 (H ionised, He neutral: 0.4 X = 0.144).

## 3. Merge
- **Blend window log T 4.0-4.2, linear in log T of log kappa.** This is the repo recipe of `merge_rosseland.py` and
  `convert_ferguson.py`. I tried and rejected the brief's 3.8-4.0. At TOPS's own nodes with rho <= 1e-8 (no T interpolation):

  | TOPS node log T | 3.764 | 3.843 | 3.968 | 4.065 | 4.162 | 4.241 |
  |---|---|---|---|---|---|---|
  | Rosseland AESOPUS-TOPS max dex | 0.80 | 0.28 | 0.09 | 0.12 | 0.14 | 0.15 |
  | Planck Ferguson-TOPS max dex | 0.72 | 0.47 | **1.28** | 0.12 | 0.08 | 0.20 |

  - TOPS's two lowest nodes deviate by up to 0.8 dex. Its Planck mean at 3.968 deviates by up to 1.3 dex.
    The coarse TOPS T spacing (0.08-0.12 dex) also interpolates across the H-ionisation ramp, where AESOPUS
    has 0.01-0.02 dex steps.
  - From 4.0 to 4.2 both sources agree to 0.15 dex (R) and 0.2 dex (P).
  - AESOPUS sits a systematic median 0.06-0.11 dex below TOPS over 4.0-4.5 (`checks/overlap.txt`).
- **Bitwise:** at log T >= 4.2 and log rho >= -14, all 42993 nodes of both tables carry the old file's strings. My
  recomputation of that region differed from the old file by at most 5e-6 dex, which is the print rounding.
  - The Fe bump is unchanged: max kR 0.645 / 0.831 / 1.327 at log rho -10 / -9 / -8, log T 5.175-5.225, as in
    tables/table_check.txt.
  - In a lookup the density node positions differ by round-off (-21 + 0.05 j vs -14 + 0.05 j, max 1.8e-15). Along
    the IC columns below 0.9 R_ph, old and new kappa differ by <= 3.0e-15 relative (`column_compare`). The
    old production run with he_opac_logd_min -21 had the same round-off shift.
- **Continuity:**
  - At the blend edges, max node-to-node jumps for rho <= 1e-8 are 0.17 / 0.06 (R) and 0.16 / 0.24 (P) dex.
  - The largest Rosseland jump overall is 0.71 dex, at log R > 6, T < 4000 K. That region is above AESOPUS's
    range and edge-filled; it is physically irrelevant.
  - Inside A's atmosphere box (log T 3.45-4.0, log rho -17..-15) the max step is 0.50 (R) / 0.52 (P) dex per
    0.025 in log T: the H-ionisation ramp.

## 4. Along the AG Car columns (`checks/kappa_columns.{txt,png}`)
"Old" is what the code actually used: the old table plus the in-code he_opac_logd_min -21 extension (first-interval slope).
| column, band | T, rho | kR old -> new | kP old -> new |
|---|---|---|---|
| IC A, r < 0.9 R_ph | interior | equal (<= 3e-15) | equal |
| IC A, 0.9-1.0 | 9.1-48 kK | max -0.08 dex | -0.32 dex max |
| IC A, 1.1-3.2 R_ph | 9051 K, 1e-16 | 0.168 -> 0.194 | **0.228 -> 13.5** (+1.77 dex) |
| colA (evolved, t 5.2e7 s), 1.5-3 R_ph | 5000-5340 K, 1e-20 | 0.097 -> 0.144 (Saha es) | **170 -> 1.4-1.9** (-2 dex) |
| IC B, 1.1-3.2 R_ph | 20.1 kK, 1e-16 | 0.250 -> 0.248 | 1730 -> 1490 |
| colB (t 4.1e5 s), 1.5-2 R_ph | 10-12 kK, 1e-16 | 0.19-0.22 -> 0.21 | 1.4-42 -> 42-90 |
| colB, 2-3 R_ph | 4.0-84 kK, 1e-18..1e-16 | down to 0.0022 at the coolest cells | -1.7..+2.3 dex |

- The old atmosphere opacities were the TOPS 5800 K or rho 1e-14 edge values with the code's slope.
- The new kP at A's floor rests on the extension model. At 9000 K Ferguson vs TOPS still disagree by up to 1.3
  dex. **The atmosphere Planck mean is uncertain to ~1 dex**, the Rosseland mean to ~0.1 dex (es dominated).

## 5. 1-D A and B, old vs new table (`run1d/`, `checks/run1d_*.txt/png`)
**Setup**
- Binary `/viper/ptmp2/jinma/fix_1009/bin/athena_cpu_c1` (md5 e67f042e), 1 rank.
- Inputs: the production inputs, A = viper `geos/shake/agcar_shakeA_ge.athinput`, B = Raven
  `files/agcar_shakeB_ge.athinput` with viper paths. Both have implicit_flux blend / tau_f.
  - Added: the four fixbundle swap lines (vet_source_noesrc, implicit_face_opac_n, implicit_pos_floor_solve,
    implicit_bcg_fallback best) and an output-only `eos_table_dump`.
  - Command line: `mesh,meshblock nx2 = nx3 = 8`, `problem/he_seed=0` (identical columns), tlim 1e5.
  - The two arms differ only in the two table paths.
- Smokes (CPU, nlim 3): `run1d/smoke2A_new`, `run1d/smoke2B_new`, rc 0, no FATAL.

**Results at t = 1e5 s**
| | A (80 cycles) | B (435 cycles) |
|---|---|---|
| interior r < 0.9 R_ph max \|drho/rho\|, \|dT/T\|, \|dE/E\| | 4.0e-6, 1.1e-5, 4.0e-5 | 9.5e-7, 4.4e-6, 9.5e-7 |
| 0.9-1.0 R_ph T new/old (median, range) / E | 0.999 (0.973-1.038) / 0.998 (0.898-1.159) | 1.000 / 1.000 |
| 1.0-1.1 R_ph T / E | 1.000 (1.000-1.053) / 1.001 (0.999-1.233) | 1.000 / 1.000 |
| 1.1-3.0 R_ph T / E | 1.000 / 1.000-1.001 (T 3740-6600 K) | 1.000 / 1.000-1.001 (T 8300-14800 K) |
| L/L* at 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 2.9 R_ph old | 0.387 0.286 1.000 0.519 0.999 0.999 1.000 | 0.941 0.655 1.000 0.572 0.998 0.998 0.999 |
| same, new | 0.387 0.286 1.001 0.518 0.999 0.999 1.001 | 0.941 0.655 1.000 0.567 0.998 0.998 0.999 |
| R_ph(tau_R 2/3, own table) old / new | 0.9989 / 1.0000 R_ph | 1.0124 / 1.0124 |
| dt mean (min) old / new | 1252.6 / 1252.6 s | 231.0 (201.7) / 230.6 (201.5) s |
| Picard mean (max) old / new | 10.61 (15) / 10.38 (16) | 5.03 (24) / 5.04 (11) |
| NEWTON-FALLBACK lines old / new; FATAL; NON-CONV | 424 / 443; 0; 0 | 3345 / 3335; 0; 0 |

- In 1-D the atmosphere T barely moves. It sits at radiative equilibrium with the dilute E (A 3740-8700 K), and
  the gas-radiation coupling is fast at either kP. The A photospheric layer (0.9-1.1 R_ph, the blend window)
  moves T by -2.7/+5 % and E by -10/+23 %. B is unaffected.
- The L(1.1 R_ph) dip (0.52 / 0.57) is the same in both arms, so it is not caused by the table.
- Cost: dt and Picard counts are equal within noise. The CPU walls (A 3:47 / 3:36, B 12:52 / 12:45) come from
  4 concurrent login-node runs and are **not** a timing result. No GPU timing was done.
- The 1-D test cannot show the 3-D effect of a ~2 dex kP change at the floor (thermal relaxation time of floor
  gas, sponge, FOFC). Check that in the swap smoke and in the first 3-D output.

**General EOS range:** eos_logt_min 3.4 (A) / 3.6 (B) and eos_logd_min -21 cover tfloor 3000 / 5000 K and
dfloor 1e-20. The new tables cover log T 2.6-8.0 and log rho -21..0, which encloses the EOS box. The EOS
dump of the A run (`run1d/A_new/eos_dump.txt`, `-21 0.02 3.4 0.004`) confirms the A range. The 1-D A
atmosphere minimum is 3741 K, above tfloor.

## 6. ext2 = FERGUSON PAIR (user follow-up): low-T Rosseland AND Planck from Ferguson et al. 2005 GS98
**Files** (md5 in `MD5SUMS`)
- `rosseland_ext2_gs98_x0.36_z0.02.txt` (md5 **b28e97c5**)
- `planck_ext2_gs98_x0.36_z0.02.txt` (md5 **dc482452**)
- Build command: `scripts/build_ext.py fergR ext2`. Grid and format are as in ext; header line "ext2 = FERGUSON PAIR".

**Source**
- Rosseland: Wichita State `f05.gs98.tar.gz` (md5 96a1b1db): `g98.35.02.tron` and `g98.5.02.tron`.
- Planck: `f05.g98.pl.tar.gz` (md5 69d8b6d3): `g98.pl.35.02.tpon` and `g98.pl.5.02.tpon`.
- The files in bench/m1_opac/ferguson05 are byte-identical to the tarball members. All four headers read
  "Grevesse & Sauval 1998 ... X= 0.35/0.5, Z= 0.02".
- **Mixture = GS98, the same as TOPS and AESOPUS, so there is no mixture mismatch at the seam.** The earlier ext
  Planck already came from this GS98 set and was not changed in source.

**Interpolation**
- Inside the Ferguson box (log T 2.70-4.50, log R -8..1) both means use convert_ferguson.py's own interpolation
  (`bilin`, bilinear in log T and log R). They match it to 0 dex.
- X is interpolated linearly to 0.36 (w 0.0667). This is our step: convert_ferguson only interpolates in Z.
- The ext Planck had used per-row log-rho interpolation, which differs from `bilin` by at most 0.024 dex. That is
  the only change to the Planck values in ext2.

**Same as ext:** the blend (log T 4.0-4.2), the low-rho extension and the bitwise interior. In both tables all
42993 nodes at log T >= 4.2, log rho >= -14 equal the old table (`checks/table_check_ext2.txt`). Grains are in
both means. They matter only below ~2000 K, under both tfloors.

**Caveat:** Ferguson stops at log R = 1, AESOPUS at 6. For log T < 4.0 and log R > 1 the ext2 values are
edge-filled, i.e. dense cool gas. No AG Car cell is there: on all four columns the max log R at log T < 4.2 is
-4.5 (IC A), -4.6 (colA) and -6.3 (colB).

**(2) Source uncertainty: Rosseland Ferguson - AESOPUS** (`checks/ferg_vs_aes.txt`), at Ferguson nodes where both
are data, log R -8..-3:

| log T | 3.4-3.6 | 3.6-3.7 | 3.7-3.8 | 3.8-3.9 | 3.9-4.0 | 4.0-4.1 | 4.1-4.2 | 4.2-4.5 |
|---|---|---|---|---|---|---|---|---|
| median dex | +0.02 | -0.08 | +0.02 | +0.03 | +0.02 | +0.03 | +0.04 | +0.03 |
| max dex | 0.12 | 0.18 | 0.08 | 0.07 | 0.03 | 0.05 | 0.07 | 0.19 |

Over all log R from -8 to 1 the maximum is 0.51 dex, at log T 3.4-3.6, high R. Along the columns kR(ext2)/kR(ext)
is within **-0.02..+0.05 dex**.

**kP/kR along the columns** (median, with the range in brackets; old = old table + code extension):
| column, band | old | ext | ext2 |
|---|---|---|---|
| IC A 1.0-1.1 R_ph | 34 (1.35-2230) | 207 (70-1420) | 197 (69-1340) |
| IC A 1.1-3.2 R_ph (9051 K, 1e-16) | 1.35 | 69.8 | 69.1 |
| colA 1.1-1.5 (4.8-10 kK, 1e-20) | 631 (**0.010**-1760) | 7.0 (1.66-99) | 6.9 (1.68-99) |
| colA 1.5-3.2 (5.0-5.3 kK, 1e-20) | 1760 | 12.8 (9.5-12.9) | 12.8 (9.5-12.8) |
| IC B 1.1-3.2 (20.1 kK, 1e-16) | 6910 | 6000 | 6000 |
| colB 1.1-1.5 | 397 (196-1570) | 621 (432-1590) | 613 (425-1570) |
| colB 1.5-3.2 | 16 (**0.27**-8700) | 209 (9.4-7320) | 206 (9.5-7320) |

- With ext and ext2, kP/kR >= 1.66 in every column cell. So **eps = min(kP/kR, 1) = 1 everywhere** and the VET source
  weight does not depend on the choice.
- The old table gave eps < 1 in some cells (colA min 0.010, colB min 0.27).
- The gas coupling (kP) is the same in ext and ext2 to 0.024 dex.

**(3) 1-D check, ext2 vs old** (same setup as section 5; smokes `run1d/smoke4{A,B}_ext2` rc 0, no FATAL;
`checks/run1d_{A,B}_ext2.txt`):
| | A | B |
|---|---|---|
| interior max \|drho\|, \|dT\|, \|dE\| | 1.7e-6, 6.6e-6, 1.2e-5 | 9.5e-7, 4.4e-6, 9.5e-7 |
| 0.9-1.0 R_ph T / E new/old (range) | 0.991-1.051 / 0.965-1.218 | 1.000 / 1.000 |
| 1.0-1.1 R_ph T / E | 1.000-1.064 / 1.001-1.283 | 1.000 / 1.000 |
| 1.1-3 R_ph T / E | 1.000 / 1.000-1.001 | 1.000 / 1.000 |
| L/L* at 0.5 ... 2.9 R_ph | old and ext2 equal to 0.001 | equal to 0.001 (L(1.1) 0.572 -> 0.568) |
| dt mean, Picard mean (max), FATAL / NON-CONV | 1252.6 s, 10.23 (15), 0 / 0 | 230.5 s, 5.05 (11), 0 / 0 |

ext2 and ext differ by about 1 % in the A photospheric layer and are identical in B.

**Recommendation: the Ferguson pair (ext2).**
- kR and kP come from one code: same line data, EOS, molecules and grain treatment, same GS98 mixture as TOPS.
  That is the consistency kP/kR needs.
- Its Rosseland mean is within 0.18 dex of AESOPUS (median <= 0.08) at the densities AG Car reaches. Along the
  columns ext2 and ext differ by <= 0.05 dex in kR. eps is 1 in both, and the 1-D results agree.
- The only drawback, the edge-fill at log R > 1 for T < 10^4 K, is far from every AG Car state.
- The dominant uncertainty is the same for both choices: the thin-atmosphere Planck extension (Flags).

## Proposed swap (NOT applied; production inputs and tables untouched); recommended = ext2
In `<problem>` of `geos/shake/agcar_shakeA_ge.athinput` (viper) and of the Raven `files/agcar_shake?_ge.athinput`
(copy the two files to Raven `files/` first):
```
he_opac_table   = /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/rosseland_ext2_gs98_x0.36_z0.02.txt   # md5 b28e97c5 (Ferguson pair, tables_ext 10-09)
he_planck_table = /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/planck_ext2_gs98_x0.36_z0.02.txt      # md5 dc482452
he_opac_logd_min = -21.0         # now a no-op (the table starts at -21); keep or delete
```
(Alternative: ext, AESOPUS-R + Ferguson-P, md5 9720fb29 / a88f5373.) Before any release, the smoke rule applies:
an apudev/gpudev smoke with the exact script, input and keys.

## Flags / open
- **The Planck mean in the thin atmosphere is model-dependent by ~1-2 dex.** The extension is below every
  available table, and Ferguson and TOPS differ by up to 1.3 dex at 9000 K. A physical fix would be an LTE
  line-list Planck mean at rho 1e-20..1e-15 and T 3-10 kK, which is outside this task.
- The Rosseland mean there is es + small absorption and is robust (Saha vs TOPS within 6 %).
- The AESOPUS web job for exact X 0.36 (output409103337749) was still computing when I stopped. Its output is
  deleted 2 h after completion; the X-interpolation error is 0.004 dex anyway.
- `checks/kappa_columns.txt` shows the in-code first-interval extension of the OLD table made kP at A's IC
  atmosphere 0.228. That is 1.1 dex below the old table's own -14 edge value (2.91).
