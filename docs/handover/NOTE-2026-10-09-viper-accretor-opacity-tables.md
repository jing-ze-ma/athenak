# NOTE 2026-10-09 viper: X 0.70 Z 0.02 GS98 opacity tables for the accretor (answers TASK-2026-10-08-viper-accretor-opacity-tables)

From: viper. To: Caltech. Docs and data only; no code was changed and no run was made. Files are in `docs/handover/accretor-tables-1008/`
(md5 in its `MD5SUMS`).

## Tables
| | file | md5 |
|---|---|---|
| Rosseland (includes electron scattering) | `rosseland_ext2_gs98_x0.70_z0.02.txt` | **ad516c7e** ae541fd3aa34d6048c2c2c0a |
| Planck (absorption only) | `planck_ext2_gs98_x0.70_z0.02.txt` | **99d668f6** 9b16a502f50f1ec555fb0cd1 |

- Format: as the He giant / AG Car tables (he_star_m1 `HsReadOpacityTable`). It is written by the same writer as the AG
  Car ext2 pair, which runs.
- Grid: `217 421 2.6 0.025 -21 0.05`, i.e. log T 2.6-8.0 (dlog T 0.025) and log rho -21..0 (dlog rho 0.05), T slowest,
  log10 kappa.
- Because the table starts at -21, an `*_opac_logd_min` extension key below -21 is a no-op.

## Source and mix
- **Mix:** GS98 scaled solar, X 0.700, Y 0.280, Z 0.020. Number fractions come from `planck_tools/mixspec.py 0.70 0.02`
  (`mixture_x0.70_z0.02.txt`). Species below 1e-7 by number are dropped (Co here); realised Z is 0.01999.
- **High T: LANL TOPS/ATOMIC** via `fetch_tops.py`.
  - T 0.0005-10 keV, in two queries joined at 1 keV, where the rows are identical (`merge_tops_T.py`).
  - Densities: 71 log points over 1e-14..1, plus a separate low-density query of 36 points over 1e-21..1e-14.
  - TOPS data therefore reach log T 8.065, so the whole grid up to 8.0 is data rather than an edge fill. This is the
    only change to the AG Car procedure.
- **Low T: Ferguson et al. 2005, GS98.**
  - Rosseland `g98.7.02.tron` from `f05.gs98.tar.gz`; Planck `g98.pl.7.02.tpon` from `f05.g98.pl.tar.gz` (Wichita
    State).
  - X 0.7 exists directly, so there is no X interpolation. Both files match their tarball members.
- **Procedure: exactly the AG Car ext2 one.** The script is `/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/scripts/build_ext.py fergR`,
  copied as `scripts/build_ext.py`. Only paths, file names, X and headers were changed.
  - Ferguson is used inside its box (log T 2.7-4.5, log R -8..1), with convert_ferguson's bilinear interpolation.
  - The blend is linear in log T of log kappa over **log T 4.0-4.2**.
  - At log T >= 4.2 and log rho >= -14 the values are bitwise the TOPS-only table: `tables_tops_only/*` from
    `convert_tops.py`, 42993/42993 nodes equal in both tables.
- **Physical low-density extension** below each source's own floor rho_e(T); it is not an edge fill.
  - TOPS floor: TOPS returns clamped values below a T-dependent density. For X 0.70 the floor in log rho is -15.8 at
    5800 K, -14 at 1e5 K, -12.6 at 1 MK, -11.2 at 8 MK and -9.4 at 116 MK. Ferguson's floor is log R -8.
  - Rosseland below the floor: kR = Saha e-scattering (`saha_es.py`, LTE: H, He I/II and first ionisation of the metals)
    plus the edge absorption scaled by (rho/rho_e)^s.
  - Planck below the floor: kP = kP_e (rho/rho_e)^s.
  - The slope s is the local log slope over 0.6 dex, clamped to [0, 1] and boxcar-averaged over ±0.1 dex in log T.
  - This extension is a model, not data, and the table header says so.
- **Edge-filled (not data):** log T < 2.7 and, at log T < 4.0, log R > 1 (dense cool gas; e.g. rho > 1e-5 at 1e4 K).

## Checks (`checks/table_check_x070.txt`, `checks/blend_jump_locations.txt`)
- All values are finite.
- **Electron-scattering limit (0.2(1+X) = 0.34):**
  - TOPS kR = 0.3306 at log T 7.0 (log rho -10 and -18), 0.3331 at (6.83, 5e-16) and 0.3376 at (6.0, 1e-20).
  - The Saha es model at full ionisation gives 0.3347. TOPS is about 1 % lower (its own plasma/Compton treatment); the
    repo OPLIB X 0.7 table gives 0.3307.
  - At higher T Klein-Nishina lowers it: 0.314 at log T 7.5 and 0.274 at 8.0, both as in OPLIB (0.3143 / 0.2740).
  - The old dev pair edge-fills 0.3296 above log T 7.065.
- **Fe bump**, max kR at log T 4.9-5.6:

  | log rho | new Z 0.02 (log T) | dev TOPS Z 0.008 | repo OPLIB X 0.7 Z 0.014 |
  |---|---|---|---|
  | -12 | 0.551 (5.150) | 0.481 | 0.658 |
  | -10 | 0.761 (5.175) | 0.605 | 0.718 |
  | -9 | 0.967 (5.175) | 0.719 | 0.907 |
  | -8 | 1.522 (5.225) | 1.020 | 1.358 |
  | -7.2 | 4.249 (4.925) | 3.127 | 3.774 |
  | -6 | 53.5 (4.925) | 40.1 | 46.9 |

  - The bump sits at the same log T in all three tables.
  - The Z 0.02 table is 26-49 % above the Z 0.008 dev pair and 6-14 % above OPLIB Z 0.014 at log rho >= -10.
  - At log rho -12 it is 16 % below OPLIB (TOPS vs OPLIB at low density).
- **Continuity at the blend** (max |dlog k| per 0.025 step in log T, log rho <= -8):
  - Rosseland: 0.18 (4.0->4.025), 0.10 (4.1), 0.09 (4.175->4.2), 0.10 (4.2->4.225).
  - Planck: 0.18, 0.15, 0.26 and 0.24 at the same steps.
  - Context: in the H/He ionisation ramp over log T 3.5-4.6 the table's own steps reach 0.59 (R) and 0.60 (P).
  - Larger steps in the window (R 0.65 at 4.0, 0.44 at 4.1) occur only at log R > 4.7 (log rho > -1), which is the
    Ferguson edge-fill region.
  - Along density the largest step is 0.11 (R) and 0.17 (P) per 0.05 dex.
- **Ratio kP/kR:**
  - Log T 3.5-5.0, log rho -18..0: kP >= kR everywhere. log kP/kR has a median of +2.1 to +2.4 and a minimum of
    +0.01 to +0.39.
  - Log T 5-6: kP < kR in 22 % of nodes. Log T 6-7: 81 %. Log T 7-8: 100 %.
  - Examples: (5.2, 1e-10) 62; (6.0, 1e-6) 1.96; (6.83, 5e-16) 3.9e-10. In the last, the hot ambient is pure
    scattering, so the absorption mask is right.
  - Consequence: eps = min(kP/kR, 1) = 1 below log T 5 and < 1 in the hot gas, as it should be.
- **Ext2 vs the TOPS-only table below log T 4.2** (log rho -14..-8; median / max |d| in dex):
  - Log T 4.0-4.2: R -0.008 / 0.11; P -0.004 / 0.80.
  - Log T 3.76-4.0: R -0.008 / 0.41; P +0.03 / 1.04.
  - TOPS's coarse T nodes interpolate across the H ramp; this is the same finding as for AG Car.

## Differences vs the old dev pair (TOPS-only X 0.7 Z 0.008, bsg_1001/opac, md5 0cf58c08 / 8255eae0)
- **Z 0.008 -> 0.02:**
  - In the log rho -14..0 overlap, kR rises by median +0.07 dex at log T 5.0-5.5 (up to +0.26) and is otherwise
    +0.0-0.03 dex.
  - kP rises by median +0.36 to +0.40 dex at log T >= 5 (metal absorption, about x Z).
- **Low T:** below log T 4.2 it is Ferguson (molecules, H-, grains below ~2000 K), blended. The dev pair is TOPS down to
  3.764 and edge-filled below.
- **Low density:**
  - The dev pair stops at log rho -14 (the code's first-interval slope extension takes over below).
  - The new tables use real TOPS data down to its own floor (e.g. -15.8 at 5800 K), then the Saha-es + absorption
    extension down to -21.
  - Example: kR at (6.0, 1e-20) = 0.338.
- **High T:** TOPS data up to log T 8.065, where the dev pair edge-fills above 7.065. kR at log T 8 is 0.274 vs 0.330.

## Flags
- **Extension below the floors is a model.**
  - At log T ~5 and rho 1e-20, kR = 0.40 > kes: the absorption slope fitted at the TOPS edge is ~0.05 (bound-bound),
    so line absorption barely drops over 6 dex.
  - The Planck mean in the thin gas below the floors is uncertain to ~1 dex, as in the AG Car NOTE.
- The TOPS-only bitwise region and both ext2 tables were built on the viper login node with nice. No sbatch, no build,
  and no code edit were involved.
- `scripts/build_ext.py` keeps viper absolute paths (provenance). Inputs are under `/viper/ptmp2/jinma/accretor_tables_1009/`.

## Rebuild
```
python3 scripts/mixspec.py 0.70 0.02                       # line 2 -> mixture_x0.70_z0.02.txt
raw/fetch_all.sh                                           # 4 TOPS queries (aphysics2.lanl.gov)
python3 scripts/merge_tops_T.py raw/tops_gs98_x0.70_z0.02.dat raw/tops_hiT_x0.70.dat raw/tops_gs98_x0.70_z0.02_10keV.dat
python3 scripts/merge_tops_T.py raw/tops_lo_x0.70.dat raw/tops_hiT_lo_x0.70.dat raw/tops_lo_x0.70_10keV.dat
python3 scripts/convert_tops.py raw/tops_gs98_x0.70_z0.02_10keV.dat tables/rosseland_tops_gs98_x0.70_z0.02.txt rosseland "<hdr>"   # and planck
python3 scripts/build_ext.py fergR ext2
python3 scripts/table_check_x070.py ext2 tables
```
