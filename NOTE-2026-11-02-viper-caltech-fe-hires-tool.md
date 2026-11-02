# NOTE 2026-11-02 viper -> Caltech (USER): fe_hires tool pushed + Seidel+2026 observable checklist

Reply to TASK-2026-11-02-caltech-viper-fe-hires-tool.

## 1. The tool: `docs/handover/scripts/fe_hires/` (this branch)
- **What it is.** A portable copy of the viper pipeline. Start with `README.md`. The original
  write-up is `FE_HIRES_w121_results.md`.
- **Physics unchanged.** The WASP-121b defaults reproduce the 10-07 result.
- **No /viper paths.** Everything comes from CLI arguments plus a planet JSON:
  ```
  python fe3d.py --planet planets/P.json --eos <run eos_table> --out OUT <dump.bin> [key=val]
  python ana.py  --planet planets/P.json --out OUT <tag>
  ```
- **Per-planet hook** (`planets/*.json`; keys left out take the WASP-121b values):
  - system: R*, Rp, a/R*, b, P, omega, Kp0 (derived or set), quadratic limb darkening (default uniform = baseline);
  - instrument: LSF FWHM or R;
  - exposures: nph, phase range, RM-mask phase;
  - composition: eos_xh, eos_yhe, comptab table, a_fe, Fe top cut.
- **Inputs in git** (`data/`, each <= 5.2 MB):
  - the Fe I line list (Kurucz gfall, 9620 lines), Fe I / Fe II partition functions;
  - the pRT 3.4 3000 K template (raw; the LSF is applied at run time), H2/He Rayleigh at the continuum nodes;
  - composition tables for **3x** (WASP-121b MHD runs) and **1x** (X 0.7381, Y 0.2485, [M/H] 0: your runs).
  - The 1x table was built with `comptab.cpp` against AthenaK `src/eos/eos_composition.hpp` (unchanged since 8bc1aea0d). The 3x table built the same way is byte-identical to the 10-07 one.
- **Not in git:** the Kurucz raw files and the 7.3 GB pRT lbl Fe file. They are needed only to *rebuild* the cached inputs. Size, md5 and provenance are in the README.
- **Run-time dependencies:** numpy, scipy, numba. No petitRADTRANS, FastChem or dhjcs.py; the EOS inversion, geometry and bin reader are vendored.
- **Ready for your runs:** `planets/wasp121b_1x_caltech.json` + `--eos .../anchor/eos_1x.txt`.
  - Other planets: copy that JSON and set the system keys.
  - Dumps must be version 1.2, i.e. carry `mb_x1f`. This is needed because your decks use the `f_stretch_r_b*` bump terms.

**Reproduction check (viper login node, CPU, 16 processes, 24 min wall, RSS 3.0 GB):**
WASP-121b s10_e15 dump 00166 (10 G MHD, rot 318). The log is `run_w121.log`.

| km/s | 10-07 | portable 11-02 |
|---|---|---|
| Kp-Vsys Gaussian peak at Kp0 | -4.26 | **-4.26** |
| centroid / single Gaussian | -3.78 / -5.07 | -3.78 / -5.07 |
| map max (Kp, v) | 206.2, -3.60 | 206.2, -3.60 |
| morning / evening segment | -3.43 / -6.78 | -3.43 / -6.78 |

The per-exposure spectra are bit-identical (max |diff| = 0).

**Cost per dump:**
- 31 exposures at about 450-670 s each: about 10-12 min per dump on 31 cores.
- 11 dumps x 5 planets: about 280 core-h.

## 2. The observable vs Seidel+2026 Fig. 1

| step | Seidel+2026 | fe_hires now | reproduces? | missing / cost |
|---|---|---|---|---|
| template | Kitzmann+2023 (Mantis) Fe I, 3000 K isothermal, hydrostatic + chemical equilibrium | pRT 3.4 lbl Fe I, 3000 K isothermal, solar Fe all neutral, H2/He continuum | **partly**: same T and species, different line list/weights; sub-km/s offset untested | Mantis 3000 K Fe template -> `data/`, rerun ana only. ~0.5 day; ana ~2 min/dump |
| Kp-vsys + Kp marginalisation | stacked Kp-vsys (nights S/N-weighted, per-night vsys), offset "marginalised over Kp" | Kp-vsys map (equal weights, vsys = 0 exactly); offset = Gaussian at the Kp0 row; map max also given | **partly**: no Kp-marginalised number. s10: Kp0 row -4.26 vs map max -3.60 at Kp 206, so the choice matters at ~0.7 km/s | Kp-marginalised offset (projection over a Kp window, Gaussian fit) in ana.py. ~1-2 h, no fe3d rerun |
| transit phases | real exposures per night; per-night phases not published (ED Table 1 has dates and S/N only) | 31 uniform exposures T1-T4, contact points dropped; `phase_range`, `nph`, `phases=` hooks | **no** for the real coverage (e.g. the WASP-121b 4-UT nights are partial: 2023 ingress, 2018 egress) | per-night phase lists from ESO/Gemini archive headers + a phase list per night. ~0.5-1 day per planet; fe3d rerun per dump |
| RM / Doppler shadow | corrected (StarRotator, quadratic LD), not masked | model spectra have no shadow; optional mask \|phi\| < 0.017 (the `rmcut` set; Seidel+25); limb-darkening hook added, default off | **yes in effect**: their correction aims at what the model already is. Residuals not modelled | per-planet LD (Claret) = 1 h + fe3d rerun. Shadow residuals: not planned |
| instrument LSF | ESPRESSO 1-UT 2.14, 4-UT MR 4.28, MAROON-X 3.53 km/s | one Gaussian per planet (`lsf_fwhm_kms`), applied to spectra and template (4.28 for WASP-121b) | **yes** for single-instrument planets (WASP-178b 2.14, TOI-1518b 3.53, WASP-121b 4.28) | WASP-76b and WASP-189b mix ESPRESSO + MAROON-X: per-instrument ana + weighted stack. ~2 h code; weights only approximate (ED Table 1 gives S/N ranges) |
| Fe top cut | n/a (data). Paper CF: 1-100 mbar, peak 30 mbar | Fe removed at p < 1e-6 bar (`pcut_bar`, key `pcut=`). Our CF peaks at ~1e-4 bar | **assumption**: s10 evening -6.78 with the cut, -3.96 without | run each dump with `pcut=1e-6` and `pcut=0` and report both (2x cost). Real fix (Fe photoionisation at the top) not planned |

**Further differences (not in the TASK list):**
- no exposure smearing;
- Fe/H = 3.2e-5 (solar), not scaled with [M/H] (fine for your 1x runs);
- no Fe condensation (the `tcond=2000` variant moved WASP-121b by <= 0.14 km/s);
- one dump per number unless you average (the 10-07 dump spread was not measured);
- model line widths are 1.5-3x the observed.

**Verdict:**
- fe_hires gives the same kind of number as Fig. 1: the Fe I CCF Kp-vsys offset at the planet Kp.
- It is not yet the paper's exact estimator. The Kp marginalisation and the real phase coverage are the two gaps that can move the result by ~1 km/s.

## 3. The paper's data table
- **Supplementary Table 1 is not available to us.** It is not in the ESO PDF.
- **No per-night phase ranges exist** in the paper.
- **What we have** is the digitised Fig. 1 (seidel2026.md sec. 2) and **Extended Data Table 1**:
  - the per-night instrument, programme and S/N for each planet;
  - now pushed as `docs/handover/scripts/fe_hires/seidel2026_ed_table1.md`, with the LSF per night.
- **Correction to our seidel2026.md:** it says WASP-121b used "E-4UT x2 (1-UT nights not used)" and WASP-76b "E x2". ED Table 1 lists:
  - WASP-121b: 3 x 1-UT + 2 x 4-UT nights;
  - WASP-76b: 3 x E + 3 x MX.

## 4. fe_hires number for your WASP-121b hydro p1 (rot ~300)
- **viper cannot read Caltech disks.** Two options:
  - **(a) Caltech runs the pushed tool** with `planets/wasp121b_1x_caltech.json` and `--eos /resnick/.../w121r1_1101/anchor/eos_1x.txt` on the last dumps. About 10-12 min per dump on 31 cores.
  - **(b) Copy the last 11 dumps** (~160 MB) of WASP-121b p1 to viper. We run them and reply with the numbers.
- **Comparison point:** our 3x C32 hydro h151 (rot 300) gives -6.03 at the Kp0 row.
