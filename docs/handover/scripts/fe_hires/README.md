# fe_hires: high-resolution Fe I transmission + CCF / Kp-Vsys for AthenaK cubed-sphere dhj dumps

Portable copy (11-02) of the viper pipeline
`/viper/ptmp2/jinma/w121_obs_1003/mhd_final/fe_hires/` (10-07). The method, its validation and
the WASP-121b results are in `FE_HIRES_w121_results.md` (the original write-up; its paths are
viper paths).

**Physics is unchanged.** The WASP-121b defaults reproduce the 10-07 numbers; the check is below.

## Files

| file | role |
|---|---|
| `fe3d.py` | 3-D ray tracer: dump -> per-exposure, per-limb excess-absorption spectra (R = 1e6, 380-790 nm) -> `OUT/spec_<tag>.npz` |
| `ana.py` | CCF with the Fe I template, Kp-Vsys map, Seidel+25 morning/evening segments, contribution function -> `OUT/ana_<tag>.npz`, `OUT/ana_summary.txt` |
| `fe_core.py` | Fe I line list, partition functions, Saha Fe I/Fe II, H- / Rayleigh continuum, template and LSF |
| `csio.py` | dump reader, EOS-table inversion (T, p), cubed-sphere geometry, radial faces |
| `planet.py` + `planets/*.json` | **per-planet hook**: system, instrument, exposures and composition |
| `comptab.cpp` | builds the composition table (n_e, H2, H+, mu) from AthenaK `EOSCompositionModel` |
| `template.py` | regenerates the CCF template with petitRADTRANS 3.4 (not needed to run) |
| `bin_convert.py` | AthenaK `vis/python` reader (v1.1 + v1.2), vendored |
| `run_example.sh` | the WASP-121b reproduction command lines |
| `seidel2026_ed_table1.md` | the paper's per-night instruments, with the LSF of each |

## Requirements

- Python 3.9+, numpy, scipy, numba.
  - Tested with numpy 2.1.3, scipy 1.15.2 and numba 0.61.0.
  - **No petitRADTRANS and no FastChem are needed at run time.**
- g++ and the AthenaK `src/` tree, only to build a new `comptab`.

## Per-planet hook (`planets/<planet>.json`; any key left out takes the WASP-121b value)

| key | meaning | WASP-121b |
|---|---|---|
| `R_star_Rsun`, `R_p_RJ`, `M_p_MJ` | R*, Rp; M_p (static test only) | 1.461, 1.742, 1.17 |
| `a_over_Rs`, `b` | a/R*, impact parameter | 3.7844, 0.160 |
| `P_orb_d`, `omega` | orbital period; rotation rate [1/s] (null = 2 pi/P) | 1.27492504, 5.704026e-5 |
| `Kp_kms` | Kp0 (null = 2 pi a/P sqrt(1 - (b/a)^2)) | null -> 219.2 |
| `limb_darkening` | quadratic u1, u2; [0, 0] = uniform disc | [0, 0] |
| `lsf_fwhm_kms` or `resolution_R` | Gaussian LSF for the spectra AND the template | 4.28 (ESPRESSO 4-UT MR) |
| `nph`, `phase_range` | exposures, equally spaced; null range = T1..T4 | 31, null |
| `rm_mask_phase` | \|phi\| below this left out of the `rmcut` map and the segments | 0.017 |
| `eos_xh`, `eos_yhe`, `comptab` | composition; must match the run's `eos_xh`, `eos_yhe`, `eos_metal_mh` | 0.7188, 0.2420, `data/comptab_3x.bin.gz` |
| `a_fe` | Fe/H by number (not scaled with [M/H]) | 3.2e-5 |
| `pcut_bar` | Fe removed at p < pcut (model top); 0 = no cut | 1e-6 |

- Paper LSFs: ESPRESSO 1-UT HR 2.14, 4-UT MR 4.28, MAROON-X 3.53 km/s.
- `planets/wasp121b_1x_caltech.json` is the 1x composition (eos_xh 0.7381, eos_yhe 0.2485,
  metal_mh 0) for the Caltech 1x hydro runs.
- New composition: build the table, then point `comptab` at it.
  ```
  g++ -O2 -I <athenak>/src comptab.cpp -o comptab
  ./comptab my.bin <xh> <yhe> <metal_mh>
  ```

## Run

```
export OMP_NUM_THREADS=1 NPROC=31          # one process per exposure; ~450 s per exposure
python fe3d.py --planet planets/P.json --eos <run eos_table.txt> --out OUT <dump.bin> [key=val]
python ana.py  --planet planets/P.json --out OUT <tag>
```
- `tag` defaults to the dump basename without `.bin`.
- fe3d keys: `pcut=` (bar), `tcond=2000` (Fe removed below 2000 K), `afe=`, `nph=`,
  `phases=a,b,c`, `dth=`, `bsub=`, `ds=`, `taumin=`, `static=1 [rot=1 Tiso=]` (validation).
- Dumps: AthenaK cubed-sphere `.bin` with `dens`, `velx..z`, `eint` (hydro_w or mhd_w_bcc).
  - Version-1.2 dumps carry the radial faces.
  - Version-1.1 dumps need `--deck run.athinput`, and only for decks with the polynomial stretch
    (no `f_stretch_r_b*`).
- Cost (WASP-121b C32 x nx1 112, 31 exposures, 16 processes):
  - 1000 s wall on 10-07; 1453 s on 11-02 on a busy login node.
  - Peak RSS 3.0 GB.
  - ana.py takes about 1-2 min.

## Inputs

**In git** (`data/`):

| file | size | md5 | what |
|---|---|---|---|
| `fe1_lines.npz` | 386 kB | 1b609d9590a868832d990ffc3e33a789 | 9620 Fe I lines 380-790 nm (vacuum), gf, E_low, gamma_rad, gamma_vdW, from Kurucz `gfall_2600.all` |
| `pf_fe0.npz`, `pf_fe1.npz` | 10 kB each | 4bfbd0eb..., 62b70435... | U(T), Fe I / Fe II, direct level sums over Kurucz `gf2600.pos` / `gf2601.pos` |
| `template_raw.npz` | 5.2 MB | 2b9d56fb43e68c7f49ade1563a81eb46 | pRT 3.4 lbl Fe I transit-depth template (3000 K isothermal, solar Fe, all neutral), unbroadened; = 'raw' of the viper `template.npz` |
| `rayleigh_nodes.npz` | 1 kB | e0ab2893fce28fdd49b18dc00c042fac | pRT 3.4 H2 / He Rayleigh cross sections at the 9 continuum nodes |
| `comptab_3x.bin.gz` | 5.2 MB | deb6f1375c4ec66c545f3ffaf5befb5a | composition table, X 0.7188 Y 0.2420 [M/H] 0.4771 (gunzip = viper `comptab.bin`, md5 fcba1f393b4d5d0b494702923a305129, byte-identical) |
| `comptab_1x.bin.gz` | 5.2 MB | 8fdae489249ef37f8cd8ab055e090e2c | the same, X 0.7381 Y 0.2485 [M/H] 0 (Caltech 1x runs) |

- The template is planet-independent: the CCF divides by sum(template). `fe_core.template_hp`
  applies the continuum removal and the planet LSF at run time.
  - For 4.28 km/s this is bit-identical to the viper `tbhp`.

**Not in git** (only needed to REBUILD the cached inputs above; all on viper):

| file | size | md5 | provenance |
|---|---|---|---|
| `gfall_2600.all` | 6,993,290 B | b572940e38ab3e7c047a0d4588ca7909 | Kurucz Fe I line list with lab gf (kurucz.harvard.edu, linelists/gfall, `gf2600.all` in gfall format). Viper: `/viper/ptmp2/jinma/prt_data/kurucz/`; env `FEHIRES_KURUCZ` |
| `gf2600.pos` | 20,332,368 B | 0f79db1804b5b74454914c8992073d67 | Kurucz atoms/2600 (Fe I, computed lines, all levels), same dir |
| `gf2601.pos` | 20,069,294 B | 586ce7821af8761ab9c261176357277f | Kurucz atoms/2601 (Fe II), same dir |
| `gf2600.all` | 20,332,368 B | 5d39bbd19bd014dcb8458548a36571b7 | Kurucz atoms/2600, not used by the final pipeline |
| `56Fe__Kurucz.R1e6_0.3-28mu.xsec.petitRADTRANS.h5` | 7,294,188,600 B | c7eb0791fc7c36cd6032086030fb780a | pRT 3.4 input data, lbl Fe (MPIA Keeper / Seafile, fetched with curl through the Seafile API on 10-07). Viper: `/viper/ptmp2/jinma/prt_data/input_data/opacities/lines/line_by_line/Fe/56Fe/`. Used only by `template.py` and the pRT opacity validation |
| WASP-121b 3x EOS table `eos_table.txt` | 3,435,515 B | 3a903c5fcf277c642f54f8b56bf7093b | viper `/viper/ptmp2/jinma/w121_obs_1003/diag/3x/eos_table.txt` (the run's general-EOS table dump). Each run must use its own (Caltech 1x: `anchor/eos_1x.txt`) |

- The Kurucz download log on viper is empty. The provenance column is therefore the files'
  standard Kurucz-site location, not a recorded URL.

## Reproduction check (11-02, viper login node, this portable copy)

WASP-121b s10_e15 dump 00166 (10 G MHD, rot 318), the 3x EOS table, `planets/wasp121b.json`;
command lines in `run_example.sh`; log in `run_w121.log`.

| quantity [km/s] | 10-07 viper | portable 11-02 |
|---|---|---|
| Kp-Vsys Gaussian peak at Kp0 (full) | -4.26 | -4.26 |
| Kp-Vsys centroid / single Gaussian | -3.78 / -5.07 | -3.78 / -5.07 |
| map maximum (Kp, v) | 206.2, -3.60 | 206.2, -3.60 |
| morning / evening segment | -3.43 / -6.78 | -3.43 / -6.78 |

- The spectra (`spec_s10_00166.npz`) are **bit-identical** to the 10-07 file (max |diff| 0).
