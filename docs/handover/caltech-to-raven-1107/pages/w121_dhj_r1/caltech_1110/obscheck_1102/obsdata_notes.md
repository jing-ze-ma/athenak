# obsdata.json: where the numbers come from, and caveats (2026-11-02)

Builder: `build_obsdata.py` (it re-creates `obsdata.json` from the numbers transcribed in it plus the data files in `lit/`).
Sources: `lit/src_<arxiv>/` (arXiv e-print, .tex or PDF converted to text.txt), `lit/zenodo_*`, `lit/nature_me22/` (Nature Astron. Source Data).

Conventions used in the JSON:
- Phase: 0 = mid-transit, 0.5 = mid-eclipse.
- Fp/Fs: ppm.
- `hotspot_offset_deg_east`: positive means the peak comes before eclipse (eastward). The original sign convention is stored per entry.
- Velocities: km/s, in the planet rest frame, negative = blueshift.
- Values: `{v, e}`, or `{v, ep, em}` for asymmetric errors.
- Status flags: tabulated / derived / digitised.

## WASP-121b emission
| id | status | provenance / caveats |
|---|---|---|
| MikalEvans2023 (2301.03209) | tabulated | Text and abstract, white light only, preliminary dipole fit. Superseded by Evans-Soma 2025. Nightside = dipole minimum. |
| EvansSoma2025 (2506.01771) | tabulated (Zenodo 15377434) | The dayside spectrum is phase bin n_phi=34 (0.426-0.451) and the nightside is n_phi=18 (0.925-0.950), as in their Fig. 1. These are not exactly phases 0.5 and 0. The band-averaged phase curve is derived here as the unweighted channel mean per detector. White-light offsets (2.95/2.39 deg E) come from their Supplementary Table 1 as quoted by Splinter 2025. Stellar parameters are in their Supplementary Table 9, which was not read. All 36x349 arrays are in `lit/zenodo_15377434`. |
| Splinter2025 NIRISS SOSS (2509.09760) | white light tabulated, spectra digitised | The white-light nightside is not tabulated. The nightside spectrum depends on the reduction (exoTEDRF vs NAMELESS differ by up to about 250 ppm below 1.5 um). The digitised spectra are coarse; for the dayside, use the Kahle Zenodo NIRISS spectrum instead. |
| Kahle2026 (2606.21855) | tabulated (Zenodo 20767846) | Panchromatic eclipse spectrum: NIRISS (Pelletier 2026 reduction), NIRSpec (Evans-Soma bins converted to eclipse phase), MIRI LRS 5-12 um. The MIRI white-light depth is 6188 +- 36 ppm. Their retrieval allowed instrument offsets of about 240-260 ppm. Preprint, not checked for acceptance. |
| Sing2025 (2501.03844) | tabulated | Planet Kp from NIRSpec emission lines: NRS1 212.7 +- 1.8, NRS2 218.2 +- 1.0 km/s. The authors read the difference as a hint of winds. This entry also holds the newest stellar parameters. |
| MikalEvans2022 HST (2202.09884) | tabulated (Source Data) | Spectra are at 16 phases, plus Tb vs phase from Supplementary Table 6. The broadband offset is pi - c2 = 6.4 (+1.7/-2.3) deg E from the preferred 2nd-order fit; the 1st-order fit gives 4.7. Some errors are garbled in the arXiv PDF. The binned broadband light-curve points near transit and eclipse are not planet flux (they are flagged). |
| MikalEvans2020 (2005.09631) | tabulated | G141 eclipse spectrum from 5 eclipses (includes the Evans 2017 visit, which was not transcribed separately). |
| Changeat2024 (2401.01465) | no numbers | Only a model-dependent hotspot (fixed grid; best 30 deg E, 50 deg wide). Noted inside the MikalEvans2022 entry. |
| Morello2023 Spitzer (2307.00669) | tabulated | The paper's offset sign is + = after eclipse, so it was flipped here (3.6 um: -5.9 deg, i.e. west). The nightside minimum at 3.6 um sits 31 +- 6 deg before transit. IRAC band edges are nominal and not taken from the paper. |
| Davenport2025 (2503.12521) | tabulated | Same Spitzer data, offsets about 0. "Amplitude" is the semi-amplitude. The temperatures in the table do not match the text (2779 vs 2543 K). The stellar-parameter table shows R* = 0.91, which is probably a typo. Use the fluxes. |
| Dang2025 4.5 um (2408.13308) | tabulated | Uniform reanalysis. Values are rounded to 100 ppm. |
| Daylan2021 / Bourrier2020 TESS | tabulated | Daylan's offset is reported as westward and was sign-flipped. |

## WASP-76b emission
- **May & Stevenson 2021 (2107.03349):** eclipse depths, phase amplitudes, offsets and Tb are tabulated. The nightside Fp/Fs is **derived** as eclipse - 2 x amplitude, assuming the amplitude is a semi-amplitude, which gives about 930 and 800 ppm. This agrees with Dang's F_min of 900 +- 300 at 4.5 um, but verify before tight use. The paper does not state the sign of the offset (SPIDERMAN lo0). The authors call the 3.6 um data less reliable.
- **Dang 2025:** the 4.5 um depth of 3400 ppm is well below May's 3729 ppm, possibly because of the companion-dilution correction (unverified).
- **HST G141 eclipse:** both Fu 2021 (Table 7, dilution-corrected, half-width 0.0093 um) and Edwards 2020 (Table 3, larger dilution factor, depths about 10-20 % higher) are included. The two disagree systematically, so treat the difference as a systematic floor.
- **JWST:** no published WASP-76b eclipse or phase curve was found. GO 5268 (NIRSpec G395H phase curve) is unpublished.

## Fe I high-resolution transmission
WASP-76b:
- **Ehrenreich 2020:** the summary numbers are tabulated. The time series (14 binned points, velocity and differential depth) is **digitised** from Extended Data Fig. 7, with a digitisation error of about 0.4 km/s and 30 ppm. Mid-transit (-0.2 to +0.7 h) is excluded. Kp = 196.52.
- **Kesseli & Snellen 2021 (HARPS):** first half +0.04 +- 1.47, second half -6.9 +- 0.75 km/s.
- **Kesseli 2022 (ESPRESSO):** Fe I RV1 -3.89 +- 0.48 and RV2 -10.26 +- 0.24 km/s; amplitudes 440 +- 47 and 833 +- 83 ppm.
- **Gandhi 2022 and 2023:** retrieved dVsys per quarter (the wind shift of the whole spectrum, not Fe-only CCF centroids).
- **Pelletier 2023:** qualitative only (kinked trail; Kp 180.7, Vsys -7.0 from a fit to a kinked trail). No per-phase numbers.

WASP-121b:
- **Seidel 2025 (4-UT):** Fe I template morning -4.12 +- 0.15 and evening -6.90 +- 0.11 km/s (table). The per-exposure points are **digitised** from Fig. 4a. The colour-to-template mapping was inferred from the segment means, because the caption has no key.
- **Borsa 2021:** first/second half -2.80 / -7.66 km/s and CCF contrast 1050 / 970 ppm. This uses a stellar mask that is about 76 % Fe I lines, not an Fe-only template.
- **Gandhi 2023:** dVsys per phase range.
- **Prinoth 2025:** no Fe numbers (only Ti and V are time-resolved).
- **Hoeijmakers 2024 (Mantis IV):** this is dayside Fe I *emission*, not transmission. The peak-offset sign was converted on an assumption from their eq. 3.

Not found or not done:
- Seidel 2026 (2606.22455, magnetic fields) was skipped because it is derived, not new observed data.
- The KPF WASP-76b paper (2511.14175) and Costa Silva 2024 (WASP-76b dayside Fe emission) were outside the list and were not extracted.
- Bourrier 2020 HARPS Fe halves for WASP-121b were not extracted.
- Gapp 2026 (JWST rotational transits, morning/evening limbs) is transmission and was not extracted.
