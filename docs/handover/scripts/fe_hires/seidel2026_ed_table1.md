# Seidel et al. 2026, Extended Data Table 1 (observations), the five Caltech planets

Source: text layer of the ESO PDF (viper /viper/ptmp2/jinma/dhj_planets_1101/src/eso2606a.txt,
p. 52), copied verbatim; LSF column added from the paper's Methods 5.1.3 ("broadened to the
approximate line-spread functions of the instruments: 2.14 and 3.53 km/s for ESPRESSO (1UT-mode HR)
and MAROON-X ... 4UT MR ... 4.28 km/s"). E = ESPRESSO 1-UT HR, E-4UT = ESPRESSO 4-UT MR,
MX = MAROON-X. S/N = minimum-maximum per night.

**Not in the paper (or not available to us):** per-night phase ranges / exposure lists, the
per-night stacking weights, and Supplementary Table 1 (planet parameters; not in the ESO PDF).

| planet | date | inst. | LSF [km/s] | PI, programme | S/N | Teq [K] |
|---|---|---|---|---|---|---|
| WASP-76b | 2018-09-02 | E | 2.14 | Pepe, 1102.C-0744 | 13.1-83.7 | 2229 |
| | 2018-10-30 | E | 2.14 | Pepe, 1102.C-0744 | 10.5-67.5 | |
| | 2019-10-18 | E | 2.14 | Gibson, 0104.C-0642 | 29-55 | |
| | 2020-09-03 | MX | 3.53 | Pelletier, GN-2020B-Q-122 | 20.0-68.4 | |
| | 2020-09-12 | MX | 3.53 | Pelletier, GN-2020B-Q-122 | 17.6-63.1 | |
| | 2021-10-28 | MX | 3.53 | Debras, GN-2021B-Q-138 | 29.4-99.0 | |
| WASP-121b | 2021-01-26 | E | 2.14 | Gibson, 0106.C-0516 | 27-37 | 2358 |
| | 2021-03-04 | E | 2.14 | Gibson, 0106.C-0516 | 21-36 | |
| | 2019-01-06 | E | 2.14 | Pepe, 1102.C-0744 | 21-36 | |
| | 2018-11-30 | E-4UT | 4.28 | 60.A-9128 (Comm.) | 23.1-163.3 | |
| | 2023-09-23 | E-4UT | 4.28 | Seidel, 111.24J8 | 26.6-174.5 | |
| WASP-178b | 2021-05-03 | E | 2.14 | Pepe, 1104.C-0350 | 11.3-48.6 | 2398 |
| | 2021-07-09 | E | 2.14 | Pepe, 1104.C-0350 | 9.3-37.7 | |
| TOI-1518b | 2022-08-13 | MX | 3.53 | Parmentier, GN-2022B-Q-128 | 23.4-64.9 | 2491 |
| | 2023-10-19 | MX | 3.53 | Parmentier, GN-2023B-Q-127 | 27.4-80.2 | |
| | 2024-06-26 | MX | 3.53 | Parmentier, GN-2024A-Q-130 | 24.0-65.1 | |
| WASP-189b | 2021-06-04 | E | 2.14 | Prinoth, 107.22QF | 107.7-432.6 | 2638 |
| | 2022-04-03 | MX | 3.53 | Pelletier, GN-2022A-FT-208 | 84.1-231.6 | |
| | 2022-06-02 | MX | 3.53 | Pelletier, GN-2022A-FT-208 | 84.7-233.5 | |

Notes:
- WASP-121b: the table lists 3 ESPRESSO 1-UT nights and 2 4-UT nights, but Methods 5.1.1 (p. 12)
  says: "For WASP-121 b, we only used the data taken in 4-UT mode". So only the two E-4UT nights
  (LSF 4.28) enter Fig. 1 (corrected 11-02; the earlier note here was wrong).
- WASP-76b: Methods 5.1.1: "we discarded the night of 18 Oct 2019 for WASP-76 b" -> 2 E + 3 MX used.
- WASP-76b and WASP-189b mix ESPRESSO and MAROON-X nights (two LSFs). The paper stacks nights in
  Kp-Vsys after a per-night system-velocity correction, "based on S/N per night".
- The 2018-11-30 / 2023-09-23 4-UT nights are partial transits (egress / ingress; Seidel+2025).
