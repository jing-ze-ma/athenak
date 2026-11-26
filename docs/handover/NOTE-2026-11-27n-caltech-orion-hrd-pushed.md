# NOTE 2026-11-27n caltech -> orion: dark HR diagram pushed (the USER: "push the dark version")

docfig-movies-1126, commit 32e826192 (only these four files added; nothing else changed):

| file | bytes |
|---|---|
| hrd_3d_dark_vet.png (2300 x 1520) | 447,656 |
| hrd_3d_dark_vet.pdf | 86,545 |
| hrd_3d_dark_vet_transparent.png | 429,916 |
| sims_3d_vet.csv | 14,202 |

Made with a copy of viper's hrd_3d/make_hrd.py and background (Marchant & Bodensteiner 2024 Fig. 2): only our four stars plotted (code AthenaK; no code legend, no other codes' points), labels M20 BSG, M35 MMS, M35 ZAMS, M56 LBV. In the csv the three older AthenaK rows (B-star box, He-star box, He4 planned) have plot = no and our four rows are added at the end; the other rows are viper's, unchanged. The light version and hr_tracks_zoom are not pushed (USER: only the dark version).

Caption text: "Time means over the run window of each Caltech AthenaK simulation: log L = log10 of the mean of the luminosity through the top of the domain (hst L_top x 4 pi/Omega); T_eff = mean over dumps of the median over columns of the comoving-flux effective temperature at each column's own tau_R = 2/3. Windows: M20 BSG 30-57 d, M35 MMS and M35 ZAMS 5-10 d, M56 LBV 17.25-22.25 d (still evolving). Time spreads (standard deviation) are smaller than the markers (at most 0.0065 dex in L and 0.0023 dex in T_eff)."

Values (log L, T_eff [K]): M20 BSG 5.176, 15646; M35 MMS 5.444, 35124; M35 ZAMS 5.238, 41870; M56 LBV 5.830, 22490 (NOTE 2026-11-27b2).

-- Caltech
