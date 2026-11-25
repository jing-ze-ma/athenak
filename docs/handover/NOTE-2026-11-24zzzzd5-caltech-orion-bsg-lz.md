# NOTE 2026-11-24zzzzd5 Caltech -> orion: bsg merge fall-back and L_z (zzzzd5)

Read-only, hst + logs + three hydro_w dumps. Files: caltech_results/BSG_LZ/ (lz.txt: every 500-cycle window and MLATDIAG line; step.txt per-cycle; groups.txt; scripts). Wedge values (theta 60-120, phi 0-60 deg), not scaled, g cm2/s.

**1. Lz_tot.** The hst has Lz_tot and Lz_hflx, no other component. With residual = change of Lz_tot + integral of Lz_hflx dt (source he_star_vet.cpp: d Lz_tot/dt = -Lz_hflx, the L_z leaving through the boundaries; the column is 0 without its module), the column accounts for the drift to 1e-3 - 3e-3 of it in cycles 0-6000:
| window (cycles) | c1 (95813cf85) dLz | c1 residual | T2 (262affd6f) dLz | T2 residual |
|---|---|---|---|---|
| 0-500 | -7.43e36 | +2.3e34 | -7.39e36 | +2.3e34 |
| 4500-5000 | -1.29e37 | +2.6e34 | -1.24e37 | +2.4e34 |
| 5500-6000 | -1.96e37 | +1.5e34 | -1.74e37 | +1.2e34 |
| 6000-6500 | +4.02e38 | +4.15e38 | +5.46e38 | +5.54e38 |
| 6500-7000 | -1.69e39 | -1.66e39 | -1.73e39 (to c6902) | -1.71e39 |
| 7000-7500 / 7500-8000 | -2.0e38 / -1.2e38 | -1.1e38 / +6e35 | run ended | |
| 8000-8500 / 8500-9000 | -4.9e38 / +6.25e39 | -3.7e38 / +6.31e39 | | |
| 9000-9500 / 9500-10000 | +8.3e38 / -6.8e38 | +8.3e38 / -7.2e38 | | |
| 10000-10500 / 10500-11000 | -1.6e38 / -4.8e38 | -1.8e38 / -4.5e38 | | |

Per cycle, c1's residual first exceeds 1000 x its c1000-5000 median (3.2e31) at **c6050**, T2's at **c6061**; none before c6000 (max 6.9e32 / 1.1e33). **nfb in c1 is 0 at c6000 and 101260 at c7000, 105044 at c8000, 575785 at c9000, 599403 at c10000, 667181 at c11000; T2's is 0 at c6000 and its last line; T2 ends at c6902, so no T2 nfb after 6000 exists.** MLATDIAG prints every 1000 cycles, so nfb cannot be dated inside a window. The c1 windows with the larger nfb increase (6000-7000: 1.0e5; 8000-9000: 4.7e5) have the larger residuals, but T2, whose nfb after c6000 is unknown, jumps in the same window by the same size. This table doesn't show the fall-back as the cause.

rc3 (af836f5a3; chunk 4 from c39675 on dflfix 5b6159ad0): the mean Lz_hflx is 1e29 in c0-500 and 0 in c500-39500, so the residual there is just dLz; it is nonzero from the c39500-40000 window on (residuals -2.1e38, then -3.9e36, -4.3e37, -4.5e36, +6.1e36 per 500 cycles). Lz_tot changes: c0-500 +7.8e33; c5000-5500 +5.2e37; 10000-10500 -1.0e38; 20000-20500 +1.2e38; 30000-30500 -7.0e38; 40000-40500 -1.3e38; 50000-50500 -2.1e38. Its nfb (counted per process, reset at every restart) is 1.95e6 at c1000, 1.25e7 at c2000, 6.9e7 at c6000, 2.9e9 at c38000; 5.1e7 in the last rc3 chunk's print at c39000; **3.8e5 at c40000 and 4.0e5 at c55000 after the dflfix restart**.

**2. Scale (dumps: c1 00037, c9601; T2 00027, c6902; rc3 00114, c55547).** Lz_tot first row 1.636e35; c1 extremes 4.9e39; net Lz_tot c1 3.4e39 at the last row (c11575), T2 -1.31e39 (c6902), rc3 3.5e39 (c55547). Sum of |rho r sin(theta) v_phi dV| over cells: 1.54e48 (c1), 1.33e48 (T2), 2.10e48 (rc3); **the residual jumps (up to 6e39) are 4e-9 of it, and the net Lz is also 1e-9 of it** (so a relative drift against the initial Lz is not meaningful). Merge shells: r < 2.435e12 gross 6.73e47 (c1), net -4.9e45; 3.826-5.566e12 gross 4.06e44, net -1.4e44, mass 9.1e25 g of 2.90e31.

**3. MLATDIAG.** Two lines per print (states inside the cycle; values differ in the 3rd digit); L_x, L_y, L_z are Cartesian cell-centre sums and L_z there is not Lz_tot of the hst (c6000: -3.3e40 against -1.39e38; different arm/sum, cause not checked), so no drift is quoted from it. **chk is 0 in every line: merge_lat_check is not set, so the group sums before/after the average are not measured.** nfb at c6000 0 / c7000 101260.

**4. Where are the fall-back groups.** No per-cell print exists; only the count nfb (source hydro_mergelat.cpp, incremented in Average, kind 0, i.e. every stage and again after the radiation update; not per group). Source fact: the code falls back when erot > 0.5 (E_group - E_kin,bulk - E_pot), and that bracket still contains the group's own rotation and non-rigid energy; so it fires at about erot > E_th + residual, not at erot > 0.5 E_th (the header comment says half of the thermal energy).
Dumps are the POST-merge state (rho is identical inside every group: spread 0), where a fall-back group has zero rotation, so a past fall-back is not visible in them. What they show (2x2 lateral groups, f = 2 shells only; Erot = 0.5 omega.L, E_th = e_int V):
| dump | groups r < 2.435e12 | 3.826-5.566e12 | Erot > 0.5 E_th | Erot > E_th | rho / e_int range of the 0.5 set |
|---|---|---|---|---|---|
| c1 c9601 | 786432: 0 | 1310720 | 930 (r 3.84-4.40e12) | 0 (max 0.998) | 6.2e-14-9.9e-11 / 0.15-394 erg/cm3 |
| T2 c6902 | 0 | | 5811 (r 4.10-4.42e12) | 0 (max 0.999) | 7.6e-14-5.7e-11 / 0.19-116 |
| rc3 c55547 | 0 | | 103 (r 3.84-3.95e12) | 0 (max 0.950) | 1.0e-13-6.0e-12 / 0.22-12 |

No group in r < 2.435e12 (f = 2 inner shell) passes 0.5 E_th in any of them; all sit in the outer f = 2 shell just above its inner edge. Groups with all 4 cell velocities equal to float32 (zero rotation): 1 / 3 / 0 (not a count of fall-backs).

**Version.** An uncorrected version of this NOTE was public in f57377f32 and replaced in 40853393d; do not use it. Wrong there: it said the sign of Lz_hflx was not looked up (the source gives d Lz_tot/dt = -Lz_hflx), that the closure is 1e-3 (1e-3 - 3e-3), and that rc3's Lz_hflx is 0 until chunk 4 (it is 1e29 in c0-500, 0 in c500-39500, nonzero from the c39500-40000 window on). The tables were not affected.

-- Caltech
