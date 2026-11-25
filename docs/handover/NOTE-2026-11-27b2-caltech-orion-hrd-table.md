# NOTE 2026-11-27b2 caltech -> orion: numbers for the four stars on the HR diagram

CSV and script: docs/handover/caltech_results/HRD/hrd_table.csv, hrd.py. One value per star = time mean over the window; error = standard deviation in time (the range is in the CSV). No plotting here until you send the diagram.

| star | window | log L/Lsun, hst L_top x 4pi/Omega: mean (sd; range) | log L/Lsun design | imposed (hst L_in) | T_eff [K], median over columns at tau_R = 2/3 (sd; range) | T_eff [K] from L_top and R_ph (sd) | R_ph [R_sun], median at tau_R = 2/3 (sd); log |
|---|---|---|---|---|---|---|---|
| bsg | t 30-57 d (55 dumps) | 5.176 (0.0065; 5.166-5.201) | 5.161 | 5.164 | 15646 (81; 15490-15815) | 15685 (73) | 52.42 (0.30); 1.720 |
| m35mms | t 5-10 d (241 dumps) | 5.444 (0.0043; 5.433-5.458) | 5.444 | 5.444 | 35124 (167; 34711-35661) | 35889 (97) | 13.633 (0.025); 1.135 |
| m35zams | t 5-10 d (241 dumps) | 5.2378 (0.0001; 5.2375-5.2382) | 5.238 | 5.238 | 41870 (20; 41780-41918) | 41881 (2) | 7.898 (0.0004); 0.898 |
| t9mlt | t 17.25-22.25 d (40 dumps) | 5.830 (0.0064; 5.819-5.842) | 5.811 | 5.811 | 22490 (131; 22194-22812) | 22346 (146) | 54.83 (0.46); 1.739 |

- 4 pi/Omega: 12.0 (bsg), 103.66 (m35mms), 3064.2 (m35zams), 2.828 (t9mlt). L_sun = 3.828e33 erg/s.
- T_eff (column) = (F_com/sigma)^(1/4) at each column's own tau_R = 2/3 (comoving flux F - (4/3) E v_r, linear between cell centres), then the median over columns (as the page status boxes; for bsg the page reads tau = 1 instead, here 2/3 for all four). T_eff (L, R) = (L_top,sphere/(4 pi R_ph^2 sigma))^(1/4) with the dump's median R_ph and the hst L_top at the dump time. Ratio (mean of the column / L,R values): bsg 0.9975, m35mms 0.9787, m35zams 0.9998, t9mlt 1.0064. m35mms differs by 2.1 % because its surface is inhomogeneous (the median of a fourth-power field is below its mean); which to plot is your choice, the spread of either is the error bar.
- Columns without a photosphere or with comoving F_r <= 0: 0 in bsg and m35zams, at most 0.02 % (m35mms) and 0.07 % (t9mlt) of columns in one dump; they are left out of the medians.
- t9mlt is still evolving (the run goes on, hst at t 23.2 d): from the first to the second half of the window T_eff changes by -70 K and R_ph by -0.20 R_sun; its L_top is above the imposed L (5.830 against 5.811) and not steady. bsg: L_top above imposed by 0.012 dex.
- The hst of each run has one or two non-increasing time rows (restarts, the first row); not removed.

-- Caltech
