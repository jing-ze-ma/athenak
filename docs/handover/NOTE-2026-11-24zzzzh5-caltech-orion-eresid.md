# NOTE 2026-11-24zzzzh5 Caltech -> orion: per-cycle energy and mass residual of the production histories (zzzzh4, zzzzh5)

Read-only. Files: caltech_results/ERESID/ (run.log = all tables, outliers_<star>.txt, eres.py, inject.txt, inject.py).

**Formula** (source a000d803c, he_star_vet.cpp HeStarHist comment): d Etot/dt = L_in - L_top + P_esrc - P_sponge - E_hflx (Etot = gas IEN incl. kinetic and static potential + radiation; E_hflx = IEN through the physical x1 faces); d M_tot/dt = -M_hflx. Per row c: R = Etot(c) - Etot(c-1) - dt x rate(c), RM = M_tot(c) - M_tot(c-1) + dt M_hflx(c), dt = time(c) - time(c-1). L_bot - L_top (your wording) is not used: the source uses L_in (imposed wall flux); L_bot is a diagnostic.
**Rate row.** The rates at row c are the step that ended there (the dt column of row c equals the time difference on every row of every star; 0 rows with dt <= 0). Taking them at row c-1 or the mean gives scatter 1e-12 to 1e-11 of Etot in bsg, m35mms, m35zams, t9mlt. (For t9l62 the row c-1 scatter came out exactly 0: the running median of a smooth monotone series is its own value, so that number is meaningless and was not used.)
**Baseline and rule:** running median of 51 rows inside each chunk; outlier = |R - median| > 10 x MAD (MAD of |R - median| over the star). Rows are every cycle (dcycle 1, dt 0), all consecutive; chunks cut at the headers.
**Precision.** hst format %24.16e = 17 significant digits (round-trips a double), so the print does not limit a 1e-10 bound. The limit is the double sum over the cells: MAD 1e-15 relative.

| star (chunks) | rows used | t range [s] | MAD of R/Etot | E outliers (10 MAD) | > 1e-10 | largest \|R - med\|/Etot | MAD of RM/M_tot | M outliers | largest |
|---|---|---|---|---|---|---|---|---|---|
| bsg (1; sdorder3 95813cf85, c0-23850) | 23859 | 0 - 3.96e6 | 1.76e-15 | 0 | 0 | 1.12e-14 | 1.71e-15 | 0 | 1.21e-14 |
| t9l62 (2; c4 a000d803c, c4speed d979d5268 from c3714) | 4739 | 0 - 4.49e5 | 1.25e-15 | 0 | 0 | 6.9e-15 | 1.24e-15 | 1 (c1) | 4.0e-14 |
| m35mms (2; sdorder3, c4speed from c16341) | 18363 | 0 - 2.14e5 | 1.61e-15 | 0 | 0 | 1.04e-14 | 1.70e-15 | 0 | 1.10e-14 |
| m35zams (2; sdorder3, c4speed from c50817) | 55452 | 0 - 1.93e5 | 7.7e-16 | 0 | 0 | 5.2e-15 | 7.4e-16 | 0 | 5.5e-15 |
| t9mlt (1; c4speed) | 561 | 0 - 5.4e4 | 1.21e-15 | 0 | 0 | 5.4e-15 | 1.23e-15 | 1 (c1) | 2.6e-14 |

- **No row of any star has an energy or mass step above 1e-13 of the total; the only outliers are the first step in t9l62 (-4.0e-14 M_tot, c1) and t9mlt (+2.6e-14, c1).** Those rows have retested 2.5e7 / 5.1e7 (medians 4.0e7 / 5.2e7); the floor counters of the same cycle are 0 in t9l62 and at or below typical in t9mlt (floor_pay 9.3e25 against 1.28e26).
- Rows are still being added (bsg, t9l62, m35mms, m35zams, t9mlt running): numbers are at the time of the read (11-24, shortly after 21:07 PST).
- **Sensitivity checked by injection** (inject.txt; 20 random single-row jumps in copies of chunk 0): 1e-9, 1e-10, 1e-12 and 1e-13 are flagged 20/20 in E and M for every star; 1e-14 is at the edge (0-3 of 20; m35zams 20 of 20). So KEYTRIM's 1e-9 to 5e-8 jumps would be seen in these files; none is.
- What it cannot see: a per-cycle jump below 1e-13 of Etot (bsg: 1.6e33 erg), or jumps that cancel inside one row. A seam error that conserves the total Etot and M_tot is invisible here.
- **Switch rows** (first row of the new binary minus the last old row, the rates of the first new row): energy residual t9l62 +1.8e-15 (c3713 -> 3714), m35mms +2.1e-15 (c16340 -> 16341), m35zams -9.3e-16 (c50816 -> 50817); mass +1.1e-15, -4.7e-15, +3.7e-16. All round-off (your 3.7e-16 of m35zams reproduced). The raw dEtot/Etot across them is one ordinary step (t9l62 -1.6e-7, m35mms +2.8e-9, m35zams -8.3e-11), not a jump. bsg and t9mlt have one chunk.
- floor_from_ekin = true in all five decks (and fofc, fofc_retest, floor_pay = rad_column).

**FOFC counts per cycle.** The logs have no `fofc_report` lines; they have `fofc_retest: retested N unflagged M` (floor_audit), cumulative in the process, differenced per cycle here. Source: retested = flagged trial cells re-tested, summed over the stages of the cycle, so it counts cell-stages, not distinct cells; unflagged is 0 in every cycle of every star.
| star | retested per cycle: median | p99 | max | other counters (median per cycle) |
|---|---|---|---|---|
| bsg | 7.59e6 | 8.66e6 | 9.26e6 | c2p_dfloor 5.3e5, skip 2.6e5, floor_pay 0 (45 cycles > 0, max 1.4e17 erg) |
| t9l62 | 4.05e7 | 4.20e7 | 4.72e7 | floored_update 1.1e6, floor_pay 3.3e25 erg, pos_gas 1.2e6, floorcap v max 23 (175 cycles > 0) |
| m35mms | 5.14e6 | 6.29e6 | 6.44e6 | floored_update 0 (8109 cycles > 0, max 65), floorcap v 0 (2382 cycles > 0, max 2.5e5) |
| m35zams | 1.12e6 | 1.41e6 | 1.42e6 | all 0 (c2p_dfloor max 2) |
| t9mlt | 5.15e7 | 5.27e7 | 5.72e7 | floored_update 8.1e6, floor_pay 1.3e26, pos_gas 8.6e6 |

There are no outlier rows whose counts could be compared with the typical values (only the two c1 rows); floor counts at those rows are in outliers_*.txt.
floor_pay energies (1e25-1e26 erg per cycle in t9l62 / t9mlt) are 1e-22 of Etot and do not affect the residual.

-- Caltech
