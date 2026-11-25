# NOTE 2026-11-24zzzzj8 Caltech -> orion: photospheric radius and T_eff across the rel-f restarts (your zzzzj6 item 2 / zzzzj7)

Read-only on the arms of zzzzj4 (b1_1124/rc/<star>/{R1,R25a,R50,C50}) and the production dumps. **Radius passes the 0.05 % (5e-4) gate at the switch in all five stars.** L_top / L_bot: NOTE zzzzj6-caltech-orion-lcont, not here. Files: caltech_results/RCONT/.

**Which dumps exist.** Arm dumps are the END dumps only: c+1 (R1, rel-f), c+25 (R25a, rel-f), c+50 (R50 rel-f, C50 running binary). There is no c+10 dump, and no c+1 / c+25 dump of the running binary. The production dump nearest before the restart has the SAME cycle as the rst in m35mms (c22144), bsg (c29127), t9l62 (c6332), t9mlt (c2683), so R1 - P there is exactly the first rel-f step. m35zams has none at c50816: its nearest dump is c50683, 134 cycles before R1.
**Tools.** (a) your rsg_rph/rph.py (md5 77340fc0), ported to the spherical wedge (rph_wedge.py; changed: vet instead of m1 names, wedge solid angles, hst L_top x 4 pi/Omega, table header line found by pattern, more digits; kappa(T_rad) and tau from the cell inner faces are unchanged). Columns: R23_mean and T_eff(R23_mean). (b) our reducer (rph_caltech.py; kappa(T_gas), tau at cell centres): RA = radius where the solid-angle mean tau is 2/3, RM = median over the columns of their own tau = 2/3 radius (wedge page), TM = median column T_eff from the comoving flux. Both agree on every conclusion. The bsg page's definition (mean tau = 1) is in ana.txt.

| star (bins) | first step R1 - P: R23_mean | same, RA | mean rel-f rate per step c1-c25 (R23_mean) | R50 - C50: R23_mean / RA / TM | R50 - P (= C50 - P): R23_mean | production dump-to-dump rms / max (R23_mean) |
|---|---|---|---|---|---|---|
| m35zams (float32) | -8.3e-7 (134 cycles, not one step) | -3.5e-6 | -4.6e-8 | -3.6e-10 / -3.5e-10 / +1.6e-7 | -3.3e-6 | 3.5e-5 / 6.7e-5 |
| m35mms (float32) | -2.3e-5 | -2.2e-5 | -2.4e-5 | -6.6e-9 / +4.3e-9 / -1.2e-4 | -1.22e-3 | 1.8e-3 / 3.5e-3 |
| bsg (float64) | -6.0e-6 | -2.2e-6 | -4.3e-6 | -6.1e-8 / -1.1e-8 / +8.1e-7 | -1.4e-4 | 3.9e-3 / 8.6e-3 |
| t9l62 (float32) | +1.6e-4 | +1.3e-4 | +1.6e-4 | +4.5e-8 / +5.3e-9 / -4.1e-5 | +7.9e-3 | 3.1e-2 / 4.1e-2 |
| t9mlt (float32) | +2.2e-4 | +2.2e-4 | +2.2e-4 | +3.9e-9 / +6.5e-9 / -8.4e-6 | +1.10e-2 | 1.3e-2 / 2.7e-2 |

- **Switch (first step):** the largest |change| over all radius definitions is 2.3e-4 (t9mlt Rsh23), below 5e-4. For R23_mean the first step equals the star's own mean rate over c1-c25 to a few per cent in m35mms, t9l62, t9mlt (table). The other definitions scatter more around that rate (t9mlt RM: +9.9e-5 against +2.5e-4 per step; bsg RA: -2.2e-6 against +2.1e-6), always below the gate. m35zams has no dump at the rst cycle: R1 is -8.3e-7 (R23_mean) from the dump 134 cycles earlier and +7e-7 from the production curve. First-step T: at most 1.2e-4 in any of rph.py T_eff, RA-flux T and TM.
- **R50 against C50** (same rst, same keys): radius equal to <= 6.1e-8 in every star and radius definition, T within 1.2e-4 (m35mms TM; rph.py T_eff within 4.5e-6). Differences below ~1e-7 are at the single-precision rounding of the float32 dumps and cannot be read as physics, except bsg (float64).
- **R50 against the production dump (50 cycles later)** is above 5e-4 in m35mms, t9l62 and t9mlt. It is the evolution over the 50 cycles: C50 of the running binary gives the same number to 1e-7. It is not a jump at the switch. (The dump-to-dump numbers of t9l62 and t9mlt are mostly their expansion, not noise; m35mms and bsg are variable.)
- **R50 against the production curve** (cubic through the two dumps before and two after, production binary): R23_mean m35zams -1.0e-6, m35mms -5.0e-5, bsg -5.0e-4 (RA -5.6e-4), t9l62 -4.6e-6, t9mlt +4.8e-5; C50 gives the same values. Other definitions reach -8.7e-4 (bsg RM) and +7.7e-4 (t9mlt RM), and -2.0e-3 in T (m35mms RA-flux T), again the same for C50 (T to 1.2e-4); the leave-one-out rms of the cubic at the full dump spacing is 1e-3 to 4e-3 for those series (bsg radius 2.3e-3), so these residuals are within the interpolation error, not a difference between the binaries.
- **Not tested:** continuity inside c2-c24 and c26-c49 (no dumps), and the dump-to-dump noise of the 50-cycle window; whether R50 stays equal to the production beyond 50 cycles.

-- Caltech
