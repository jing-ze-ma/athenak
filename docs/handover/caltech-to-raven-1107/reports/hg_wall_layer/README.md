# HG wall layer (L_bot 4.7 % above L_in, +0.23 %/100 d, odd-even v_r in the first cells)

Source: `$J/hg_relax_1106/` (read-only analysis of HG prodG dumps + hst, 11-06; summary in `HANDOVER-2026-11-06-pm-caltech.md` "HG relaxation", "Open: wall layer").
All files here are unchanged copies (figures, results*.txt, analysis scripts; the .npz/.npy caches were not copied).

Where each number is:
- L_in = 1.4037e38: `results2.txt` line 1.
- L_bot = 1.4695e38 (= 1.047 L_in, i.e. +4.7 %), trend +0.230 +- 0.021 %/100 d (|s|/SE 10.8), 4 x 53-d windows from T0 = 328 d to 573 d: `results2.txt` line 31 (a second fit with another T0 at line 52: 1.4655e38, +0.334 %/100 d).
- Odd-even mean v_r in the first wall cells: `results3.txt` "WALL:" block (lines 166-177): column `mvr` (cm/s) for i = 0..7: -2.29e3, +1.53e3, -8.5e2, +1.65e3, +3.2e2, +1.19e3, +2.9e2, +7.2e2 (alternating, i.e. +-10-20 m/s), with the per-cell Lr / Lrcom / Lfl / LE split (lab Lr also alternates 0.92 / 1.76 / 1.17 / 1.57 e38). Produced by `an3.py` ("# wall layer" section).
- Budget: energy closes to 0.09 % of L over 328-573 d (HANDOVER); no source for the wall excess found.
- Figures: fig1_series.png (hst series), fig2_profiles.png, fig3_flux.png (flux split vs r), fig4_timescales.png, look_ts.png.

Keys tried: NONE. No run or key was tried against the HG wall layer on Caltech; it is listed as open (untested idea: numerical diffusive flux at depth, from the dump mean flow disagreeing with continuity).
Related but different: `$J/ke_growth_1106/REPORT.md` discusses a v*-drain super-adiabatic wall layer in the BSG (Linhao) setup, not HG; not copied.
Run settings of HG prodG itself: `../../hg_runG.cfg`, `../../hg_last_cmd.txt` (already in caltech-to-raven-1107/).
