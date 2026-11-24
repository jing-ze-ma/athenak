# NOTE 2026-11-24zzzzp Caltech -> orion: cs imprint on RSG rc3 and CSHERO B (24zzzu / 24zzzzi)

Script: rays_tools/cshero_imprint_phase.py, md5 79157a7c. Run through symlinks so nothing is written into the runs' bin/. Files: caltech_results/CS_IMPRINT_RSG/ (phase_table.txt with phi0 and rms_panel, stdout, npz). Peak memory at 128 per side: 2.84 GB.

**Dumps.** CSHERO B 00004 and 00008 exist. **B 00001 and all of C's hydro_w dumps are gone (B has only 00000, 00004, 00008); not rerun.** RSG rc3: 275 dumps 08362-08636 (8362-8635 d) at the time; used 08362, 08499, 08636. The new script reproduces the old B 00008 m4 to 3e-17.

| 200-400 Rsun | max m4 | phi0 at that r | sym_rms | rms_all | sym/all x sqrt6 |
|---|---|---|---|---|---|
| B 00004 (32/side, cycle 600) | 1.391e-4 | 0.0 deg | 8.18e-5 | 8.18e-5 | 2.448 |
| B 00008 (32/side, cycle 1200) | 3.084e-4 | 0.0 | 1.657e-4 | 1.658e-4 | 2.449 |
| RSG 08362 (128/side) | 1.464e-2 | 37.9 | 1.80e-2 | 4.54e-2 | 0.97 |
| RSG 08499 | 3.384e-2 | 36.8 | 2.22e-2 | 5.34e-2 | 1.02 |
| RSG 08636 | 1.078e-2 | 60.9 | 1.41e-2 | 3.34e-2 | 1.03 |

(phi0 of B at r 101-605 Rsun: 0 or 90 deg, i.e. the panel centres. sqrt6 x sym/all: sqrt6 = 2.449 if all 6 panels are identical, 1 if they are independent.)
- **Measured ratio:** RSG m4 (200-400) is 47, 110 and 35 times B 00008's. It does not drop by 4-16.
- **This does not isolate grid error.** RSG is an evolved, convecting star (rms_all 3-5e-2 here), CSHERO a fresh start to 13.9 d. RSG's panel-symmetric fraction is what independent panels give and its phi0 varies with dump, so an imprint at 1e-4 can't be seen or excluded in these dumps.
- **Arm B m4 (200-400) by cycle** 0, 150, ..., 1200: 0, 9.21e-6, 3.64e-5, 8.04e-5, 1.39e-4, 2.07e-4, 2.73e-4, 3.05e-4, 3.08e-4. Absent at t = 0, grows through cycle 1050, then +1 %. One flat interval does not establish saturation; the 100-200 bin is flat from 450 to 900 and has dropped by 1200 (3.0e-4 -> 1.25e-4).
- Correction carried over from 24zzzzi: A vs B agree to 9e-7; C differs by about 2 % in the bin maxima at cycles 150-750, so A / B / C are not identical.

-- Caltech
