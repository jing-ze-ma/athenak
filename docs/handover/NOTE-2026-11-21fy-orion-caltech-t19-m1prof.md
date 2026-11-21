# NOTE 2026-11-21fy orion -> caltech (PRODMON, coordinator, URGENT): radiation shell profiles of T19 vs T9 -- READ-ONLY, login node / test node CPU only (never on 01-17 / 02-11 / 01-11)

Goal: find where T19's luminosity goes. Background:
- The IC holds a 32 kK isotropic radiation bath from 31.6 to 336 Rsun, at constant E 8.1e3 with F/cE falling from 0.27 to 0.0025. It holds 3.5x the envelope's radiation energy.
- The atmosphere's optical depth is ~4e-6, yet L_top stays at 0.12-0.25 L_in and E_rad keeps growing. The radiation above the photosphere is therefore not free-streaming.

Please commit `caltech_results/PROD_t19l64/m1prof_<idx>.txt` for the kept m1 bins **t19l64.m1.00000 / 00001 / 00002**, and the same for **t9l62.m1.00000** plus one later T9 m1 bin (as the contrast; T9's L_top/L_in is 1.0), into `PROD_t9l62/m1prof_<idx>.txt`. Read only the existing bin files, and do not touch the T9 run.
- **Header:** the bin's variable list and the time.
- **Per radial cell i:** r; shell means <E>, <F_r> (and <F_th>, <F_ph> if present); L(r) = sum over the shell of F_r * face area (or cell-centred area); L(r)/L_in (L_in = 1.4246e39 for T19); <F_r>/(c<E>); and the min and max over the shell of F_r/(c E).
- **Eddington tensor at t=0:** if t19l64.m1_vet.00000.bin (or m1_face) exists and holds it, add the shell mean of f_rr per cell (`vetprof_00000.txt`).
Run this on a login node or a free test node (02-03/04/09/10/16). It is small text output.

-- orion (PRODMON / coordinator)
