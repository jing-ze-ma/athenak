# NOTE 2026-10-21 raven -> deltaai: RSG E link 1 — background check (check-before-prod rule)

At the end of link 1 (or after the first ~2 h interactive link), report from the hst / first dumps:
1. Floor / atmosphere: fraction of cells at dfloor / tfloor / pfloor vs time; max |v_r| and mass flux in the atmosphere above R_ph
   (any infall onto the photosphere?); 2. sponge: radial sponge active (r0 = 0.9 x1max), energy/mass it removes per cycle;
3. L_top / L_core, Picard mean/max + NON-CONVERGED summary, dt history, dt_retry count, BE-fallback (stage fallbacks) count;
4. peak GB/GH200 after the first dumps (52.9 at start).
Stop rules unchanged (summary NON-CONVERGED=[1-9]/FATAL/DIVERGED/NaN/rc).
