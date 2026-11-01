# NOTE 2026-11-01 Caltech -> viper (USER): status: integ gates queued, hydro production binary ready, sponge verdict, vet_gd twin A/B prepared

Reply to 44cf9fd1 (final integ commit e1ecd210).

## 1. Integration gates
- **e32a104e:** 97/131 sources built (~10 min). Gates are armed to start when the binary lands:
  - (a) key-free bitwise vs 8cf6db62 (binary ready) from keep/off_0p50.rst, ~150 cycles;
  - (b) fresh no-key r1m, 0.1 rot, plus a 100-cycle no-key vs explicit-keys bitwise pair.
- **Then e1ecd210:** incremental build; (a2) key-free bitwise vs our d3f878be binary; a 50-cycle fresh no-key start to confirm the defaults block shows relax_drag = false.
- **One NOTE covers all of it.** On PASS, e1ecd210 (pushed to fin51-r-fix) is our r1m launch build.

## 2. Hydro production binary ready (fin51-r-fix d3f878be)
- athena_gpu_dhj_d3f878be_fmadoff, md5 1b99f14b. The CPU binary from the same head is bitwise vs the GPU-gated relaxdrag build.
- **Final hydro setup (USER):** top_bc outflow1 for all five planets; roche_wb_tan true; relax_drag false; top sponge ON (radial); cs_seam_flux shared; am_hst_split false; BVALS_CC_SPLIT=1.
- Fresh starts from t=0, 1 H200 per planet, to 300 rot. Launch waits for the USER's go.

## 3. Top sponge: keep it (provisional, W121 to 5 rot; W189 and the 10-rot W121 running)
- Sponge-off is identical to sponge-on below 1e-3 bar: differences at 1 bar and deeper are 20-200x below the 1-rot variability. Waves reflecting from the top do not reach p > 1e-3 bar measurably.
- **But sponge-off has:**
  - dt dips (to 5.1 s in rot 0-1);
  - 8x density-floor events;
  - top-cell v_r up to 9.8 vs 2.6 km/s;
  - 25x Mdot_top (still negligible);
  - mean dt even 0.7 % lower.
- Nothing is gained by dropping it.

## 4. vet_gd twin on HG/RSG: A/B prepared, not yet run
- From the accretor work: on cs (cs48, HG-like config geo 2 + rotate_avg 8, twin on vs off), the pole panel is identical to 4 digits and L differs by 0.15 %. Twin-off is closer to the 1-D p-ray reference in the thin region.
- **Accretor recommendation: vet_gd_twin = false** (it removes the tip bug at its root). vet_gd_twin_srw remains a bitwise-safe fix for runs that keep the twin. The cs96 mass check at t 0.015 is still running.
- **HG/RSG A/B is prepared** (twinab_1101): 4 arms, one node each, HG 3 sim d and RSG 25 sim d.
  - The OFF arm must also set vet_gd_twin_fuse = false and vet_gd_twin_state = mean (both FATAL without the twin), and vet_source_noesrc = false (the code says noesrc is only meant to be used with the twin).
  - It waits for free nodes after the hydro launch.
- **Observation:** HG production WITH the twin already shows a cubed-sphere imprint. In the 3-day window-mean T_eff, panel vertices are 4.9 % hotter than panel centres, R_ph is 6 % lower at the vertices, and 4.8 % of the T_eff variance repeats panel to panel (random ~2 %). So the twin is not removing it.
