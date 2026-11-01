# NOTE 2026-11-01 viper -> Caltech (USER, IMPORTANT before the hydro launch): roche_wb_tan breaks restart continuation in the BASE code (d3f878be)

## Bisection on fin51-r-fix d3f878be (Raven CPU, C8x112 r1m-type deck, ohmic eos + STS, 150 + 150 vs 300 cycles, rst data)
- Arms WITH roche_wb_tan = true DIFFER: all keys; no rkg_lean; no rkg_dt_safety; no top_bc; no cs_seam_fc_nolimit.
- Arms WITHOUT it are IDENTICAL: no keys; all keys except roche_wb_tan.
- So **roche_wb_tan alone makes a restart non-continuous in base code.** Every chain link of an MHD production with roche_wb_tan + resistivity is affected.
- The likely mechanism is the probe re-running at the restart and clobbering eta_b/ad_coef (your d0d1a738 fix). Our gate on 2a6b5e9f (which includes d0d1a738) will show whether that is all of it.

## Please check BEFORE launching the 5 hydro productions
Hydro has no eta, so the eta clobber cannot be the cause there. But if the probe re-runs at every restart, g_t is re-measured on the restart state instead of being carried over. Then hydro chain links are not continuous either, and the tide forcing changes slightly at every link.
- **Test:** d3f878be, one hydro planet deck with roche_wb_tan = true, from t=0: 150 + 150 vs 300 cycles, cycle-matched, rst payload cmp.
  - IDENTICAL: the hydro launch is fine as is.
  - DIFFER: the probe result (g_t, or what it froze) must be written to and read from the restart instead of re-probed. Fix that before the launch; viper can help.
- Your continuation pair on the integ commit 2a6b5e9fc (MHD, all defaults) would close the MHD side.

## Other gates
- Viper GPU 150-cycle key-free gate on e1ecd210 vs d3f878be: rst data + mhd.hst IDENTICAL (floor columns at atomics noise, Efloor_rt 1.65e-12).
- The same gate on 2a6b5e9fc is running.
