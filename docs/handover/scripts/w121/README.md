# WASP-121b analysis scripts (viper, 09-29/30)

Copied from the run directories under /viper/ptmp2/jinma (paths inside the scripts point there).
Results they produced: RESULTS.md in the same-named run subdirectories; summary in
docs/handover/SESSION-2026-09-30-viper.md section 3.

| dir | source dir | what |
|---|---|---|
| ana_rot300, ana_rot300_10x | w121prod_0929/ana_rot300[_10x] | rot-300 profiles: deep drift, adiabaticity, jet, day-night (1x / 10x) |
| budget_rot300, budget_rot300_10x | w121prod_0929/budget_rot300[_10x] | per-shell energy/mass ledger (budget_dt output), top-cell T |
| synth_rot300, synth_rot300_10x | w121prod_0929/synth_rot300[_10x] | synthetic phase curves, NIRSpec night flux, eps, limb RVs |
| deepmix_c32 | w121prod_0929/deepmix_c32 | deep enthalpy flux mean/eddy split (C32 baseline) |
| met3 | w121prod_0929/met3 | 3x ck table (ln k/X interp in log Z) and RCE IC |
| oddeven_0930, x1phi_0930 | w121prod_0929/{oddeven_0930,x1phi_0930} | lhllc radial odd-even mode diagnosis and fix comparison |
| c256 | w121_c256_0929 | C256 remap helpers (head/tail), dt estimate, smoke analysis |
| mhd_eta | w121_mhd_0929/deltaai_pkg | MHD max_eta scan analysis (ana_eta.py <armdirs>) |

LITERATURE.md: WASP-121b observables from the literature (obs_lit/).
