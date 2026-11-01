# NOTE 2026-11-01 viper -> Caltech (USER): integ-test-1101 = 2a6b5e9fc (your eta fix + one more save/restore); please re-run the (b) E_B check and a continuation pair

- **2a6b5e9fc** = e1ecd2101 + cherry-pick of your d0d1a738 (1708c7a8a: the probe restores eta_b and ad_coef) + 2a6b5e9fc: ProbeT also saves and restores MHD::xe_init_pending.
  - Without that, a probe that runs before the first live conversion would use up the pending x_e Saha start on the probe state. When the start is already done it is a no-op.
- **Please run on 2a6b5e9fc:**
  - (b-recheck) the fresh no-key r1m deck, 3 cycles: E_B at cycle 1 should be ~1.8e31;
  - (d) one restart-continuation pair with the defaults on (roche_wb_tan on): 150 + 150 vs 300 cycles, cycle-matched, rst payload cmp.
- **Viper:** Raven gates (a2/b2/c2/d2 vs d3f878be) and the MHD suite on 2a6b5e9f are queued, and the viper GPU build is running.
- We push 2a6b5e9fc to fin51-r-fix when both sides pass.
