# NOTE 2026-11-01 viper -> Caltech (USER): thanks; we cherry-pick wbtan-etafix d0d1a738 into integ-test-1101. It probably also explains our restart-continuation failure.

Reply to 0e39821b and cfb2b9807. Very good catch.
- **Plan:** we cherry-pick d0d1a738 into the integration branch and NOTE the new integ-test-1101 hash. You do not need to merge wbtan-etafix-1101 into fin51-r-fix separately: when integ goes to fin51-r-fix it carries the fix. (If you merge it first anyway, the integration merge will resolve trivially.)
- **Probably the same bug:** our Raven gate on e32a104e failed the restart-continuation test (150 + 150 vs 300 cycles, all defaults on: rst data DIFFER). The tracers + AD part alone passed it.
  - If the wb_tan probe also runs at a restart, it clobbers eta for the first cycle after every restart. That would make every r1m chain link non-continuous.
  - Our bisection (rkg_lean / dt_safety / outflow1 / roche_wb_tan / cs_seam_fc_nolimit, each off) on d3f878be is running. We re-run the continuation test on the fixed commit first.
- **Please re-run on the fixed integ commit when we NOTE it:**
  - your 3-cycle (b) E_B check (expect ~1.8e31 at cycle 1);
  - if possible, one restart-continuation pair on the r1m deck with roche_wb_tan on: 150 + 150 vs 300 cycles, cycle-matched, rst payload cmp.
- Still holding the push until these pass.
