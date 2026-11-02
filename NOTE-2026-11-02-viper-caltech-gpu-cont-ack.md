# NOTE 2026-11-02 viper -> Caltech (USER): push held; we diagnose the GPU continuation on MI300A too

Reply to ac096095b.
- **Held:** the fin51-r-fix push stays held until GPU continuation passes on at least one machine.
- **Viper:** our GPU fresh-r1m continuation (20 + 20 vs 40 cycles, ohmic eos + STS, all defaults) on af982e409 is running on MI300A. If it differs, we bisect on GPU: ambipolar none, xe_tracer false, floor_tracer false, top_bc hse, rkg_lean false, one at a time.
- **ad_coef: a likely culprit in our code.** The AD EMF averages eta_A over neighbouring cells (incl. ghosts and the is-1 / ie+1 layers). If SetAmbipolar fills only part of the range, or ad_coef is allocated without initialization, the top or ghost values are garbage on GPU.
  - Also, rho_i -> 0 in top/ghost cells (x_e scalar ~0) makes ad_coef = p/(gamma_in rho_i rho_n) huge, and coef * B^2 could give inf*0 = NaN where B ~ 0.
  - We are checking the set/read ranges and adding a safe rho_i floor or full-range initialisation as needed.
- **Your initcheck/memcheck run will settle it fastest.** Please send the kernel/array name if it flags one.
- **Your g_t-in-rst robustness change (4029b621):** fine with us as an addition. It avoids the probe at restarts entirely.
