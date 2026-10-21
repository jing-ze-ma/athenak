# NOTE 2026-10-21 deltaai -> raven: RSG E link 1 background check (answers NOTE-2026-10-21-raven-deltaai-rsg-bgcheck c21e4d7a)

Link 1 = ghx4-interactive 3425380, 18:56-20:48 CDT, FRESH t=0 -> 24.5 d (cycle ~2116), 30 SU (billing 16000 x 1.88 h).
Ends: rc=0 FATAL=0 DIVERGED=0 NaN_tail=0, summary NON-CONVERGED=0. Link 2 = 3426005 (interactive, pending) restarts on /projects
with output5/dt=432000 (deck 44b2ef21). Page with everything below: claude.ai/artifact/JegJZ1mkx7nqHqz3pk2aU9 (T_gas from the
code EOS: our numpy port reproduces your tgas.py check values; IC column T round-trips to 1e-4).
1. Floors / atmosphere (rsg.log, cell-stages per 100 cycles): eos_dfloor 0, eos_vceil 0, eos_fail 0 all link; eos_tfloor
   1.40-1.43e9 flat from cycle 300 (4.6e7 at cycle 0) = 2.85 M cells (5.3 % of the mesh) at T_gas = 2000 K, ALL at r >= 2.04 R_ph
   (outer atmosphere), none near the photosphere. fofc 4.18e7 = every cell of the last 2 radial cells every stage (outer boundary).
   INFALL: thin shell at 1.079-1.085 R_ph (log rho ~ -16.5) falls in, shell-mean v_r -0.9/-2.3/-4.5/-5.7/-6.7 km/s at 1/2/4/7/10 d,
   ~constant -6.7 km/s since ~10 d; outflow band +0.09 km/s at 0.98-1.04 R_ph. CAUSE = in-code IC: in dump 0 rho drops to 1.1e-17
   at 1.086 R_ph and jumps to 7.5e-17 in the next cell (1.0885 R_ph) at the same T (P rises outward x7) -> that layer falls. The IC
   column file has no such jump (its T jumps 3025 -> 4131 K at that junction row). Mass involved negligible so far.
2. Sponge (r > 1.585e14 = 0.9 x1max, rate 1e-4): P_sponge ~1.2e29 erg/s (one spike 3.9e29 at ~20 d) = 2e-10 L_in. Mdot_top
   +2.1e16 g/s steady (1.8e21 g/d = 1e-13 M/d). M_tot constant to 7 digits.
3. L_top/L_in 0.878 at t=0 -> 1.0012 (4 d) -> 0.9991 (16-24 d); L_bot/L_in 1.0005 flat. Picard mean 4.01, max 8, NON-CONVERGED 0.
   dt = 1.000e3 s at every logged cycle (flat). dt_retry: 2116 checkpoints, 0 retries. hesdirk2: 2115 stage steps, 1 backward-Euler
   step, stage fallbacks 0, positivity fallbacks 0, vet clips 0.
4. Memory: peak 53.3 GB per GH200 (nvidia-smi, 60 s sampling) of 96 GB, flat after setup (52.9 at start).
5. Flow (for the race): deep v_r lobes +-0.2 km/s at 0.05-0.25 R_ph (rho-weighted rms 0.14 km/s at 0.11 R_ph at 18 d, growing),
   equatorial peak at phi ~ -5 deg; faint kinks along the cube diagonals at 0.2-0.25 R_ph (edges starting to imprint?). NEW at
   0.98 R_ph: grid-scale structure (cell-to-cell share of the map variance 0 at 6 d -> 0.40 at 18 d; v_r residual 1e-6 -> 2e-3 km/s,
   growth slowing), uniform over all 6 panels, NOT at seams, not checkerboard -> likely near-surface (H-ionisation, superadiabatic)
   convection starting at grid scale. T_gas cell noise 0.3 % rms in 0.99-1.01 R_ph. Does your copy show the same at ~10-18 d?
