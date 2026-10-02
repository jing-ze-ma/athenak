# NOTE 2026-10-02 ~05:00 (viper): latest state for all sites (Caltech, DeltaAI)

rt-integration on the fork = this commit. Everything below is merged and gated (keys-off bitwise, CPU suite, GPU smoke).

## BSG (Ma+2026 reproduction)
- **Keep building/running BSG from 30bf6c03** (= NOTE-2026-10-02-bsg-code-ready). Viper's fresh BSG arms run the same
  commit (binary athena_he_gpu_30bf6c03), so all sites stay comparable. Later he_star_m1 changes (below) are small and
  should NOT be switched in mid-run.
- **Input change (bitwise for BSG):** drop `<rad_m1>/implicit_eos_cache = true` from bsg3d_arm2.athinput. The gate
  showed it never hits for BSG (1 miss per cell-pass): identical results, 4.7 % faster without it.
- Picard-mean flip (2.994 vs 3.000): resolved - a lucky first linear solve at viper's step 78 (bound 1.5e-10), seen also
  CPU-vs-GPU on one machine at Caltech. Report-only in the gates.
- Viper status: fresh arm 2 (low top) and arm 1 (80-Rsun top) queued with guarded short (6 h) + long (24 h) jobs sharing
  one run dir (whichever starts first; others exit). Not started yet (apu queue ~10 h).

## New defaults since 30bf6c03 (fresh runs only; restarts whose file lacks the key keep the old behaviour)
- he_star_m1 / implicit M1 (bf3eb6db): positivity set on (implicit_g0_exchange, g0_limit 1, pos_gas, pos_floor,
  opac_newton_guard_mode 6; neutral within the 1-ulp twin spread on the BSG column); <hydro>/fofc true in he_star_m1 and
  box_convection pgens; <rad_m1>/implicit_resid_fatal = 1e2 (FATAL on a diverged solve; He runs that finished reached
  NC residuals up to ~1, so 1e-2 would be too strict).
- dhj (WASP-121b), 3f3b0414: <problem>/ck_wall_flux_exact = true - faces 5/6 injected 7.2 % too much internal flux at
  the bottom wall (4/(2+sqrt3) sigma T_int^4; effective T_int 577.5 K instead of 567.4 K). Face 6 now gives
  1.0000045 sigma T_int^4. ANY face-5/6 run started before 3f3b0414 has the 7 % excess; restarts keep it unless the key
  is named. ck_sph_face = 6 is the dhj default since d26eeda2 (fresh runs).
- WASP-121b ck cadence block (ck_impl_every 4, once, jreuse 0.2, pred): FAILS a noise gate vs the exact settings
  (OLR 50x, upper-atmosphere band T 500x the twin spread; 2.8x cheaper). Fine for the 300-rotation relaxation; use the
  exact settings for steady-state science.
- implicit_eos_cache: recommended off everywhere (neutral, never faster).
- Known pre-existing test failures: test_in_cshock2d_mpicpu, test_sbox_{hydro,mhd}shwave_mpicpu,
  test_rad_lwave2d_amr_mpicpu (segfault also at 835a06bc).

## Viper runs queued now
BSG arms 1+2 (above); He presn ext 128^2 scout + acc (fresh, 30bf6c03); WASP-121b fresh face 6 + wall fix
1x/3x/10x (athena_f6w, 3f3b0414). Reports and pages:
BSG https://claude.ai/artifact/WTWSe9cPXwivWyihbsaG3f, He https://claude.ai/artifact/HcgYDVQ3CozhxXavZppqcU.
