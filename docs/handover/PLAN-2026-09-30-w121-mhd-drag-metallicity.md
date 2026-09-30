# PLAN (user 09-30): WASP-121b magnetic drag at 3x and 10x metallicity (MHD drag x metallicity)

## Why
- Our no-drag hydro arms: 1x -> offset 17-18 deg E, eps 0.35; 10x -> 6 deg, eps 0.13; observed ~3-5 deg, eps 0.25
  (synth_rot300*, viper). Limb Fe RVs barely change with metallicity and stay ~2x too fast on the evening limb.
- Literature (obs_lit/METALLICITY_VS_DRAG.md): every published WASP-121b GCM-vs-data comparison is 1x solar and varies
  only drag (Wardenier+2024, Frazier+2026, Davenport+2025), while retrievals give ~10x (Pelletier+2026, Evans-Soma+2025).
  A matched metallicity x drag set is new. Discriminators: Doppler winds (Seidel+2025/2026, Vaulato+2025 dKp),
  offset vs wavelength (Mikal-Evans+2023, Frazier+2026, Yang+2026).

## Arms (6): {3x, 10x} x {0 G hydro control, 3 G, 10 G}  (field = polar field at the inner wall, problem/bbot_gauss)
- Start states: 10x = viper w10x rst dhj.00600 (rot 300, exists). 3x = rot-300 state of the 3x arm chosen after the
  ck solver A/B (Caltech 3640302 with 3x keys, or viper w3xk with dtmax 0.25 + maxit 24) - decide when both exist.
- Route: identity remap hydro -> MHD (dhj_remap.py, as the 1x MHD package; HOWTO /viper/ptmp2/jinma/remap_0929/HOWTO.md),
  fresh start -i; the 0 G control = the same hydro restart continued for the same time (separates drift from drag).
- Keys: lhlld; ohmic_resistivity = eos; max_eta = the cap chosen by the max_eta scan (DeltaAI 12 arms + viper held arms;
  offline estimate: 1e13 no spurious drag below 1e-4 bar at 3 G, one decade stricter at 10 G), use_rkg_sts if the
  cap makes Ohmic dt > ~2x the CFL dt; dfloor 1e-13 (user 09-29; affects only p < 3e-8 bar); hlld_bx_zero_tol original;
  ck_impl_conserve = 1; metallicity keys and ck tables as the respective hydro arm; bbot_gauss (not code units).
- Binary: rt-integration >= b6a6eff3 (ck-conserve, units/bbot_gauss, Newton guard, c2p_track CUDA, seam MPI leak, FOFC
  inline fix). On CUDA machines this is mandatory for MHD (c2p_track).
- Duration: >= 10-20 rot per arm (09-25 finding: MHD relaxes from a hydro spin-up within ~6-10 rot at 3 G code units);
  check saturation of KE, jet speed, Ohmic heating before scoring.

## Cost (C32 x nx1 76/74/76, 1 node x 2 GPUs; from the 1x MHD numbers)
- dt ~7 s (3 G) / ~3 s (10 G) with dfloor 1e-13, or the Ohmic dt of the chosen cap if smaller (1e13: 1.8 s; STS lifts it).
- ~21-30 ms/cycle on 2 MI300A (smoke), H200 ~1.6x faster. 10 rot at dt 3 s ~ 370k cycles ~ 2.5-3 h (MI300A);
  at dt 1.8 s ~ 5 h. The 6 arms fit in single 24 h jobs on viper apu or Caltech.

## Scoring (same pipeline as synth_rot300 / synth_rot300_10x + deepmix)
- Phase curves NRS1/NRS2/SOSS/bolometric: day/night Fp/F*, peak offset, and offset vs wavelength/pressure.
- T_day, T_night, eps (Splinter+2025 0.246); Bond albedo is imposed (0.277).
- Limb LOS velocities Fe (Seidel+2025: -4.12 morning, -6.90 evening km/s), Na/H-alpha, dKp (Vaulato+2025 -15 +- 3 km/s).
- Magnetic diagnostics: Ohmic heating / L, drag time eta/v_A^2 vs advection per layer, jet u(lat,p).
- Deep: deepmix.py flux split and deep drift (as the C32 baseline, deepmix_c32/).

## Dependencies / decisions before launch
1. max_eta (and STS) from the scan.  2. Which 3x state (ck solver A/B).  3. Machine (viper apu vs Caltech queues).
4. User GO for the 6 arms.
