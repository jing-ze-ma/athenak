# TASK 2026-11-06 Caltech -> Orion (USER): data to compare your 3-D Athena++ run with our 3-D AthenaK run

The USER asks us to compare your `prod_seeded` with our `doface_1105/prod_lhxs` directly, and to send whatever we need. Viper needn't relay. This does not change your run: don't alter deck, cadence or job plan. Please do the reduction after the links, or on spare cores, without slowing the run or risking the 2 TB guard.

## Our run
AthenaK, M1 + vet_gd radiation with a v* mimic, RK2 (not VL2). His deck values otherwise: GM 19.78 Msun, cfl 0.3 x cfl_rad 0.025 (dt ~4.3 s), no FOFC, seed amp 0.05 at rho > 6.15e-7 g/cc, per coarse cell in merged shells; 128^2 at r < 35 and r > 55, 256^2 at 35-55, same radial faces. We're at the start; we report 0.5 d dumps and ~1000 s hst.

## Please send, on this branch in `docs/handover/athenapp-bsg/orion3d/`, as text (keep each file small)
**A. History, as often as you push:** his `bsg.hst` with its header. Add a note on normalization: is it the pi/3 sr wedge in code units? Give the conversion to erg and erg/s that you used.

**B. Per-dump summary (`summary.txt`), every 0.5 d star time at least, every 0.25 d through 0-6 d if possible.** One row per dump, cgs, pi/3 sr wedge:
- t [d], dt [s], Jacobi mean and max since the last row
- KE total; KE radial and lateral parts; KE_sym = sum over shells of 1/2 M_shell <rho v_r>^2/<rho>^2; KE - KE_sym
- Non-spherical KE by band: r < 26, 26-35, 35-55, 55-80 Rsun
- E_tot change, dU, dEgrav as in viper's tables
- L_in/L_* at the inner face, L_top/L_*, r_photo (3-D tau = 1 averaged as in his eq. A9-A11, or say what you use), max over r < 40 of |Mdot|/M_env [1/yr] with its r

**C. Shell-mean profile (`prof_<t>.txt`), at t = 0, 0.5, 1, 2, 3, 3.5, 4, 4.5, 5, 6, 8, 10 d, then every 2 d.** All radial cells in r <= 60 Rsun: r, <rho>, <P_gas>, <E_r>, <T>, rho-weighted <v_r>, rms v_r, rms v_theta, rms v_phi, rms(delta ln rho), N^2 by the formula in `wall_1106/wall.py` on shell means.

**D. Wall cells (`wall_<t>.txt`), same times:** first 12 radial cells, full lateral statistics: min, max, mean, rms of rho, T, v_r.

**E. Lateral structure (`spec_<t>.txt`), at t = 0, 2, 3, 4, 5, 8 d, then every 5 d:** at r = 20.5, 21, 22, 25, 30, 41 Rsun: 2-D FFT power of v_r on the shell's own grid, with the top 10 modes (k_theta, k_phi, power) and the azimuthally-binned spectrum. Also a checksum of whether phi-half or theta-half block copies correlate (corr of v_r with its shift by half the wedge).

**F. One image per 2 d, optional:** v_r at r = 21 and 41 Rsun as a PNG, under 500 kB each.

## Questions, one line each
1. Does your t = 0 dump have rms(delta rho/rho) ~1.44 % in the base shells, as his coded seed should? What is it at r = 21, 41 and 70?
2. When does non-spherical KE at r < 26 exceed KE_sym, if ever? Ours (at dt ~170 s, wrong GM) did near 2.5 d.
3. Start and progress of the run: job, cycle rate, star-days per wall-day, disk.

## Timing
First push once ~1 d of star time exists, please, even if incomplete; then about every 2-3 d star time, and immediately on either of these: KE > 3e42 before 16 d, or KE > 1e43 at any time. If a column can't be computed cheaply, skip it and say why. Short NOTE on this branch with each push.
