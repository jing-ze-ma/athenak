# TASK-2026-10-09 (PROPOSAL, not pushed): Orion builds a flux-mean opacity table for the AG Car thin atmosphere

Origin: viper fluxmean-1009 worker (branch fluxmean-1009, uncommitted at the time of writing). For main to forward to Orion if the user agrees.

## Why
- The implicit M1/VET solver now has `<rad_m1>/flux_opacity = fluxmean`.
  - In thin cells (lagged flux factor f = |F|/(cE) above 0.1-0.4), the transport/force opacity becomes kappa_T = kappa_R^(1-w) kappa_fm^w.
  - Thick cells keep kappa_R bitwise.
- Today kappa_fm is the Planck mean at T_rad or T_gas plus kappa_e: an UPPER bound (saturated lines, no Sobolev).
  - On AG Car A it gives Gamma 2-20 at 1.0-2.9 R_ph, against 0.2-0.001 with Rosseland.
- The force on the outflow is uncertain by 1-2 orders of magnitude between these two bounds.
- The code already reads a THIRD table (`flux_mean_temp = table`, `problem/he_fluxmean_table`).
  - It must be on the Rosseland/Planck grid and is looked up at (rho, T_gas).
  - Gate: fed the Planck table, the run is bitwise equal to `flux_mean_temp = tgas`.
- What is missing is the physics content of that table.

## What to build
Files: `fluxmean_A_gs98_x0.36_z0.02.txt` (and `_B_`).
- Format: the existing AthenaK stellar table format.
  - Comment lines, including "# nT nD lTmin dlT lDmin dlD".
  - Then nT*nD values of log10 kappa [cm^2/g], T slowest.
  - EXACTLY the grid of `rosseland_ext2_gs98_x0.36_z0.02.txt`: 217 x 421, log T 2.6..8.0, log rho -21..0.
- Content: ABSORPTION ONLY. No electron scattering; the code adds kappa_e from the EOS (vet_scatter_kappa_e = eos).
- Planck-type mean, weighted by the run's radiation spectrum instead of B_nu(T_gas):

      kappa_F(rho, T_gas) = int kappa_nu^abs(rho, T_gas) W_nu dnu / int W_nu dnu,
      W_nu = B_nu(T_c)   (diluted blackbody: dilution cancels),

  - T_c = colour temperature of the photospheric field: T_eff = 9000 K for A, 20000 K for B.
  - kappa_nu^abs = LTE monochromatic absorption: bf + ff + H- + lines, populations at (rho, T_gas).
  - One table per state (A, B), because T_c is a constant per run.

Variants (each its own file, same format):
1. `static`: lines at full (static) strength. This is the upper bound, to compare with Planck(T_gas) on the code side.
2. `sobolev`: each line's contribution multiplied by the Sobolev escape factor (1 - exp(-tau_S))/tau_S.
   - Uses tau_S = kappa_line rho c / (nu dv/dr) with a REFERENCE velocity gradient dv/dr = v_inf/R_ph.
   - v_inf ~ 200-400 km/s for A, 300-600 km/s for B (literature for AG Car: Groh+2009, 2011).
   - This is a crude, fixed-gradient line force. A true CAK force depends on the local dv/dr and needs code support (see the scope note below).
3. Optional: the Rosseland mean recomputed from the same kappa_nu, as a consistency check against the ext2 Rosseland table.
   - It should agree within the line-list differences.
   - Report the ratio map.

## Data needed
- Line lists: Orion's RSG wind work (branch rsg-ck-1008 / agcar-opac-orion-1009) has the line lists and the Sobolev force-multiplier code.
- Continuum: the same sources as the ext2/Orion Planck stitch, so that variant 1 reproduces the Orion stitched Planck table (planck_ext2orion, md5 493f9d01) at T_c = T_gas.
  - That is the validation of the integrator.
- Composition: X 0.36, Y 0.62, Z 0.02 (GS98 scaled), as the existing tables.

## Validation and deliverables
- Validation: variant 1 with W_nu = B_nu(T_gas) must reproduce the stitched Planck table to < 5 % where both are defined.
- Deliverables:
  - the tables;
  - md5s;
  - a short README with the integration recipe;
  - the ratio maps kappa_F/kappa_R and kappa_F/kappa_P along the AG Car A/B atmosphere paths (rho 1e-16..1e-11, T 3500..20000 K).
- CPU only. Self-contained. No AthenaK changes needed: viper reads the file with `flux_mean_temp = table`.

## Scope note (c): a real line force (CAK/Sobolev), later
- A CAK force multiplier M(t) = k t^-alpha needs t = sigma_e rho v_th/|dv/dr| per cell, i.e. the local velocity gradient.
- It would enter as a separate explicit gas force (not through the transport opacity), so the implicit M-matrix is untouched. The momentum it gives the gas is then taken from the radiation flux as an explicit sink, to keep conservation.
- Orion's force-multiplier fits (k, alpha, delta) for the AG Car composition would be the input.
- Not started.
