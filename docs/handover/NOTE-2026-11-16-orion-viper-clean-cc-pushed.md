# NOTE 2026-11-16 Orion -> viper: fork/clean-cc advanced dc67dba7b -> 9a80e1a67 (milestone)

Fast-forward push after a regression on Orion. Contents (cell layout):
- Defaults: clean_closure = extrap1 (ONE formal solution per step, order 2); clean_consistency + clean_consistency_weight = true on
  Cartesian meshes (off on curvilinear until their face sums are ready): face flux (1-w) f(U) + w [F_ray + F_kin(U - Ur)] with ray face
  fluxes from the SC sweep itself and a kinetic upwind operator for the deviation from the rays' moments (stage-1 accepted on cb/cbg/ws66,
  every later BE realizable, 3-20x fewer BiCGStab iterations).
- Clean state no longer clipped by ApplyClosureLimits (the clip leaked 1-3 % radiation momentum / created energy per BE step);
  round-off-aware Newton exit (hard-capped); BiCGStab restart on rho breakdown; time2_k1_redist off for clean.
- Multi-MeshBlock / MPI (any split incl. x1, multi-pass band fill for rays wider than a block; answer = single block to solver
  tolerance), 3-D (x3 face sums), opt-in vet_path = cells with consistency.
- Spherical polar: your clean-cc-geo merged (+ k+ columns of the sphere closure: thin T6a order 0.7 -> 1.7), vet_gd ray closure for clean
  (closure = vet_col + vet_col_lat + vet_gd), reflecting-theta fix with low-order faces; output index-check fix (m1_vsc..mhd_eta).
- Regression (production tolerance): rw space self-convergence 1.8-2.6, time 1.7-2.3, all P x tau; IDORT: cb 0.460/0.953/0.490,
  cbg 1.308/1.423/1.441, ws66 1.167/1.334/1.396, hcyl 0.0118/0.0128/0.0136, shd 0.0214/0.0152, shdg 0.0318/0.0295, nltek1 9.09e-5,
  nltek4 1.327e-3 (x0.5/x1/x2); multi-block = single block 1e-10..1e-12; 0 CLEAN WARNING.
- Known open: 3-D thin gas with feedback (blob/sphere) stalls (thin-gas momentum kick; fix next: conservative momentum deviation + steady
  radiation pre-solve); sp consistency face sums + conservative vet_gd transport + polar boundary + radial vet_gd split in the next push.
