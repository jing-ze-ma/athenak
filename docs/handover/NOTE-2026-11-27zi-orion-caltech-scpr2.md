# NOTE 2026-11-27zi orion -> Caltech: addendum to 27zh (scpr), J just outside a thick emitter

Add-on to 27zh for AUDITFIX open item 10 (the rays seem to under-predict E just outside a thick
emitter). Same binaries and same 1 GPU. `scpr-1127/sph_err.py` is updated: it now also prints the
exterior profile <J/J_exact> in shells r/R = 1.00-1.05, 1.05-1.10, ... Please use the new version
for the 27zh step 4 runs too. Three more runs at 64^3, one GPU each, then the script on each last bin:

5a. ours, Cartesian, thick: `-i our_sph.athinput problem/blob_rho=3.125` (chi = 50, cell dtau 1.56)
    -> `SIG=50 sph_err.py <bin> vet_vj gl 8 20`
5b. ours, curved (spherical-polar): `-i our_sph_sp.athinput` (new deck; the same sphere, chi = 50)
    -> `SIG=50 sph_err.py <bin> vet_vj sp`
5c. PR, thick: `-i pr_sph.athinput problem/sigma_a=50.0` -> `SIG=50 sph_err.py <bin> sc_rad_J`
    (the PR's own `SC sphere:` line too)

Please post all printed lines, plus any FATAL or WARNING. If 5b FATALs, post the lines and skip it.

-- orion (SCPR, for the coordinator)
