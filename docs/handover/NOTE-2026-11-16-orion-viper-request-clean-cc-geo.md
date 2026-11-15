# NOTE 2026-11-16 Orion -> viper: please push clean-cc-geo (spherical-polar cell layout) to the fork

Orion is moving to the geometry features of the clean scheme (order set by the user: multi-MeshBlock (done on Orion),
3-D (in progress), spherical polar, cubed sphere). Your HANDOVER says clean-cc-geo (51b860d62, local on viper) holds the
spherical-polar support for the cell layout with results in progress/wt_geo_PROGRESS.md. Please push it to the fork as
fork/clean-cc-geo (no rebase needed; we will merge it ourselves) and reply with the head hash and any notes beyond
wt_geo_PROGRESS.md (open issues, decks used, reference results).

Status on Orion (local integration branch clean-cc-int, not pushed yet): cell layout with extrap1 (1 formal solution/step,
order 2), unclipped clean state, consistency (HOLO) term with ray face fluxes from the SC sweep, a kinetic upwind operator
for the deviation from the rays' moments in streaming faces (stage-1 now accepted on cb/cbg/ws66, every later BE step
realizable, 3-20x fewer BiCGStab iterations), round-off convergence exit, BiCGStab breakdown restart, multi-MeshBlock/MPI
(single-block answer to solver tolerance). Gate at production tolerance passes (space self-convergence 1.8-2.5, time 1.7-2.3).
