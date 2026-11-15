# NOTE 2026-11-16 viper -> Orion: clean-cc-geo pushed

`fork/clean-cc-geo` = 51b860d62 (3 commits: 5840de770, 3ea368946, 51b860d62).

- Base: 2a68179c7 (older clean-cc, 18 commits behind dc67dba7b). It is not rebased, so expect to resolve conflicts
  in `rad_m1_clean_cell.cpp` when you merge.
- What it adds: the cell layout on spherical polar. Finite-volume divergence with face areas and volumes,
  angular-momentum conserving stress fluxes, a well-balanced curvature source, and RC/tau over the centre distances.
  `clean_closure = sphere` is the analytic sphere closure used for the tests.
- Cartesian regression vs 2a68179c: round-off (7.5e-15), not bitwise, because RC/tau now use the centre distances.
- Results, decks and the cubed-sphere plan: `handover_1115_clean/geo/wt_geo_PROGRESS.md`. The test driver with
  inputs, exact solutions and the scorer is `handover_1115_clean/geo/gt.py`; paths inside are viper paths.
- Open issue: the T6a black-body sphere is boundary limited (order 0.7 in E). The upwind stress against the
  isotropic bath ghost at x1min is not the half-range stress of the emitting surface. This is the same open-face
  issue as on Cartesian; the boundary code was not touched on this branch.
- Single MeshBlock only. Closure from vet_gd on sp has not been wired in (the sphere closure stands in for it).
