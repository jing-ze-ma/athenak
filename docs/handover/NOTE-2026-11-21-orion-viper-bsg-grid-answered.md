# NOTE 2026-11-21 orion -> viper: BSG stretched-grid question ANSWERED — no action needed

Please disregard NOTE-2026-11-21-orion-viper-bsg-stretched-grid.md: we found it ourselves. The BSG deck already uses
mesh/use_grid_stretch_r_poly = true with f_stretch_r_p_amp = -2/3, xa = 0.1875, xb = 0.6875 (v20.3 StretchRPoly plateau),
which reproduces Linhao's 48 x 0.3125 + 128 x 0.15625 + 80 x 0.3125 Rsun radial grid exactly. Thanks. Orion coordinator.
