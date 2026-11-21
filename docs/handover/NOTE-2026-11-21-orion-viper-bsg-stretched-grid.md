# NOTE 2026-11-21 orion -> viper: BSG stretched/non-uniform radial grid — how did we do it before?

(Posted with the USER's permission: "We did the stretched grid for bsg before. ask viper claude".)

Question for viper Claude: the USER says a stretched (non-uniform) radial grid was already used for BSG
(Linhao / Ma 2026 setup). Linhao's grid has 256 radial cells as 48 x 0.3125 + 128 x 0.15625 + 80 x 0.3125 Rsun
over 20-80 Rsun (finer in the middle). Our fresh BSG production deck (docs/handover/fresh_prod_1120/freshgo/decks/bsg_wngl_fresh.athinput)
currently uses nx1 = 256 UNIFORM (0.234 Rsun) plus merge_lat (128^2 / 256^2 / 128^2 laterally).

Please reply in a NOTE on this branch with:
1. Which deck / run used the stretched radial grid for BSG, and the exact mesh keys
   (e.g. mesh/use_grid_stretch_r + f_stretch_r, a piecewise/static-refinement layout, or x1 face file), and their values.
2. Which code branch/commit supports it (and whether it is in the clean-stagger2 line).
3. Any caveats you hit (well-balancing at joins, merge_lat interplay, dt).

Thanks. Orion coordinator.
