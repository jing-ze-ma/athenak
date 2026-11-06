# NOTE 2026-11-06 Caltech -> viper (USER): Linhao's actual parameter file; correction to TASK-2026-11-06c

The USER received his input file: `docs/handover/athenapp-bsg/linhao_input/athinput.bsg_star.rusty` on this branch. Please use it for the multi-D runs.

## It settles the density unit
`density_unit = 6.15245e-5` g/cc ("cgs value of density at the bottom"). So his seed threshold, rho > 0.01 code units, is **6.15e-7 g/cc**. The base (5e-5 g/cc) is well above it, so **his base IS seeded**, assuming the rho > 0.01 line in his bsg_star.cpp is the one used for this deck. The hypothesis in TASK-2026-11-06c (nothing seeded) is wrong on this evidence.

## Changes to TASK-2026-11-06c
- **SEEDED is now the priority**, with this threshold and his grid. It's the direct test of why his 3-D is quiet at 3-5 d.
- **NO SEED is optional:** run it only if it costs little once SEEDED is set up. It no longer mimics his run; it only measures how much the seed sets the timing.

## Settings to use, from the file (differences from the arms so far may matter)
- nmu 4, angle_flag 1, **nzeta 5, npsi 6** (A_af1 used nzeta 4), taucell 4, error_limit 1e-6, reduced_factor 1, cfl_rad 0.025 (we assume it's not used by the implicit solve; say if it is).
- 192 x 128 x 128 root grid with 16^3 blocks on 20-80 Rsun; SMR level 1 at 38-55 Rsun, which covers 35-55 after block rounding. Uniform in r. cfl_number 0.3, xorder 2, rad_xorder 2.
- T_unit 1.07887e6, length_unit 6.957e10, molecular_weight 0.62, GM 260.795 code, gamma 5/3, dfloor 1e-9, pfloor 1e-12, tfloor 1e-4 (code).

Please say if this file conflicts with what your arms use, or with his code.
