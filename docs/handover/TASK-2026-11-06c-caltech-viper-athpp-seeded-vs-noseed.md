# TASK 2026-11-06 (3rd) Caltech -> viper (USER): multi-D real Athena++, seeded vs no seed

You offered a seeded multi-D Athena++ run with this wall (NOTE-2026-11-06b). The USER asks for that plus the same run with no seed.

## Why
Our 3-D AthenaK v* arms have a linear, seed-dependent, non-spherical instability at the inner wall. Amplitude growth is ~1.8/d, KE jumps 64x over 2.5-4 d and saturates near 1e43 erg, about 8e42 above Linhao's level.
It then shuts off the wall drain, and we get no collapse. His Fig. 2 stays quiet to ~3 d and has coherent radial bands at 3-10 d, then a sudden inflow at 17-19 d: it behaves like our 1-D.

**Hypothesis (USER):** his `density_unit` puts his seed threshold (rho > 0.01 code units) above the base density (5e-5 g/cc), i.e. density_unit >= ~5e-3 g/cc, so nothing is seeded. Growth from round-off at 1.8/d would add ~15-18 d, which is close to his collapse time. This is our estimate only.
We are running the same test on Caltech (`doface_1105/prod_noseed`, kick-fixed, rms of lateral rho exactly 0 at t = 0).

## Arms (real Athena++, his code + BC, same wall as A_af1 as far as multi-D allows)
1. **SEEDED:** rho *= 1 + 0.05 (U - 0.5) wherever his condition holds, on the 3-D wedge or nearest affordable multi-D setup.
2. **NO SEED:** identical, but seed off (threshold above the base density). Nothing else changes.

Use his grid near the wall if affordable (48 x 0.3125 / 128 x 0.15625 / 80 x 0.3125 Rsun, 128^2 at 20-35 and 55-80, 256^2 at 35-55; CFL 0.3). If not affordable, use a reduced lateral wedge or 2-D, but say so: a coarser lateral grid changes the column scale, and in our runs the columns are ~5 lateral cells wide. State the angular setup (nmu, angle_flag).

Priority is reaching at least 5 d for arm 1 and arm 2 together, because that is enough to answer whether they differ. Continue to ~20-25 d if affordable, at least for the arm that stays quiet: the collapse is the real target.

## Report per arm, against time (0.5 d steps at least; hst cadence if available)
1. Total KE and non-spherical KE (KE_sym = 1/2 M <rho v_r>^2/<rho>^2 per shell, summed; report KE - KE_sym), integration range and normalization, same pi/3 sr as before.
2. rms lateral rho at the wall shells and v_r rms per shell: where does non-spherical KE sit (r < 26, 26-35, 35-55 Rsun)?
3. Shell-mean profile in the first 24 cells at 0.5, 1, 2, 3, 4, 5 d, same columns as `wall_1106/wall_profiles.txt`.
4. L_wall/L and r_photo against time; Mdot/M_env max for r < 40 if easy.
5. Any FATAL, NaN or time-step collapse, and cost.

## What we'd conclude
- Seeded breaks up by ~3-5 d, no-seed stays quiet: his published run was effectively unseeded, and the hypothesis holds.
- Both quiet: the seed isn't the answer; the difference is numerical (VL2, no FOFC, DO stencil, wall treatment).
- Both break up: he has something neither arm has, or his noise level is high enough regardless of seed. Say how early it is.

Please say first what you can afford, with cost and start time, as a short NOTE on this branch, then give results as they come. Don't wait for the end to report the first 5 d.
