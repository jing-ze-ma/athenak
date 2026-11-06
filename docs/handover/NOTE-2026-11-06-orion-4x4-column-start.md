# NOTE 2026-11-06 Orion -> viper, Caltech: 4 x 4 Athena++ column job 227389 is RUNNING on Orion, so please cancel Raven's 31980425 duplicate (TASK-2026-11-06q)

Smoke test: dt 1.83317e-3 code, as on Raven. The run uses 64 ranks with `meshblock/nx1=4` (64 blocks) instead of 32 ranks and 8, because the USER allowed more cores. bsg.hst was identical at 32 and 64 ranks over 400 cycles, to its 6 printed digits, and 64 ranks is 1.7x faster, so ~30 min to 25 d if the rate holds. Radial faces are explicit, with no level boundaries, unlike `his1d_g`'s SMR; the wall (r < 35) is uniform in both.
