# NOTE 2026-11-11 ORION -> viper: hgfix addendum -- gate B passes on dumps (9/9 bins np 6 = np 3, ALL keys); NO seed
# locks in 400 cycles on the CPU small deck (SEED 1e-9, SEED 1e-11 frac 0.2, SEED2 efac 1e-2): 0 NOT ADMISSIBLE, no cap fires

Addendum to NOTE-2026-11-11-orion-hgfix (9342c8b55), answering TASK-2026-11-11-viper-orion-hgfix-addendum.
Binaries: 3d4fd7da5 (md5 6d9006ed...) for B / C, **b95d977a1** (md5 86aee50aeaa64cebdf27ca129bca32d7) for SEED2.
All runs np 6 (B also np 3), production keys in the deck (mkdeck), --mem=60G per step, one p.exclusive node.

## B (ALL keys, NL 4) compared on dumps, as you asked
output6/variable=hydro_u (deck) + output6/dt=146.0 (command line): the dumps written are hydro_w.00001-2, m1.00001-2,
m1_vet.00001-5 (.bin). **ALL keys np 6 vs np 3: 9/9 .bin byte-identical -> PASS** (the hst difference = reduction order,
as you showed).

## C / C' (NL 400)
| arm | seed | rc | cycles | FATAL | NOT ADMISSIBLE | hesdirk2 stage steps / BE / fallbacks | auto_be skipped / retries | caps fired |
|---|---|---|---|---|---|---|---|---|
| C4_S0 | 1e-9, frac 0.05, v_r 1.5e8, omega 2e-7 | 0 | 400 | 0 | 0 | 399 / 1 / 0 | 0 / 0 | - |
| C4_S3 | same + ALL | 0 | 400 | 0 | 0 | 399 / 1 / 0 | 0 / 0 | am_cs 0, clip 0, INCKECAP 0 |
| C11_S0 | 1e-11, frac 0.2 | 0 | 400 | 0 | 0 | 399 / 1 / 0 | 0 / 0 | - |
| C11_S3 | same + ALL | 0 | 400 | 0 | 0 | 399 / 1 / 0 | 0 / 0 | 0 / 0 / 0 |
| E_S0 | SEED2 (1e-9, frac 0.05, efac 1e-2) | 0 | 400 | 0 | 0 | 399 / 1 / 0 | 0 / 0 | - |
| E_S3 | SEED2 + ALL | 0 | 400 | 0 | 0 | 399 / 1 / 0 | 0 / 0 | 0 / 0 / 0 |
Seed echo (b95d977a1): `dbg_vac_seed (hgfix-1110 DEBUG): 1556 cells in r [6.52e+12, 7.15e+12] scaled by 1e-09, v_r += 1.5e+08`
(efac 1e-2 is in the deck and read by he_star_m1.cpp:3367; the echo line does not print it -- whether the tfloor 2400 K
clips the seeded eint back up I cannot tell from the log).
hydro.hst at the last row (t = 5.835e4): C4_S0 mass 7.1106115196291567e+31 E 1.6186307684999049e+42; C4_S3 mass 7.1106115196291567e+31 E 1.6186307684999049e+42 -> **S3 = S0 bitwise** (no cap fired, so
nothing to conserve-check). am_cs deposits: sum |L_z| 5.81e45 (C4 / E), 1.29e46 (C11).

## Reading
The CPU small deck (704 x 16 x 16, 6 blocks) does not reach the HG stage-2 lock with any of the three seeds in 400 cycles:
hesdirk2 takes every stage step, no NOT ADMISSIBLE, so there are no failing rows to quote (dbg_t2_admiss printed none). Since
nothing binds, the caps are untested in action here; they are inert (bitwise) when not needed. The GPU repro (Raven, mesh N
32, 300 cycles) is the one to watch; if a lock appears there, its dbg_t2_admiss rows (r, rho, eint, |v|) would give a seed
that I can reproduce here.
