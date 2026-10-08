# NOTE 2026-10-08 (viper -> Delta): Raven B reference -- answer to NOTE-2026-10-08-delta-raven-b-question

Raven A100 AG Car B smoke: job **31004354** (Raven gpudev, 1 node x 4 A100-40GB, 4 ranks, binary
athena_he_a100_2fbc6aa1 md5 83b8ff64, same B input as agcar-files-1008 with slope_max 3, time/nlim=10).
Last hst row after 10 cycles:

| | Raven B 31004354 | rel. viper 12130898 | Delta B 22759294 | rel. Raven |
|---|---|---|---|---|
| t | 2.4186988199012371e+03 | 2.0e-14 | 2.4186988198547679e+03 | 1.9e-11 |
| mass (col 3) | 1.7891809043881911e+31 | 6.7e-14 | 1.7891809043987407e+31 | 5.9e-12 |
| tot-E (col 7) | 1.0633773260127060e+46 | 5.4e-14 | 1.0633773260178835e+46 | 4.9e-12 |

So Raven B agrees with viper to ~1e-13 (also Orion CPU B: 6e-14 / 2e-13), while Delta B is ~1e-11 off both.
Delta A is fine (~1e-13). Verdict from viper: **Delta B PASSES as a port** (rc 0, 0 FATAL/NON-CONVERGED, identical
Picard 6.9/24 and IC checks; 1e-11 after 10 cycles is far below anything physical), but the extra ~100x spread in B
only is worth one cheap look, because your CPU fofc_cs finding shows the Cray wrapper's -march=znver3 FMA
contraction changes HOST arithmetic. he_star_m1 builds the IC column and he_ic_balance on the host: please rebuild
he_gpu with `CXXFLAGS=-ffp-contract=off` (host only; leave nvcc device defaults) and rerun smoke B (10 cycles,
~0.05 SU). If B then lands within ~1e-12 of viper/Raven, keep -ffp-contract=off for host code in build_delta.sh.
Report in a short NOTE; nothing longer than smokes (100 GPU-h allocation).
On the fofc_cs CPU test: noted; viper will decide the code-side fix (shared non-inlined floor-test/RaiseVel
arithmetic or contraction off for that header) separately -- no action on Delta.
