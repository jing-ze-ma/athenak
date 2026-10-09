# NOTE 2026-10-08 (Delta -> viper): Raven reference for AG Car B smoke?

From the NCSA Delta session (TASK-2026-10-08-delta-bringup). Delta has no push credentials yet; the user relays.

Delta smoke, 1 node x 4 A100-40GB, binary rt-integration 081dd60b (CUDA 12.9, gcc 14, cray-mpich 9.1.0),
`time/nlim=10`, same inputs as `agcar-files-1008` (MD5SUMS OK), job 22759294:

| | A rel. to viper 12131991 | A rel. to Raven 31006215 | B rel. to viper 12130898 |
|---|---|---|---|
| t | 1.3e-13 | 1.4e-13 | 1.9e-11 |
| mass (col 3) | 4.9e-13 | 5.5e-13 | 5.8e-12 |
| tot-E (col 7) | 4.2e-13 | – | 4.8e-12 |

A is within the ~1e-12 pass band; B is ~10x outside it (all else matches: rc 0, 0 FATAL, 0 NON-CONVERGED,
Picard 6.9/24 as viper, IC T-check 5.9e-14, he_ic_balance 1.07898e-07).
Delta B last hst row: t=2.4186988198547679e+03 mass=1.7891809043987407e+31 tot-E=1.0633773260178835e+46.

**Question:** is there a Raven A100 run of case B (10 cycles), and what are its t / mass / tot-E after 10 cycles
(or its `.hydro.hst`)? If Raven B differs from viper by a similar ~1e-11, the Delta B result is the usual
cross-GPU spread for B and passes; if Raven B is within 1e-12 of viper, Delta B needs a closer look.
Delta files: `/work/nvme/bivj/jma20/delta_1008/smoke/smokeB.22759294/` (run.log, agcar3d.hydro.hst).
