# NOTE Caltech -> raven (OT): open-ot ENTHBE — QUEUED (10-13 11:45 PDT)

## Binary
- open-ot-1013-ht @ c685e6dd, plain CPU (same flags as NPALL: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI, Release,
  built_in_pgens): `athena_cpu_omp_c685e6dd` md5 70199d65354ada08349aaea2db2c4400

## Runner fix (please check)
The bundle's run_enthbe.sh sets ONLY rad_m1/implicit_cs2_vimp_npass = 3 and writes to <dir>/full. The Caltech copy
(en_one.sh) sets all three keys the TASK names: rad_m1/time2_enth_be = true, rad_m1/implicit_cs2_vimp_npass = 3,
rad_m1/implicit_cs2_vimp = true, and writes to <dir>/enthbe. Otherwise the same (both shockref paths, two-stage A then -r).

## Run
Caltech debug QoS, 6-node rounds shared with OD mbband (mbband has queue priority). The 3 long runs (cart/cloud/FIN3_r12g2,
ddf/FIN3_d, ddiff/FIN3_d) and rs_m2/rs_m3 FIN3_n1024 each run alone on a node with 32 OpenMP threads, first in round 1
(thread count differs from Raven's 4; nt.txt per dir). Jobs: 4454667 -> 4454668 -> 4454669 (more added as needed).
