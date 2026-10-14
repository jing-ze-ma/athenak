# NOTE Caltech -> raven (MG): FIN4 suite — QUEUED (10-13 17:10 PDT)

- FIN4 fin4-1014 @ d505d21f, plain CPU with the ENTHBE / mbband flags (gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI,
  Release, built_in_pgens, no Athena_M1_VGD_NMAX): `athena_cpu_omp_d505d21f` md5 **bba193085d4b6590eda1c708c7f51dad**
- OLD (part B): `athena_cpu_omp_281dd686` md5 9ae9f09aeb81841a275e8d046a0eebe4
- Decks AS IS, no keys. Per-dir step = run_fin4.sh's oneA / oneB (f4_one.sh) with the Caltech fixes: rstonce.py (one rst
  block) in B, BIT verdict FAIL when any arm rc != 0, killed runs never cached. Part B first, then A (A threads =
  enthbe/nt.txt: cloud + rs n1024 at 32, ddf / ddiff at 4 in an expansion job).
- Jobs: debug 6-node rounds 4482167 -> 4482168 -> 4482169 -> 4482170 -> 4482171; expansion 4482172 (A cart/ddf + ddiff).
- Login smoke: B ss/n32 BITWISE, A sp/relax1_T/HFfZ_1841da_1 BITWISE; both logs print time2_enth_be=true
  implicit_cs2_vimp_npass=3.

## INTERIM 10-13 18:50 PDT — part B DONE (round 4482167), part A 213 / 214 (last: cart/rs_m2/FIN3_n1024 running)
- Part B (145): 112 BITWISE + 18 RAHIST + 13 DIFF + 2 FAIL — exactly as expected. DIFF = the 13 sp Marshak dirs
  (marsh_r1_T x 8, marsh_r1_u 128/256, marsh_r1_s 64/128, marsh_r1000_u 128). FAIL = the 2 known p2 n40 G2 dirs (vet_gd too
  many directions, OLD too). Every FIN4 log prints time2_enth_be=true implicit_cs2_vimp_npass=3.
- Part A so far: 213 BITWISE vs ENTHBE (incl. ddf / ddiff 4 threads in expansion job 4482172, rs_m3 n1024 at 32), rc 0.
- Rounds 4482169-71 cancelled (not needed).
