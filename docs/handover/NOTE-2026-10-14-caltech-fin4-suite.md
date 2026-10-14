# NOTE Caltech -> raven (MG): FIN4 suite — DONE (10-13 19:05 PDT)

## Binaries
- FIN4 fin4-1014 @ d505d21f, plain CPU with the ENTHBE / mbband flags (gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI,
  Release, built_in_pgens, no Athena_M1_VGD_NMAX): `athena_cpu_omp_d505d21f` md5 **bba193085d4b6590eda1c708c7f51dad**
- OLD (part B): `athena_cpu_omp_281dd686` md5 9ae9f09aeb81841a275e8d046a0eebe4

## Run
Decks AS IS (no keys). Per-dir step = run_fin4.sh's oneA / oneB (Caltech f4_one.sh) with: rstonce.py (one rst block) in
B, BIT verdict FAIL when any arm rc != 0, killed runs never cached. Part B first, then A; A threads = enthbe/nt.txt.
Jobs: debug 6-node rounds 4482167 (B all + A 212), 4482168 (A cart/rs_m2/FIN3_n1024); 4482169-71 cancelled (not needed);
expansion 4482172 (A cart/ddf + cart/ddiff, 4 threads, 2053 s / 2049 s).
Every FIN4 log (A and B) prints time2_enth_be=true implicit_cs2_vimp_npass=3.

## Part A (fin4_cmpA.txt, 214): 214 BITWISE vs ENTHBE (bin after <par_end>, tab, hst). No DIFF, no MISSING; all rc = 0.
(cloud / rs_m2 / rs_m3 at 32 threads as in ENTHBE: 1156 s rs_m2, 865 s rs_m3.)

## Part B (fin4_bitB.txt, 145): 112 BITWISE + 18 RAHIST + 13 DIFF + 2 FAIL — as expected
DIFF (the 13 sp Marshak dirs):
  sp/marsh_r1000_u/HFfZ_1841da_128; sp/marsh_r1_s/HFfZ_1841da_{64,128}; sp/marsh_r1_u/HFfZ_1841da_{128,256};
  sp/marsh_r1_T/HFfZ_1841da_{1,2,4,8}; sp/marsh_r1_T/HFfZR0_72ca83_{1,2,4,8}
FAIL (rc 1 in OLD and FIN4, "vet_gd: too many directions", 480 > 432):
  p2/pinc_ha20_sr30_n40/HFfZR0G2_66bdd7_64; p2/pinc_ha20_sr30_n40/HFfZR0G2_d747a5_64
rc != 0 lines: only these 2 FAIL dirs (both arms).

## Files (docs/handover/caltech_fin4/): fin4_cmpA.txt, fin4_bitB.txt, fin4_light.tgz (25 MB, the TASK tar command,
md5 file). The *.bin stay at Caltech (<EB>/runs/<dir>/fin4).
