# NOTE Caltech -> raven (OD): mbband-1013 IDORT suite, two arms — QUEUED (10-13 11:45 PDT)

## Binaries (plain CPU: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on, Release, PROBLEM built_in_pgens, no Athena_M1_VGD_NMAX — same as OLD)
- OLD 281dd686: `athena_cpu_omp_281dd686` md5 9ae9f09aeb81841a275e8d046a0eebe4
- NEW mbband-1013 @ 34f79349: `athena_cpu_omp_34f79349` md5 024a5ade2a4201e670a5185577282e5f

## IMPORTANT: the raw cmprst.py on the rst says (a) DIFF in EVERY dir; that is a restart-FORMAT artifact, not physics
NEW contains 9972a56d (rotavg-rst2-1013), OLD does not. Every IDORT deck sets vet_rotate_avg = 4, so every rst carries the
M1ROTAV1 header, which has 11 int64 entries in NEW (+ i0) vs 10 in OLD; the bytes after <par_end> shift by 8.
Checked on ss/n32 with the rst kept: after dropping the 11th entry (i0 = 0 there) the header's first 10 entries and ALL data
bytes are identical; the m1 bin outputs are bitwise too.
So the Caltech runner compares with `cmprst_ra.py` (in the results tarball): BITWISE = identical once the RA header's extra
entry is dropped; RAHIST = every differing byte lies in the vet_rotate_avg history slabs (NEW stores only i >= i0); DIFF =
state differs. The raw cmprst.py verdicts are kept per dir in bit_raw.txt.

## Run
Caltech debug QoS (30 min, 1 running job, <= 10 nodes): one 6-node job per round, one worker per node, shared claim queue
(mbband first, then ENTHBE). Per-dir step = run_mbband.sh's step (mb_one.sh) with the Caltech fixes: rstonce.py (one rst
block, dt 1e30), both shockref paths (raven + viper), gtime stand-in, killed runs never cached.
Long dirs (full_b of cart/cloud/FIN3_r12g2, ddf/FIN3_d, ddiff/FIN3_d, shd/FIN3_r12g2, rs_m2/FIN3_n1024, rs_m3/FIN3_n1024)
run alone on a node with 32 threads (BIT arms of those dirs too); all others with the list threads, 8 dirs per node.
Jobs: 4454667 -> 4454668 -> 4454669 (afterany chain; more added as needed).

## INTERIM 10-13 13:05 PDT (rounds 4454667, 4456444 done; 4456533 running)
- Round-1 lessons: the BIT step with 3 x 32 threads on a 32-core node was 10x oversubscribed (cloud BIT 985 s); now each
  node has a 32-core token pool (a BIT dir reserves 3 x nt) and long dirs run BIT at 10 threads/arm, full_b at 32.
  Rounds 4454668/9 and 4455794-6 were cancelled (old scripts); chain now 4456444 -> 4456533 -> 4456535 (+ more).
- cart/ddf/FIN3_d and cart/ddiff/FIN3_d full_b need ~37 min even at 32 threads (> the 30-min debug wall): they are in the
  2-h expansion job 4456442 (4 threads), Slurm start estimate 19:10 PDT.
- bit.txt 43 / 359: arm (a) 42 BITWISE + 1 RAHIST, arm (b) 24 BITWISE + 18 DIFF (cells acting).
- The RAHIST dir: sp/dif_r1/HFfZ_1841da_128 (single MeshBlock 128x4x4). There i0 = 126, so NEW stores 1 history slab vs
  OLD 12; RA header entries other than nra (+ i0) equal; every rst byte before the history identical (except the data_size
  word, which counts nra). (cmprst_ra.py was extended for nra != nra; the first pass had called it DIFF. Not a real DIFF.)
- Raven b77ed1c1 check (bin/ and tab/ of bitO vs bitA, bin compared after <par_end>): all 43 dirs BINTAB_BITWISE, incl. the
  RAHIST one. Per-dir list: bintab_all.txt (final tarball).
- full_b: 31 rc files so far, all rc = 0.

## INTERIM 10-13 14:50 PDT (debug rounds 4461226-28 done, 4461229 running; expansion 4456442 running since 14:34)
- BIT 359 / 359: arm (a) 299 BITWISE + 58 RAHIST + 2 FAIL; arm (b) 260 BITWISE + 58 RAHIST + 39 DIFF + 2 FAIL.
  NO real (a) DIFF. All 58 RAHIST dirs (sp/dif_r1*, patm*, rlat*, nlte2_r1000_tb ...: i0 > 0 so NEW stores fewer history
  slabs) have bin/ and tab/ BITWISE OLD vs (a) (Raven b77ed1c1 check: 357 BINTAB_BITWISE + the 2 FAIL dirs).
- The 2 FAIL dirs: p2/pinc_ha20_sr30_n40/HFfZR0G2_66bdd7_64 and _d747a5_64 (GL 20 x 24 = 480 directions): every binary
  (OLD and NEW, all arms, full_b too) stops at start with "FATAL <rad_m1>/vet_gd: too many directions" (default
  M1_VGD_NMAX 432; no Athena_M1_VGD_NMAX in OLD 9ae9f09a, so not in NEW either). NOTE: CVALL's Caltech bit_all listed these 2
  as BITWISE; that was false (both arms FATAL, the compare saw the t = 0 rst). The runners now return FAIL when any arm rc != 0.
- full_b: 355 rc files, all rc = 0 except the 2 FAIL dirs (rc 1); cart/ddf, cart/ddiff in the expansion job, 2 more running.
