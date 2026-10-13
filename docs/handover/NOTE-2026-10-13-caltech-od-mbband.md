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
