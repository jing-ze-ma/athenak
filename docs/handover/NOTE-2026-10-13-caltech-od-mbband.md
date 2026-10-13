# NOTE Caltech -> raven (OD): mbband-1013 IDORT suite, two arms — DONE (10-13 16:25 PDT)

## Binaries (plain CPU: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on, Release, PROBLEM built_in_pgens,
no Athena_M1_VGD_NMAX — the same as OLD)
- OLD 281dd686: `athena_cpu_omp_281dd686` md5 **9ae9f09aeb81841a275e8d046a0eebe4**
- NEW mbband-1013 @ 34f79349: `athena_cpu_omp_34f79349` md5 **024a5ade2a4201e670a5185577282e5f**

## Jobs
Caltech debug QoS, 6-node rounds (one worker per node, shared claim queue with OT ENTHBE and WV wavefix): 4454667,
4456444, 4456533, 4461226-30, 4468250-53 (cancelled, no results used: 4454668/9, 4455794-6, 4456535, 4459852-5,
4468254/5). Expansion jobs: 4456442 (cart/ddf + cart/ddiff, 4 threads), 4472288 (cart/ddn 4 threads, cart/pipe/FIN3_c220
24 threads): these full_b runs take 35-46 min, beyond the 30-min debug wall.
Threads: list threads, except long dirs: cloud/shd FIN3_r12g2 and rs_m2/rs_m3 FIN3_n1024 BIT 10 (cloud, shd: 32) per arm
and full_b 32; pipe/FIN3_c220 full_b 24 (<dir>/full_b/nt.txt; bit.txt "bitnt=").

## BIT: rst compare (cmprst_ra.py) — NO real (a) DIFF
NEW has 9972a56d (rotavg-rst2-1013), OLD not: the M1ROTAV1 header has 11 entries in NEW (+ i0), and NEW stores the
vet_rotate_avg history only for i >= i0. cmprst_ra.py: BITWISE = all data after <par_end> identical once the extra
header entry is dropped; RAHIST = the RA header agrees except nra (+ i0) and every byte before the history slabs is
identical (except the data_size word, which counts nra) = all state bitwise; DIFF = state differs. Raw cmprst.py
verdicts: <dir>/bit_raw.txt (DIFF everywhere, from the header).
- arm (a): 299 BITWISE + 58 RAHIST + 2 FAIL — no DIFF.
- arm (b): 260 BITWISE + 58 RAHIST + 39 DIFF + 2 FAIL.
- Raven b77ed1c1 check: bin/ and tab/ of bitO vs bitA (bin after <par_end>): BITWISE in all 359 dirs, incl. all 58
  RAHIST dirs (bintab_all.txt).
- (b) DIFF dirs (39): cart cloud/FIN3_r12g2, hcyl g8x32 + n512, marsh c1 + c10, oneo b3k3/5/7, pipe c220 + c2200, sdifk,
  shd r12g2 + r36g2, shdg g8x32, slabc n32-n512 (5); p2 pinc (2), pinc_ha20_sr30 (4), pinc_ha20_sr30_n40 HFfZR0 (2);
  sp core_L16_u (3), core_u (3), marsh_cart, marsh_r1_s/64, nlte2_tb_cart, pinc (2), shd/32.
- The 2 FAIL dirs: p2/pinc_ha20_sr30_n40/HFfZR0G2_66bdd7_64 and _d747a5_64: every arm, OLD too, and full_b stop at start
  with "FATAL <rad_m1>/vet_gd: too many directions" (GL 20 x 24 = 480 > the default M1_VGD_NMAX 432). (CVALL's Caltech
  bit_all called these 2 BITWISE: false, both arms FATAL and the compare saw the t = 0 rst; the runner now says FAIL.)

## FULL_b (NEW + vet_path = cells, every dir): 359 / 359 rc files; rc != 0 only the 2 FAIL dirs above (rc 1)
full_b_lines.txt: "FULL_b dir rc time nt subdivided-steps"; the vet_gd subdivided-step counter is printed in 199 runs,
non-zero in 191.

## Files (docs/handover/caltech_mbband/)
- mbband_runs_light.tgz (30 MB, md5 file) = the TASK tar + `--exclude='*.bin'` (with bin it is far over the push limit;
  kept at Caltech /resnick/groups/carnegie_poc/jingze/mbband_1013/caltech_mbband/runs), plus the tables and scripts.
- bit_all.txt, bintab_all.txt, full_b_lines.txt, cmprst_ra.py.
