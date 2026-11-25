# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzt caltech -> orion: HG gate for rad_vet/time2_gtime_slope (your 24zzzr)

- rst hegiant.00205 md5 82058331... confirmed by the coordinator; Caltech copies of the 3 tables (md5 identical). Binary athena_sdorder_a840d7e92 md5 1e6d0b3dc87c5b0bca468d8ee9159d07, tiny gate PASS 3.686e-11.

## +1 step check (nlim 794785): ok=yes
- 2026-11-24 14:10:46 rc=0 (srun) ended: Terminating on cycle limit wall=57 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 137487; CW 0, FATAL 0; 

## Panel-edge |F| > cE (24zzzzc / zzzzg): +6-cycle vet dump, H0 and HT
```
R = 4.298e12 cm (61.78 Rsun, rst header)
== H0 *vet*.bin cycle ?
== HT *vet*.bin cycle ?
ls: cannot access '/resnick/groups/carnegie_poc/jingze/hg_gate_1124/H0/bin/*.bin': No such file or directory
ls: cannot access '/resnick/groups/carnegie_poc/jingze/hg_gate_1124/HT/bin/*.bin': No such file or directory
```
- dumps absent; only the +6 one counts; they are deleted after this.
**POST CHECK FAILED: fewer than 4 seam_vet rows**

## H0 / H1 / HT (11-24 14:12)
- H0: 2026-11-24 14:12:13 rc=134 (srun) ended: Aborted wall=52 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 140457; last log line: `srun: error: hpc-sm-02-09: task 1: Aborted (core dumped)`
- H1: 2026-11-24 14:12:16 rc=134 (srun) ended: Aborted wall=55 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 140869; last log line: `srun: error: hpc-sm-02-03: task 3: Aborted (core dumped)`
- HT: 2026-11-24 14:12:18 rc=134 (srun) ended: Aborted wall=57 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 136743; last log line: `srun: error: hpc-sm-01-17: task 3: Aborted (core dumped)`
| arm | rows | dt min/median(15-51)/max | s/cyc | s per sim s | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/lines | SDO lines, K1K2<0 mean share, max ratio, cells >0.1, >1, ef<=0 | L_top/L_in start..end |
|---|---|---|---|---|---|---|---|---|
| H0 | 1 | 90.00 / 90.00 / 90.00 | nan | nan | 0 / 0 / 0 / 0 | nan / nan / nan / 0 | 1, 0.3464, 2.366e-02, 0, 0, 0 | 0.9628 .. 0.9628 |
| H1 | 1 | 90.00 / 90.00 / 90.00 | nan | nan | 0 / 0 / 0 / 0 | nan / nan / nan / 0 | 1, 0.3464, 6.109e-01, 9433, 0, 0 | 0.9628 .. 0.9628 |
| HT | 1 | 90.00 / 90.00 / 90.00 | nan | nan | 0 / 0 / 0 / 0 | nan / nan / nan / 0 | 1, 0.3478, 3.419e-02, 0, 0, 0 | 0.9632 .. 0.9632 |
no pair H0 H1
no pair H0 HT
**POST CHECK FAILED: no verdict**

## 4 ranks do not fit: H0, H1, HT aborted at cycle 2 (11-24 14:13); the table and seam section just above are from these aborted arms and count for nothing
- srun exit 134 on all three; log: `  what():  Kokkos ERROR: Cuda memory space failed to allocate 4.254 GiB (label="cln_rk").` (peak before it 137-141 GB of ~141 per GPU). The +1 step fits (137.5 GB), cycle 2 does not. FAIL of memory, not of the method.
- Same rst, keys and binary now on **8 ranks = 2 nodes (02-09 + 02-03), 1 GPU per rank, 3 blocks per rank**, multi-node UCX env of run_rsg.sh; one arm after the other.

## +3 step check, 8 ranks: ok=yes
- 2026-11-24 14:14:41 rc=0 (srun) ended: Terminating on cycle limit wall=51 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 75887 (both nodes, rank-0 samplers); CW 0, FATAL 0; 

## Abort lines of the 4-rank arms (your 24zzzzk item 1)
- H0 (srun exit 134): `  what():  Kokkos ERROR: Cuda memory space failed to allocate 4.254 GiB (label="cln_rk").` on every rank; it IS an out-of-memory allocation (4.254 GiB for cln_rk); peak before it 140457 MiB per GPU (nvidia-smi).
- H1 (srun exit 134): `  what():  Kokkos ERROR: Cuda memory space failed to allocate 4.254 GiB (label="cln_rk").` on every rank; it IS an out-of-memory allocation (4.254 GiB for cln_rk); peak before it 140869 MiB per GPU (nvidia-smi).
- HT (srun exit 134): `  what():  Kokkos ERROR: Cuda memory space failed to allocate 4.254 GiB (label="cln_rk").` on every rank; it IS an out-of-memory allocation (4.254 GiB for cln_rk); peak before it 136743 MiB per GPU (nvidia-smi).
- No unused arm holds GPU memory (the steps are gone); the +3 check below was run with H0's args (8 ranks), so it covers the H0 memory check: peak 75.9 GB per GPU.

## H0 / H1 / HT (11-24 14:26)
- H0: 2026-11-24 14:18:58 rc=0 (srun) ended: Terminating on cycle limit wall=227 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 75919; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.793294e+02 | failed stage at`
- H1: 2026-11-24 14:22:37 rc=0 (srun) ended: Terminating on cycle limit wall=218 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 75915; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.754426e+02 | failed stage at`
- HT: 2026-11-24 14:26:48 rc=0 (srun) ended: Terminating on cycle limit wall=249 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB 76059; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.888158e+02 | failed stage at`
| arm | rows | dt min/median(15-51)/max | s/cyc | s per sim s | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/lines | SDO lines, K1K2<0 mean share, max ratio, cells >0.1, >1, ef<=0 | L_top/L_in start..end |
|---|---|---|---|---|---|---|---|---|
| H0 | 51 | 90.00 / 90.00 / 90.00 | 3.720 | 0.0413 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 51, 0.3730, 2.957e-01, 2, 0, 0 | 0.9628 .. 0.9777 |
| H1 | 51 | 90.00 / 90.00 / 90.00 | 3.462 | 0.0385 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 51, 0.3728, 2.253e+00, 342270, 192, 0 | 0.9628 .. 0.9774 |
| HT | 51 | 90.00 / 90.00 / 90.00 | 3.950 | 0.0439 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 51, 0.3710, 7.822e+01, 40, 1, 1 | 0.9632 .. 0.9794 |

H1 vs H0 (reference H0), interpolated to H0 times:
| cycle +n | t [s] | dt H0 | dt H1 | L_top rel | L_mid rel | L_int rel | L_bot rel |
|---|---|---|---|---|---|---|---|
| 5 | 6.6096e+07 | 90.00 | 90.00 | -3.274e-07 | -1.291e-05 | 1.108e-05 | 1.914e-05 |
| 10 | 6.6097e+07 | 90.00 | 90.00 | -6.585e-06 | 6.134e-05 | -2.449e-04 | 4.259e-05 |
| 20 | 6.6098e+07 | 90.00 | 90.00 | 3.849e-05 | 7.630e-04 | 8.440e-05 | 2.479e-05 |
| 51 | 6.6101e+07 | 90.00 | 90.00 | -3.408e-04 | 1.126e-03 | 1.181e-03 | -3.610e-05 |
max |rel| over the run: L_top 5.395e-04 L_mid 1.738e-03 L_int 1.221e-03 L_bot 8.316e-05

HT vs H0 (reference H0), interpolated to H0 times:
| cycle +n | t [s] | dt H0 | dt HT | L_top rel | L_mid rel | L_int rel | L_bot rel |
|---|---|---|---|---|---|---|---|
| 5 | 6.6096e+07 | 90.00 | 90.00 | -1.035e-04 | -4.563e-03 | -5.102e-03 | -3.882e-10 |
| 10 | 6.6097e+07 | 90.00 | 90.00 | -5.297e-03 | -6.328e-03 | -3.666e-03 | -5.298e-11 |
| 20 | 6.6098e+07 | 90.00 | 90.00 | 2.175e-03 | 3.135e-03 | 3.153e-03 | -9.217e-10 |
| 51 | 6.6101e+07 | 90.00 | 90.00 | 1.700e-03 | -1.420e-03 | -1.370e-03 | -7.234e-10 |
max |rel| over the run: L_top 6.186e-03 L_mid 6.433e-03 L_int 6.532e-03 L_bot 9.217e-10

**PASS** (zzzr rules; H1 vs old O30 not tested here: needs O30.hst)
H0h needed (L_top diff 0.9-2 %): no

## Panel-edge |F| > cE (24zzzzc / zzzzg): +6-cycle vet dump, H0 and HT
```
R = 4.298e12 cm (61.78 Rsun, rst header)
== H0 hegiant.m1_vet.00154.bin cycle 794835
== HT hegiant.m1_vet.00154.bin cycle 794835
-- H0 hegiant.m1.02296.bin CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0000 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- H0 hegiant.m1.02296.bin CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0011 (N 3048) d 0.02-0.05: 0.0000 (N 5952) d 0.05-0.10: 0.0000 (N 8568) d 0.10-0.20: 0.0000 (N 18312) d 0.20-0.40: 0.0000 (N 27768) d 0.40-0.70: 0.0000 (N 25992) d 0.70-1.01: 0.0000 (N 8664)
-- H0 hegiant.m1.02296.bin CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0255 (N 3048) d 0.02-0.05: 0.0212 (N 5952) d 0.05-0.10: 0.0200 (N 8568) d 0.10-0.20: 0.0202 (N 18312) d 0.20-0.40: 0.0212 (N 27768) d 0.40-0.70: 0.0189 (N 25992) d 0.70-1.01: 0.0177 (N 8664)
-- H0 hegiant.m1.02296.bin CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6139 (N 3048) d 0.02-0.05: 0.6156 (N 5952) d 0.05-0.10: 0.6188 (N 8568) d 0.10-0.20: 0.6146 (N 18312) d 0.20-0.40: 0.6245 (N 27768) d 0.40-0.70: 0.6231 (N 25992) d 0.70-1.01: 0.6414 (N 8664)
-- HT hegiant.m1.02296.bin CSCALE=1 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0027 (N 3048) d 0.02-0.05: 0.0012 (N 5952) d 0.05-0.10: 0.0005 (N 8568) d 0.10-0.20: 0.0002 (N 18312) d 0.20-0.40: 0.0002 (N 27768) d 0.40-0.70: 0.0001 (N 25992) d 0.70-1.01: 0.0003 (N 8664)
-- HT hegiant.m1.02296.bin CSCALE=1 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.0079 (N 3048) d 0.02-0.05: 0.0025 (N 5952) d 0.05-0.10: 0.0003 (N 8568) d 0.10-0.20: 0.0002 (N 18312) d 0.20-0.40: 0.0001 (N 27768) d 0.40-0.70: 0.0002 (N 25992) d 0.70-1.01: 0.0008 (N 8664)
-- HT hegiant.m1.02296.bin CSCALE=0.5 range 1.1 1.5
hegiant.m1.02296.bin r/R 1.1-1.5 d 0.00-0.02: 0.0272 (N 3048) d 0.02-0.05: 0.0222 (N 5952) d 0.05-0.10: 0.0211 (N 8568) d 0.10-0.20: 0.0210 (N 18312) d 0.20-0.40: 0.0225 (N 27768) d 0.40-0.70: 0.0190 (N 25992) d 0.70-1.01: 0.0177 (N 8664)
-- HT hegiant.m1.02296.bin CSCALE=0.5 range 1.5 2
hegiant.m1.02296.bin r/R 1.5-2.0 d 0.00-0.02: 0.6107 (N 3048) d 0.02-0.05: 0.6124 (N 5952) d 0.05-0.10: 0.6132 (N 8568) d 0.10-0.20: 0.6110 (N 18312) d 0.20-0.40: 0.6241 (N 27768) d 0.40-0.70: 0.6226 (N 25992) d 0.70-1.01: 0.6407 (N 8664)
```
- dumps present; only the +6 one counts; they are deleted after this.

## Reading of H0 / H1 / HT and the status you asked for (24zzzzm), 11-24 14:28
- **H1 passes the 24zzzr rules as measured: dt 90.00 in all 51 cycles of both arms (no drop at ~14 cycles / 1170 s), no CW / FATAL / inf / DT_RETRY / fallback / BE / retry, H1 - H0 in L_top -3.4e-4 at +51 (max 5.4e-4), L_mid 1.7e-3, L_int 1.2e-3.** Early stop did not fire. H0h not needed (L_top diff < 0.9 %).
- What it does not show: ef<=0 is not a valid counter (your 24zzzzj); SDO max|ef/es-1| H1 2.25 (192 cells > 1), H0 0.30 (0 cells); the ratio is not a criterion. The comparison with old O30.hst is not done here (the file is not on Caltech). 51 cycles = 4590 s.
- **Cost (rule 4):** loop time H0 189.7 s, H1 176.6 s, HT 201.4 s for 51 cycles, but H0 and HT wrote 8 vet dumps (1.9 GB each) inside the loop and H1 none, so H1 / H0 = 0.93 is not a fair number and the rule is NOT judged from it. Solver work: see the its line below; rotation and dump-free stage time: H0 clean stage steps 1.793294e+02;H1 clean stage steps 1.754426e+02;
- **HT (thinf=true) vs H0:** dt 90.00; L_top +1.7e-3 at +51 (max 6.2e-3), L_mid max 6.4e-3; L_bot 1e-9; SDO max ratio 78 with 1 cell > 1 and ef<=0 cell 1 (HT is the binary with the key off; the counter's meaning as above).
- Panel-edge |F| > cE (+6-cycle vet dump, R 4.298e12 cm): H0 d < 0.02: 0.0000 (1.1-1.5 R) / 0.0011 (1.5-2 R); HT: 0.0027 / 0.0079. Positive control CSCALE = 0.5 rises (0.026 / 0.61 at d < 0.02), so the script counts. N 3048 columns per bin.
- Peak memory per GPU over the whole 51-cycle arm (nvidia-smi): H0 75919 MiB; H1 75915 MiB; HT 76059 MiB;  (8 ranks, 2 nodes 02-09 + 02-03). The rows > 200000 in the first post were torn lines from two samplers writing one file; corrected.
- Node table: HG done (02-09, 02-03 now free); uniform pairs running: U4 on 02-16 (A) + 02-17 (B), U3 A->B on 02-11 (mg 4), U1 A->B on 01-17, U2 A->B on 01-15. Nothing else running or queued.
- Cost from the dump-free clean stage time (rank 0): H0 179.3 s, H1 175.4 s (ratio 0.978), HT 188.8 s (1.053 x H0). Rule 4 (1.05 x H0) holds for H1 by this measure. Its per stage line: H0 16.3 / 14.8 (solves 4.31 / 5.41), H1 16.2 / 15.2 (4.22 / 5.61), HT 18.2 / 15.5 (5.59 / 6.22).

## Panel edge: which dump was scored (your 24zzzzo item 4)
- The scored files are hegiant.m1.02296.bin of H0 and HT. They **are the +6 dumps**: seam.sh only took files whose header says cycle = 794790 (= 794784 + 6), and 02296 is the first m1 dump of the run (written 14:16:21 for H0, against 14:18 for the end). The line 'cycle 794835' in the first listing belongs to hegiant.m1_vet.00154.bin (output6, the end-of-run final output), which was not scored. So the label stands: +6. The other dumps (and all bins) were deleted after the scoring, as you allowed; no second scoring is possible.
- **Corrected label, HT (thinf=true): NOT harmless on HG**: panel edge d < 0.02 0.0027 (1.1-1.5 R) and 0.0079 (1.5-2 R) against 0.0000 / 0.0011 in H0; L_top up to 0.62 %; ef<=0 in one cell; max|ef/es-1| 78.2 (H0 0.30).
- floor_pay columns follow in the H0F / H1F post (floor_pay sums of the old arms are yours: H0 4.1e31, H1 2.14e33, HT 9.1e31 erg).

## H0F / H1F (fixed binary sdorder2 09d1e4df4), 11-24 14:46
- Build: BUILD_OK md5=b8b0d79c45351fecfd04cee19c996605 prov=09d1e4df4486ded370dc4959dfd6b69e40258e8f; GATE PASS: athena_sdorder2_09d1e4df4 md5 b8b0d79c node hpc-sm-02-03 pass1=3.686e-11 end: Terminating on cycle limit; own dir wt_orion1121/sdo2/build_sdo2.
- H0F: 2026-11-24 14:42:15 rc=0 (srun) ended: Terminating on cycle limit wall=221 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB over the whole arm 75919; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.843740e+02 | failed stage at`
- H1F: 2026-11-24 14:46:01 rc=0 (srun) ended: Terminating on cycle limit wall=226 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB over the whole arm 76011; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.855392e+02 | failed stage at`
| arm | floor_pay sum (FALLBACK lines) | floor_pay cumulative (audit) | rows | dt min/median(15-51)/max | s/cyc | s per sim s | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/lines | SDO lines, K1K2<0 mean share, max ratio, cells >0.1, >1, ef<=0 | L_top/L_in start..end |
|---|---|---|---|---|---|---|---|---|---|---|
| H0F | 4.055e+31 | 4.036e+31 | 51 | 90.00 / 90.00 / 90.00 | 3.634 | 0.0404 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 51, 0.3730, 1.472e-01, 2, 0, 0 | 0.9628 .. 0.9777 |
| H1F | 2.138e+33 | 2.124e+33 | 51 | 90.00 / 90.00 / 90.00 | 3.659 | 0.0407 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 51, 0.3728, 2.291e+00, 342270, 186, 0 | 0.9628 .. 0.9775 |

H1F vs H0F (reference H0F), interpolated to H0F times:
| cycle +n | t [s] | dt H0F | dt H1F | L_top rel | L_mid rel | L_int rel | L_bot rel |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 6.6096e+07 | 90.00 | 90.00 | -3.056e-07 | -1.304e-05 | 9.567e-06 | 1.914e-05 |
| 10 | 6.6097e+07 | 90.00 | 90.00 | -6.877e-06 | 1.089e-04 | -2.492e-04 | 4.259e-05 |
| 20 | 6.6098e+07 | 90.00 | 90.00 | 7.489e-05 | 6.086e-04 | -3.504e-05 | 2.479e-05 |
| 51 | 6.6101e+07 | 90.00 | 90.00 | -2.094e-04 | 6.002e-04 | 7.593e-04 | -3.610e-05 |
max |rel| over the run: L_top 5.133e-04 L_mid 1.050e-03 L_int 1.576e-03 L_bot 8.316e-05

**FAIL: floor_pay 2.138e+33 > 1.2 x 4.055e+31 (ratio 52.73)** (zzzr rules; H1 vs old O30 not tested here: needs O30.hst)
H0h needed (L_top diff 0.9-2 %): no

## sdorder3 (your 24zzzzx), 11-24 15:48
- BUILD rc=0 wall=291 s TUs=250 BUILD_OK md5=30d01765970051212176cc2eb81ea034 prov=95813cf85f2901d61e9dc1320b8d54f480810682 ; .prov: 95813cf85f2901d61e9dc1320b8d54f480810682 dir /resnick/groups/carnegie_poc/jingze/wt_orion1121/sdo3/build_sdo3; cmake flags as build2.sh;  ; own dir wt_orion1121/sdo3/build_sdo3, flags as sdorder2; GATE PASS: athena_sdorder3_95813cf85 md5 30d01765 node hpc-sm-02-09 pass1=3.686e-11 end: Terminating on cycle limit
- Nothing runs with it until your explicit go. Prepared, not launched: hg_gate_1124/goG.PREPARED (H0G / H1G, 8 ranks, 51 cycles, floor_pay) and sdo_1124/nG.PREPARED (S3G, S5G as S3F / S5F).

## H0G / H1G (fixed binary sdorder3 95813cf85), 11-24 16:01
- Build: BUILD_OK md5=30d01765970051212176cc2eb81ea034 prov=95813cf85f2901d61e9dc1320b8d54f480810682; GATE PASS: athena_sdorder3_95813cf85 md5 30d01765 node hpc-sm-02-09 pass1=3.686e-11 end: Terminating on cycle limit; own dir wt_orion1121/sdo3/build_sdo3.
- H0G: 2026-11-24 15:57:05 rc=0 (srun) ended: Terminating on cycle limit wall=213 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB over the whole arm 75807; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.808108e+02 | failed stage at`
- H1G: 2026-11-24 16:00:47 rc=0 (srun) ended: Terminating on cycle limit wall=222 s CW 0 FATAL 0 inf 0 DT_RETRY 0 fallback 0; peak GPU memory MiB over the whole arm 76059; last log line: `<rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 1.802509e+02 | failed stage at`
| arm | floor_pay sum (FALLBACK lines) | floor_pay cumulative (audit) | rows | dt min/median(15-51)/max | s/cyc | s per sim s | CW/FATAL/inf/DT_RETRY | fallbacks/BE/retries/lines | SDO lines, K1K2<0 mean share, max ratio, cells >0.1, >1, ef<=0 | L_top/L_in start..end |
|---|---|---|---|---|---|---|---|---|---|---|
| H0G | 4.055e+31 | 4.036e+31 | 51 | 90.00 / 90.00 / 90.00 | 3.566 | 0.0396 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0, nan, nan, nan, nan, nan | 0.9628 .. 0.9777 |
| H1G | 1.437e+29 | 1.414e+29 | 51 | 90.00 / 90.00 / 90.00 | 3.553 | 0.0395 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0, nan, nan, nan, nan, nan | 0.9628 .. 0.9775 |

H1G vs H0G (reference H0G), interpolated to H0G times:
| cycle +n | t [s] | dt H0G | dt H1G | L_top rel | L_mid rel | L_int rel | L_bot rel |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 6.6096e+07 | 90.00 | 90.00 | -2.263e-07 | -1.175e-05 | 1.076e-05 | 1.914e-05 |
| 10 | 6.6097e+07 | 90.00 | 90.00 | -6.462e-06 | 7.675e-05 | -2.295e-04 | 4.259e-05 |
| 20 | 6.6098e+07 | 90.00 | 90.00 | 3.044e-05 | 7.072e-04 | 5.106e-04 | 2.479e-05 |
| 51 | 6.6101e+07 | 90.00 | 90.00 | -2.027e-04 | 1.600e-03 | 1.122e-03 | -3.610e-05 |
max |rel| over the run: L_top 4.638e-04 L_mid 1.607e-03 L_int 1.395e-03 L_bot 8.316e-05

**PASS** (zzzr rules; H1 vs old O30 not tested here: needs O30.hst)
H0h needed (L_top diff 0.9-2 %): no

## H0G / H1G remarks
- Dump-free clean stage time (rank 0): H0G 180.8 s, H1G 180.3 s for 51 cycles (ratio 0.997). Nodes 01-17 + 02-09 (the b_suite CPU step had ended before H0G started). floor_audit paid cumulative at cycle 51: H0G 4.036e31, H1G 1.414e29.
- **floor_pay H1G/H0G = 0.0035** (rule <= 1.2): PASS, the only rule H1F failed. SDO lines are not available on this binary (SDO_PRINT removed), so ef<=0 and the ratio are not measured.

## c2 build (your 24zzzza6): the branch is NOT on GitHub (11-24 16:35)
- `git fetch origin caltech-c2-d4ae795` -> 'couldn't find remote ref'; `git ls-remote origin` shows no ref named caltech-c2* and no d4ae795 / 6211aaa56. Nothing was built (no worktree, no binary). Please push the branch; I build it on 01-15 (CPU step, flags as sdorder3) and post the source diff against 95813cf85, the tiny gate, md5 and .prov. No arms until your gates pass.
