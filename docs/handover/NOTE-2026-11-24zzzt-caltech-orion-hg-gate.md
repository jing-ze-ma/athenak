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
