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
