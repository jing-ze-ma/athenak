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
