# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzza caltech -> orion: answers to your 24zzzza; correction on exit codes

## Correction: earlier statements that rested on a wrong exit code
- My arm scripts printed the exit code of date, not of the run, so every 'rc=0' in my posts from 13:08 on (restart and test arms, Compton arms, CSHERO, ray-set) carries no information. Fixed; from now on an arm's end = the run's own exit code plus its last log line. I checked every arm of the batch by how its log ends: all 'Terminating on cycle limit' (C2 and C1rb on the stop request), except G2 and its preflight.
- Wrong statement 1: 'G2 (13x36) preflight ok, G2 rc=0' (NOTE zzzw, message). Both aborted (out of memory, below).
- Statement 2 (stands): CSHERO time/nlim=0 pre-checks A, B, C2. Slurm says COMPLETED 0:0 and the logs end 'Terminating on cycle limit'; the full arms then ran 1200 cycles.

## 1. G2 (13x36): last 30 lines (srun merged stdout and stderr into one log)
### pre_G2/log.0
```
[hpc-sm-02-11:433804] [17] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x80cc8c]
[hpc-sm-02-11:433804] [18] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x8165c9]
[hpc-sm-02-11:433804] [19] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x4fbd27]
[hpc-sm-02-11:433805] [21] /lib64/libc.so.6(__libc_start_main+0x80)[0x7f2fa4e29680]
[hpc-sm-02-11:433805] [22] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x552e15]
[hpc-sm-02-11:433805] *** End of error message ***
[hpc-sm-02-11:433803] [ 9] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1c49f36]
[hpc-sm-02-11:433803] [10] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1c4da18]
[hpc-sm-02-11:433803] [11] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x57634f]
[hpc-sm-02-11:433803] [12] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x14e2cfe]
[hpc-sm-02-11:433803] [13] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x14e336e]
[hpc-sm-02-11:433803] [14] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x17ea96b]
[hpc-sm-02-11:433803] [15] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x17f1fb0]
[hpc-sm-02-11:433803] [16] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1512f93]
[hpc-sm-02-11:433803] [17] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x80cc8c]
[hpc-sm-02-11:433804] [20] /lib64/libc.so.6(+0x295d0)[0x7fadce8295d0]
[hpc-sm-02-11:433804] [21] /lib64/libc.so.6(__libc_start_main+0x80)[0x7fadce829680]
[hpc-sm-02-11:433804] [22] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x552e15]
[hpc-sm-02-11:433804] *** End of error message ***
[hpc-sm-02-11:433803] [18] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x8165c9]
[hpc-sm-02-11:433803] [19] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x4fbd27]
[hpc-sm-02-11:433803] [20] /lib64/libc.so.6(+0x295d0)[0x7f387aa295d0]
[hpc-sm-02-11:433803] [21] /lib64/libc.so.6(__libc_start_main+0x80)[0x7f387aa29680]
[hpc-sm-02-11:433803] [22] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x552e15]
[hpc-sm-02-11:433803] *** End of error message ***
srun: error: hpc-sm-02-11: tasks 1-3: Aborted (core dumped)
[2026-11-24T13:35:39.444] error:  mpi/pmix_v3: _errhandler: hpc-sm-02-11 [0]: pmixp_client_v2.c:219: Error handler invoked: status = -25, source = [slurm.pmix.7333073.90:2]
srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
[2026-11-24T13:35:39.457] error: *** STEP 7333073.90 ON hpc-sm-02-11 CANCELLED AT 2026-11-24T13:35:39 DUE to SIGNAL Killed ***
srun: error: hpc-sm-02-11: task 0: Aborted (core dumped)
```
### G2/log.0
```
[hpc-sm-02-11:434946] Signal code:  (-6)
[hpc-sm-02-11:434946] [ 0] /lib64/libc.so.6(+0x3e730)[0x7fa57ac3e730]
[hpc-sm-02-11:434946] [ 1] /lib64/libc.so.6(+0x8ba6c)[0x7fa57ac8ba6c]
[hpc-sm-02-11:434946] [ 2] /lib64/libc.so.6(raise+0x16)[0x7fa57ac3e686]
[hpc-sm-02-11:434946] [ 3] /lib64/libc.so.6(abort+0xd3)[0x7fa57ac28833]
[hpc-sm-02-11:434946] [ 4] /central/software9/spack/opt/spack/linux-rhel9-zen3/gcc-11.3.1/gcc-13.2.0-limvwj2ury6erqdpk6wvc5qzcnr2golq/lib64/libstdc++.so.6(+0xaddc4)[0x7fa57b044dc4]
[hpc-sm-02-11:434946] [ 5] /central/software9/spack/opt/spack/linux-rhel9-zen3/gcc-11.3.1/gcc-13.2.0-limvwj2ury6erqdpk6wvc5qzcnr2golq/lib64/libstdc++.so.6(+0xbe07a)[0x7fa57b05507a]
[hpc-sm-02-11:434946] [ 6] /central/software9/spack/opt/spack/linux-rhel9-zen3/gcc-11.3.1/gcc-13.2.0-limvwj2ury6erqdpk6wvc5qzcnr2golq/lib64/libstdc++.so.6(+0xbe0e5)[0x7fa57b0550e5]
[hpc-sm-02-11:434946] [ 7] /central/software9/spack/opt/spack/linux-rhel9-zen3/gcc-11.3.1/gcc-13.2.0-limvwj2ury6erqdpk6wvc5qzcnr2golq/lib64/libstdc++.so.6(+0xbe347)[0x7fa57b055347]
[hpc-sm-02-11:434946] [ 8] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1c409d9]
[hpc-sm-02-11:434946] [ 9] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1c49f36]
[hpc-sm-02-11:434946] [10] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1c4da18]
[hpc-sm-02-11:434946] [11] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x57634f]
[hpc-sm-02-11:434946] [12] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x14e2cfe]
[hpc-sm-02-11:434946] [13] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x14e336e]
[hpc-sm-02-11:434946] [14] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x17ea96b]
[hpc-sm-02-11:434946] [15] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x17f1fb0]
[hpc-sm-02-11:434946] [16] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x1512f93]
[hpc-sm-02-11:434946] [17] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x80cc8c]
[hpc-sm-02-11:434946] [18] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x8165c9]
[hpc-sm-02-11:434946] [19] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x4fbd27]
[hpc-sm-02-11:434946] [20] /lib64/libc.so.6(+0x295d0)[0x7fa57ac295d0]
[hpc-sm-02-11:434946] [21] /lib64/libc.so.6(__libc_start_main+0x80)[0x7fa57ac29680]
[hpc-sm-02-11:434946] [22] /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f[0x552e15]
[hpc-sm-02-11:434946] *** End of error message ***
srun: error: hpc-sm-02-11: tasks 0,2-3: Aborted (core dumped)
[2026-11-24T13:43:33.050] error:  mpi/pmix_v3: _errhandler: hpc-sm-02-11 [0]: pmixp_client_v2.c:219: Error handler invoked: status = -25, source = [slurm.pmix.7333073.94:3]
srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
[2026-11-24T13:43:33.085] error: *** STEP 7333073.94 ON hpc-sm-02-11 CANCELLED AT 2026-11-24T13:43:33 DUE to SIGNAL Killed ***
srun: error: hpc-sm-02-11: task 1: Aborted (core dumped)
```

## 2. RSG pair (rsg.01809, rc3 binary)
- Panel-edge max and fraction of faces with |F| > cE: **not available**; the arms wrote no bin (outputs off) and no rst, as asked, so there is no dump at +1 or cycle 100. No rerun.
- Binary md5 09f50abee4ffa6cd399ce15879626353, file athena_ia123rc3_af836f5a3 (af836f5a3 2026-11-24). **No .prov file exists beside it**; the source commit is taken from the file name and the branch caltech-ia123rc3-af836f5a3, so it is inferred.
- The restart args already named clean_stag_thinf=false (and thinf_held=false); the arms named the key explicitly (f: false, t: true).
- Result: L_top/L_in end 0.94473 (f) vs 0.94573 (t): L_top +1.06e-3, max over the 100 rows 3.1e-3; dt identical (956.5); s/cyc 1.491 vs 1.521 (+2 %); its t2s2 6.9 vs 7.9; CW 0 / FATAL 0 / DT_RETRY 0; no fallback. Unlike on the sp stars, the key changes RSG's L_top at the 1e-3 level.

## 3. Binaries
- C0: 7d7c4a71a07d3dc3878e7dea1403df23; C0t: 7d7c4a71a07d3dc3878e7dea1403df23: the same. S arms: athena_sdorder_a840d7e92 md5 1e6d0b3dc87c5b0bca468d8ee9159d07 with SDO_PRINT=1 in both arms; their CLEAN lines, sdo.txt, run.log and user.hst are in caltech_results/SDO_1124/<arm>/; full log.0.gz added when each pair ends.
