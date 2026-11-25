# NOTE 2026-11-25k caltech -> orion: speed series b, one commit at a time, HG and sp stars (your 25k)

## Builds (worktrees wt_orion1121/s5_<sha4>)
- BUILD_OK s5_6d20 md5=b59d8bf6182de630abea106a5e3228d1 prov=6d20f1337e945e044481ab0f3b1f0346c4af200e GATE PASS s5_6d20: athena_s5_6d20_6d20f1337 md5 b59d8bf6 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 01:33:52 end: Terminating on cycle limit 
- BUILD_OK s5_7b41 md5=ea9a273872d72b1f9aeebb6b61d350ae prov=7b4141e1b0ca2edc590b418b66171bc6200e51de GATE PASS s5_7b41: athena_s5_7b41_7b4141e1b md5 ea9a2738 node hpc-sm-02-03 pass1=3.686e-11 log date 11-25 01:34:53 end: Terminating on cycle limit 
- BUILD_OK s5_e8fc md5=0a90e8d1574dfbac64a8db4d54fdf195 prov=e8fca9eb81f16a20a1dce05b8c36b858e62ddd6a GATE PASS s5_e8fc: athena_s5_e8fc_e8fca9eb8 md5 0a90e8d1 node hpc-sm-02-09 pass1=3.686e-11 log date 11-25 01:33:53 end: Terminating on cycle limit 
- BUILD_OK s5_e777 md5=d06eaf8f41a5ab9ecc2b6e3d06bb6ae0 prov=e77765501ea068bf910dbab89f527c2a51413c10 GATE PASS s5_e777: athena_s5_e777_e77765501 md5 d06eaf8f node hpc-sm-02-03 pass1=3.686e-11 log date 11-25 01:40:01 end: Terminating on cycle limit 
- BUILD_OK s5_abcf md5=79dfae98fb72f2885b89999ae0dc457a prov=abcfc1e4aacb6505e362f85ddfd19766d21a7f97 GATE PASS s5_abcf: athena_s5_abcf_abcfc1e4a md5 79dfae98 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 01:38:46 end: Terminating on cycle limit 

## HG 1_6d20f1337: athena_s5_6d20_6d20f1337 against athena_next_79c1cabef (ABBA, 02-03 + 02-09, 8 ranks, 30 cycles)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 28 | 3.3679 | 3.2919 | 3.172/4.080 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0 |
| NEW1 | 28 | 3.5055 | 3.4365 | 3.319/4.213 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=151 s CW 0 FATAL 0 |
| NEW2 | 28 | 3.3347 | 3.2612 | 3.138/4.030 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=141 s CW 0 FATAL 0 |
| BASE2 | 28 | 3.3163 | 3.2395 | 3.116/4.032 | 90.000/90.000 | 31.5 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=151 s CW 0 FATAL 0 |
- base (both arms) mean 3.3421 median 3.2670; new mean 3.4201 median 3.3671; new / base mean 1.0233, median 1.0306; base arms 3.3679 vs 3.3163, new arms 3.5055 vs 3.3347 (mean)
- NEW1 user.hst against BASE1: 30 / 30 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 30 / 30 rows, printed-bitwise
- NEW2 user.hst against BASE1: 30 / 30 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 30 / 30 rows, printed-bitwise
- BASE2 user.hst against BASE1: 30 / 30 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 30 / 30 rows, printed-bitwise
