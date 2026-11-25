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

## t9l62, 7b4141e1b against 6d20f1337 (ABBA, 60 cycles, 02-16, own keys rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=2 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; ndiag=1 so that every cycle is timed)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 58 | 2.6116 | 2.6206 | 2.509/2.652 | 54.393/56.118 | 18.6 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=225 s CW 0 FATAL 0 |
| NEW1 | 58 | 2.6097 | 2.6174 | 2.500/2.676 | 54.393/56.118 | 18.6 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=224 s CW 0 FATAL 0 |
| NEW2 | 58 | 2.6089 | 2.6182 | 2.504/2.650 | 54.393/56.118 | 18.6 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=227 s CW 0 FATAL 0 |
| BASE2 | 58 | 2.6107 | 2.6197 | 2.502/2.653 | 54.393/56.118 | 18.6 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=221 s CW 0 FATAL 0 |
- base (both arms) mean 2.6111 median 2.6200; new mean 2.6093 median 2.6177; new / base mean 0.9993, median 0.9991; base arms 2.6116 vs 2.6107, new arms 2.6097 vs 2.6089 (mean)
- NEW1 user.hst against BASE1: 60 / 60 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 60 / 60 rows, printed-bitwise
- NEW2 user.hst against BASE1: 60 / 60 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 60 / 60 rows, printed-bitwise
- BASE2 user.hst against BASE1: 60 / 60 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 60 / 60 rows, printed-bitwise

## HG 2_7b4141e1b, 7b4141e1b against 6d20f1337: stopped by 25l after NEW2 (BASE2 and the pairs 3-5 not run)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/b1_1124/sbab2.py", line 14, in <module>
    f = glob.glob(S + a + '/*.user.hst'); g = glob.glob(S + a + '/*.hydro.hst'); H[a] = [[l for l in open(x) if l.strip() and l[0] != '#'] for x in (f + g)]
NameError: name 'H' is not defined
srun: error: hpc-sm-02-03: task 0: Exited with exit code 1

## m35mms, 7b4141e1b against 6d20f1337 (ABBA, 60 cycles, 02-16, own keys rad_vet/time2_gtime_slope=true rad_vet/clean_stag_thinf=true rad_vet/clean_stag_lat_sums=true rad_vet/clean_stag_tie_gid=true rad_vet/clean_stag_opac_drag=true rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8 rad_vet/clean_pc_sweep=rb rad_vet/clean_stag_mg_levels=4 time/cfl_number=0.3 rad_vet/clean_conv_max=1e-5; ndiag=1 so that every cycle is timed)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 58 | 1.0673 | 1.0636 | 1.025/1.141 | 12.007/12.007 | 61.3 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=94 s CW 0 FATAL 0 |
| NEW1 | 58 | 1.0684 | 1.0617 | 1.027/1.146 | 12.007/12.007 | 61.3 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=99 s CW 0 FATAL 0 |
| NEW2 | 58 | 1.0671 | 1.0600 | 1.027/1.146 | 12.007/12.007 | 61.3 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=97 s CW 0 FATAL 0 |
| BASE2 | 58 | 1.0670 | 1.0599 | 1.028/1.140 | 12.007/12.007 | 61.3 | 0/0/0/0 | rc=0 (srun) ended: Terminating on cycle limit wall=98 s CW 0 FATAL 0 |
- base (both arms) mean 1.0672 median 1.0612; new mean 1.0677 median 1.0616; new / base mean 1.0005, median 1.0003; base arms 1.0673 vs 1.0670, new arms 1.0684 vs 1.0671 (mean)
- NEW1 user.hst against BASE1: 60 / 60 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 60 / 60 rows, printed-bitwise
- NEW2 user.hst against BASE1: 60 / 60 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 60 / 60 rows, printed-bitwise
- BASE2 user.hst against BASE1: 60 / 60 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 60 / 60 rows, printed-bitwise
