# NOTE 2026-11-27zh caltech -> orion: scpr, PR #798 short characteristics vs our rays on H200 (your 27zh)

## Build: PR #798 c5ac58778 (fetched read-only into its own worktree wt_orion1121/pr798, nothing pushed), its Kokkos submodule 4.7.02, CUDA 12.9 + g++13, Release, MPI, HOPPER90, -ffp-contract=off, PROBLEM=built_in_pgens; built as a CPU step on 02-17
- BUILD pr798 rc=0 wall=214 s
- BUILD_OK md5=a77bc7ee11e5e4017e48cf85532ca5f8 prov=c5ac5877899fe95ffe89f042139f4ffa486443e4 kokkos 4.7.02 6739bc623

## Runs (02-17, one rank per GPU)
- Ours (G1 A cart, G1 B sp, G4 A cart, with the timer) = the 27zg runs on this node in this session (scspeed_1127/runs/G1_*, kp_G1_*, kp_G4_A; posted in NOTE-27ze). Our G4 B FATALs (theta mirror wall).
```
S0 rc=0 wall=6 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput time/nlim=1
G1_P1a rc=1 wall=10 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 problem/radius=10.0 problem/sigma_a=16.0 problem/r_mask=0.0 problem/tol=1e30 output1/dt=1e30 time/nlim=10 nr_radiation/iter_min=1 nr_radiation/iter_max=1
G1_P2a rc=1 wall=5 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 problem/radius=10.0 problem/sigma_a=16.0 problem/r_mask=0.0 problem/tol=1e30 output1/dt=1e30 time/nlim=10 nr_radiation/iter_min=1 nr_radiation/iter_max=1 nr_radiation/sweep_kernel=tiled nr_radiation/tile_size=16
G1_P1b rc=1 wall=5 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 problem/radius=10.0 problem/sigma_a=16.0 problem/r_mask=0.0 problem/tol=1e30 output1/dt=1e30 time/nlim=10 nr_radiation/iter_min=1 nr_radiation/iter_max=1
G1_P2b rc=1 wall=4 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 problem/radius=10.0 problem/sigma_a=16.0 problem/r_mask=0.0 problem/tol=1e30 output1/dt=1e30 time/nlim=10 nr_radiation/iter_min=1 nr_radiation/iter_max=1 nr_radiation/sweep_kernel=tiled nr_radiation/tile_size=16
G4_P1 rc=1 wall=14 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 problem/radius=10.0 problem/sigma_a=16.0 problem/r_mask=0.0 problem/tol=1e30 output1/dt=1e30 time/nlim=10 nr_radiation/iter_min=1 nr_radiation/iter_max=1
G4_it rc=0 wall=10 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64 time/nlim=5 output1/dt=1e30 nr_radiation/iter_tol=1e-6
A64_pr rc=0 wall=5 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput
A64_our rc=0 wall=6 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/our_sph.athinput
A32_pr rc=0 wall=11 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=32 mesh/nx2=32 mesh/nx3=32 meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=32
A32_our rc=0 wall=10 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/our_sph.athinput mesh/nx1=32 mesh/nx2=32 mesh/nx3=32 meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=32
```
### S0: s/cycle nan over cycles ; cpu time used  = 7.190183e-02
  - `Iteration stops at niter: 2 relative error: 0.000000e+00`
  - `SC sphere: checked 6272 interior cells (r < 6.000000e-01R), rel L2 err  J = 6.501907e-03  H = 2.106668e-03  K = 3.756178e-03  (J Linf 8.910804e-03; 2 sweeps, residual 0.000000e+00)`
  - `Terminating on cycle limit`
### G1_P1a: s/cycle nan over cycles ; 
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/pr798/src/parameter_input.cpp at line 400 | Parameter 'iter_min' in block 'nr_radiation' on command line not found`
### G1_P2a: s/cycle nan over cycles ; 
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/pr798/src/parameter_input.cpp at line 400 | Parameter 'iter_min' in block 'nr_radiation' on command line not found`
### G1_P1b: s/cycle nan over cycles ; 
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/pr798/src/parameter_input.cpp at line 400 | Parameter 'iter_min' in block 'nr_radiation' on command line not found`
### G1_P2b: s/cycle nan over cycles ; 
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/pr798/src/parameter_input.cpp at line 400 | Parameter 'iter_min' in block 'nr_radiation' on command line not found`
### G4_P1: s/cycle nan over cycles ; 
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/pr798/src/parameter_input.cpp at line 400 | Parameter 'iter_min' in block 'nr_radiation' on command line not found`
  - `### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/pr798/src/parameter_input.cpp at line 400 | Parameter 'iter_min' in block 'nr_radiation' on command line not found`
  - (4 lines in all)
### G4_it: s/cycle nan over cycles ; cpu time used  = 5.481628e-01
  - `Iteration stops at niter: 4 relative error: 0.000000e+00`
  - `SC sphere: checked 12860 interior cells (r < 6.000000e-01R), rel L2 err  J = 2.957478e-03  H = 9.533103e-04  K = 1.708034e-03  (J Linf 4.101422e-03; 4 sweeps, residual 0.000000e+00)`
  - `SC sphere: checked 12860 interior cells (r < 0.6R), rel L2 err  J = 0.00295748  H = 0.00095331  K = 0.00170803  (J Linf 0.00410142; 4 sweeps, residual 0)`
  - (4 lines in all)
  - `Terminating on time limit`
### A64_pr: s/cycle nan over cycles ; cpu time used  = 1.529642e-01
  - `Iteration stops at niter: 2 relative error: 0.000000e+00`
  - `SC sphere: checked 6272 interior cells (r < 6.000000e-01R), rel L2 err  J = 6.501907e-03  H = 2.106668e-03  K = 3.756178e-03  (J Linf 8.910804e-03; 2 sweeps, residual 0.000000e+00)`
  - `Terminating on cycle limit`
### A64_our: s/cycle nan over cycles ; cpu time used  = 1.861305e-01
  - `<rad_vet> vet_sc: formal solutions=2.000000e+00 SC seconds=8.358986e-02 (4.179493e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.358986e+298`
  - `Terminating on cycle limit`
### A32_pr: s/cycle nan over cycles ; cpu time used  = 1.778433e-01
  - `Iteration stops at niter: 2 relative error: 0.000000e+00`
  - `SC sphere: checked 840 interior cells (r < 6.000000e-01R), rel L2 err  J = 1.436257e-02  H = 4.789458e-03  K = 8.301927e-03  (J Linf 2.026437e-02; 2 sweeps, residual 0.000000e+00)`
  - `Terminating on cycle limit`
### A32_our: s/cycle nan over cycles ; cpu time used  = 1.369426e-01
  - `<rad_vet> vet_sc: formal solutions=2.000000e+00 SC seconds=2.226134e-02 (1.113067e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=2.226134e+298`
  - `Terminating on cycle limit`
### Accuracy (sph_err.py)
```
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/scpr_1127/sph_err.py", line 14, in <module>
    J = np.asarray(f['mb_data'][var][0], dtype=float)
                   ~~~~~~~~~~~~^^^^^
KeyError: 'sc_rad_J'
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/scpr_1127/sph_err.py", line 15, in <module>
    x = np.asarray(f['mb_x1v'][0]); y = np.asarray(f['mb_x2v'][0]); z = np.asarray(f['mb_x3v'][0])
                   ~^^^^^^^^^^
KeyError: 'mb_x1v'
srun: error: hpc-sm-02-17: task 0: Exited with exit code 1
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/scpr_1127/sph_err.py", line 14, in <module>
    J = np.asarray(f['mb_data'][var][0], dtype=float)
                   ~~~~~~~~~~~~^^^^^
KeyError: 'sc_rad_J'
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/scpr_1127/sph_err.py", line 15, in <module>
    x = np.asarray(f['mb_x1v'][0]); y = np.asarray(f['mb_x2v'][0]); z = np.asarray(f['mb_x3v'][0])
                   ~^^^^^^^^^^
KeyError: 'mb_x1v'
srun: error: hpc-sm-02-17: task 0: Exited with exit code 1
```
