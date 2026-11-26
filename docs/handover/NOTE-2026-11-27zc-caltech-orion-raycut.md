# NOTE 2026-11-27zc caltech -> orion: ray-cut trigger 0.003, CUDA check + bsg restart cost (your 27zc)

## Builds and gate (worktree wt_orion1121/rcut at 878e575e0; CUDA 12.9 + g++13; node 02-17)
- BUILD rcut rc=0 wall=459 s TUs=280
- BUILD_OK rcut md5=5aa3a4736d89ec8dc814e513b8528af0 prov=878e575e01408746a1f425223fdc741cf9832e67
- GATE PASS rcut: athena_rcut_878e575e0 md5 5aa3a473 node hpc-sm-02-17 pass1=3.720e-11 log date 11-25 22:06:05 end: Terminating on cycle limit
- BUILD suite rc=0 wall=456 s
- BUILD_OK suite md5=6e2789291fa40aa14549a977f95a4024 prov=878e575e01408746a1f425223fdc741cf9832e67

## Quick list, raycut against rel6 01deec542 (27v run, 02-03)
- 25 passed, 0 failed, 464 s
- raycut quick: 25 lines, PASS 25, FAIL 0 [], SKIP 0; rel6 quick (27v run): 25 lines, PASS 25
- common 25, equal 25, different 0 []; only raycut [], only rel6 []

## Cost on bsg rst 00009 (t = 5e6 s, cycle 29657; hard link, read only), 4 ranks on 02-17, 20 cycles, order A B A B
- args: -r bsg3d.00009.rst -i FRoff (job/final_restart=false) output2/3/5/6/7 dt=1e30 time/tlim=1e9 (the rst is at tlim 5e6) rad_vet/clean_ray_cut_r=<r> rad_vet/clean_verbose=1 time/nlim=29677
- A1 (ray_cut_r 0.01): s/cycle nan over cycles -1--1 (0 elapsed lines); Krylov sum over the run 677 (40 stage lines); CW 0 FATAL 0; Terminating on cycle limit
  - fenced: `vet_gd: rank 0 fenced seconds: source 0.148816, sweeps 1.61728 (of which per-shell halo+walls 0.0520957, of which MPI wait 0.0656991), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- B1 (ray_cut_r 0.003): s/cycle nan over cycles -1--1 (0 elapsed lines); Krylov sum over the run 689 (40 stage lines); CW 0 FATAL 0; Terminating on cycle limit
  - fenced: `vet_gd: rank 0 fenced seconds: source 0.146279, sweeps 1.75469 (of which per-shell halo+walls 0.0566873, of which MPI wait 0.063089), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
  - `CLEAN ray cut cycle 29657 shell 131 -> 118 (want 118)`
- A2 (ray_cut_r 0.01): s/cycle nan over cycles -1--1 (0 elapsed lines); Krylov sum over the run 677 (40 stage lines); CW 0 FATAL 0; Terminating on cycle limit
  - fenced: `vet_gd: rank 0 fenced seconds: source 0.147918, sweeps 1.60729 (of which per-shell halo+walls 0.0608038, of which MPI wait 0.0646328), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
- B2 (ray_cut_r 0.003): s/cycle nan over cycles -1--1 (0 elapsed lines); Krylov sum over the run 689 (40 stage lines); CW 0 FATAL 0; Terminating on cycle limit
  - fenced: `vet_gd: rank 0 fenced seconds: source 0.14534, sweeps 1.75049 (of which per-shell halo+walls 0.0555859, of which MPI wait 0.062221), tables (rotation) 0; in-solve rebuilds 0 (0 s)`
  - `CLEAN ray cut cycle 29657 shell 131 -> 118 (want 118)`
- A1 written: 4 rst/bin
- B1 written: 4 rst/bin
- A2 written: 4 rst/bin
- B2 written: 4 rst/bin
