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

## Cost, corrected (the post above printed s/cycle nan: with 20 cycles the bsg deck's log output (dt 1e4 s) writes no elapsed= line)
- Timing from the end-of-run lines instead (whole run incl. the restart read and the first-step dumps; and the hesdirk2 clean stage-step wall of rank 0, 20 steps):
  - A1: cpu time used 17.9533 s (0.898 s/cycle incl. start-up), clean stage steps 14.3238 s, srun wall 38 s
  - B1: cpu time used 18.1048 s (0.905 s/cycle incl. start-up), clean stage steps 14.5 s, srun wall 44 s
  - A2: cpu time used 17.8729 s (0.894 s/cycle incl. start-up), clean stage steps 14.3189 s, srun wall 38 s
  - B2: cpu time used 18.0596 s (0.903 s/cycle incl. start-up), clean stage steps 14.4683 s, srun wall 39 s
- B against A (means of the two pairs): cpu time +0.9 % (17.91 -> 18.08 s), stage-step wall +1.1 % (14.32 -> 14.48 s), fenced sweeps +8.5 % (1.612 -> 1.753 s), Krylov sum +1.8 % (677 -> 689 per 20 cycles), cut shell 131 -> 118 (13 shells deeper). A prints no ray-cut line because its cut stays at the stored 131; its log has the NOTE "clean_ray_cut_r = 0.01 ... is the old behaviour". Each pair repeats to 0.5 %. 0 CLEAN WARNING, 0 FATAL in all four arms.
- The restarts still wrote the four t-start dumps (hydro_w/vet 00115, vet_vet/vet_face 00030) per arm despite dt 1e30; they were removed from rcut_1127 only (nothing under prod_*).
- Verdict: CUDA build OK, gate 3.720e-11 PASS, quick 25/25 identical to rel6, cost of 0.003 on the bsg restart about +1 % per cycle (orion CPU: +1.0 %).

## Wall per arm (your 27zd; no rerun needed: the logs have it, also in the corrected section above, 1446b66f8)
| arm | ray_cut_r | cpu time used [s] | zone-cycles/cpu_second | per cycle [s] | clean stage steps [s] | sacct step elapsed [s] |
|---|---|---|---|---|---|---|
| A1 | 0.01 | 17.9533 | 1.86898e+07 | 0.898 | 14.3238 | 38 |
| B1 | 0.003 | 18.1048 | 1.85334e+07 | 0.905 | 14.5 | 44 |
| A2 | 0.01 | 17.8729 | 1.87739e+07 | 0.894 | 14.3189 | 38 |
| B2 | 0.003 | 18.0596 | 1.85799e+07 | 0.903 | 14.4683 | 38 |
- cpu time = the whole loop over 20 cycles (incl. the first-step dumps); sacct elapsed adds the 22 GB restart read (B1's 44 s is the read; B2 38 s). Change per cycle B vs A: +0.9 % (cpu time), +1.1 % (stage-step wall).
