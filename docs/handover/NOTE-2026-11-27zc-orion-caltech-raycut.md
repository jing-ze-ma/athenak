# NOTE 2026-11-27zc orion -> Caltech: ray-cut trigger 0.003, CUDA check + cost on a 3-D star restart

Branch caltech-raycut-878e575e0 = release aa2646b92 + one commit (USER-approved): the ray cut's
trigger is the key `clean_ray_cut_r` again, default 0.003 (was the constant 0.01). The cut starts
the sweep 10 optical depths below each column's lowest cell with |grad S|/(chi S) >= r, so 0.003
sweeps more shells (star_small: L_top error 1.9e-5 -> 1.3e-7). Restart rule: a restart file that
holds `clean_ray_cut_r = 0.01` keeps 0.01 (NOTE in the log); a file without the key (every rst of
the release) takes 0.003 (NOTE: "restart file without it: the fix applies").

Orion CPU (all green): build, tiny gate tiny_bsg_rel5 PASS 3.720e-11 (unchanged), quick 25/25,
full 61/61, 0 CLEAN WARNING; bitwise restart at 0.003, layout invariance (1/2/4 MeshBlocks)
bitwise, an rst with 0.01 continues bitwise equal to the release binary. Orion cost (star3d_sp
256x64x64, 16 ranks x 4 threads, 30 cycles, arms interleaved twice): cut shell 147 -> 136
(+10 % swept shells), sweep time +4.5 %, s/cycle +1.0 % (3.75 -> 3.79 s).

Please, on one free H200 node:
1. CUDA 12.9 + g++13 build of caltech-raycut-878e575e0; gate deck tiny_bsg_rel5 (3.720e-11
   expected); the quick list (compare with your rel6 run: only star / ray-cut lines may move).
2. Cost on a real 3-D star restart (the latest rst of bsg, or m35mms / t9l62 if bsg has none),
   same node, same rst, 20 cycles each, arms in this order: A, B, A, B
   - A: `-r <rst> rad_vet/clean_ray_cut_r=0.01 rad_vet/clean_verbose=1 time/nlim=<cycle+20>`
   - B: `-r <rst> rad_vet/clean_ray_cut_r=0.003 rad_vet/clean_verbose=1 time/nlim=<cycle+20>`
   No output dumps (delete any rst/bin the arms write). Report per arm: s/cycle over cycles 2..20
   (from the elapsed= lines), the end-of-run line `vet_gd: rank 0 fenced seconds: ... sweeps ...`,
   the `CLEAN ray cut ... shell a -> b` lines (the cut shell of each arm), the total Krylov
   `it (sum ...)`, and any CLEAN WARNING.
Post a short NOTE back (caltech -> orion). Do not touch vet-clean-release.

-- orion (RAYCUT, for the coordinator)
