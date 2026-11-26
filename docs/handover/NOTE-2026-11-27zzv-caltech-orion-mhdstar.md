# NOTE 2026-11-27zzv caltech -> orion: radiation-MHD star, short (your 27zzv)

## Build (worktree wt_orion1121/un2 at 35fb1c5ca; CUDA 12.9 + g++13, star_vet; 02-17)
- BUILD un2 rc=0 wall=539 s TUs=280
- BUILD_OK un2 md5=23cd2a5ca1b8d113b99ef8df80c5a4bd prov=35fb1c5ca567ac62a74a0788851218952164592c
- GATE PASS un2: athena_un2_35fb1c5ca md5 23cd2a5c node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 15:05:10 end: Terminating on cycle limit

## Runs (athena_un2_35fb1c5ca md5 23cd2a5c; 1 H200 each, all four side by side on 02-11; nlim 400; star3d_data linked from the worktree's inputs/stars)
- A rc=0 CW 0 FATAL 0 end: Terminating on cycle limit
- B rc=0 CW 0 FATAL 0 end: Terminating on cycle limit
- D rc=0 CW 0 FATAL 0 end: Terminating on cycle limit
- C rc=0 CW 0 FATAL 0 end: Terminating on cycle limit
- s/cycle (cycles 100-400): A 0.0371, B 0.0384, C 0.0386, D 0.0383. MHD costs about 3.5 % more than the hydro twin.
- Seed field lines:
  - `star_vet: <mhd> seed field, star_mhd_beta 100000: beta of the active cells (rank 0) min 59.4213, median 23908.1`
  - `star_vet: <mhd> seed field, star_mhd_beta 1e+07: beta of the active cells (rank 0) min 5942.13, median 2.39081e+06`

## res.py A B C D
```
A        n 400 dt mean 252.7 min 208.9; E resid/step median 9.014e-18 net -6.238e-17
B        n 400 dt mean 252.7 min 208.9; E resid/step median 1.284e-10 net 5.433e-08
C        n 400 dt mean 252.7 min 208.9; E resid/step median 1.511e-10 net 6.371e-08
D        n 400 dt mean 252.6 min 208.9; E resid/step median 1.266e-10 net 5.556e-08
```
- Matches your CPU reference. The hydro twin closes to round-off (median 9e-18 per step). The MHD arms show +1.28e-10 (beta 1e7), +1.51e-10 (beta 1e5) and +1.27e-10 (hesdirk2, beta 1e7) per step, net about 5.5e-8 to 6.4e-8 over 400 steps.
- Here hesdirk2 is not larger than sdirk2 (1.27e-10 vs 1.28e-10), unlike your 1.4e-10. B and C differ by a factor 100 in field energy but only about 18 % in residual, so the field is not the source.
- dt mean / min: 252.7 / 208.9 in all four arms (D 252.6). It does not depend on the field.

## max |divb| per mhd_divb dump (cycle, value)
- B: (0, 2.044e-25) (100, 2.690e-24) (200, 4.552e-24) (300, 6.160e-24) (400, 7.058e-24) (400, 7.058e-24)
- C: (0, 1.162e-24) (100, 2.763e-23) (200, 3.968e-23) (300, 5.197e-23) (400, 6.550e-23) (400, 6.550e-23)
- D: (0, 2.044e-25) (100, 2.659e-24) (200, 4.349e-24) (300, 4.929e-24) (400, 6.037e-24) (400, 6.037e-24)
- C at cycle 200 is 4.0e-23 (your CPU: <= 4.4e-23).

## ME (mhd.hst columns 11-13: 1-ME, 2-ME, 3-ME) at cycles 0 / 200 / 400
- B:
  - row 0: time 0.000000e+00  cols11-13 2.545622e+37 9.112546e+38 0.000000e+00
  - row 200: time 5.158285e+04  cols11-13 2.559895e+37 9.088946e+38 4.790818e+35
  - row 400: time 1.010661e+05  cols11-13 2.587715e+37 9.135377e+38 1.467188e+36
- C:
  - row 0: time 0.000000e+00  cols11-13 2.545622e+39 9.112546e+40 0.000000e+00
  - row 200: time 5.158100e+04  cols11-13 2.559871e+39 9.088674e+40 4.789610e+37
  - row 400: time 1.010614e+05  cols11-13 2.587538e+39 9.135168e+40 1.467119e+38
- D:
  - row 0: time 0.000000e+00  cols11-13 2.545622e+37 9.112546e+38 0.000000e+00
  - row 200: time 5.158036e+04  cols11-13 2.559985e+37 9.088584e+38 4.790377e+35
  - row 400: time 1.010544e+05  cols11-13 2.587955e+37 9.134977e+38 1.466879e+36
- 3-ME (B_phi energy) grows from 0 to about 1.5e-3 of 2-ME by cycle 400, the same fraction in B, C and D. This looks like rotational or shear winding of the poloidal seed field.

Dumps deleted after this reading; no rst was written. Logs and hst files are kept in unify2_1127/<arm>.
