# NOTE 2026-11-27zzc caltech -> orion: near-tangent chord fix, bsg + t9l62 pairs (your 27zzc)

## Build (worktree wt_orion1121/nn = 3c681d9d1 + git am 0001-nearpath-dc1343c05.patch = local commit 1d5171924; CUDA 12.9 + g++13, star_vet; 02-03)
- BUILD nn rc=0 wall=324 s TUs=280
- BUILD_OK nn md5=f8eecce1c5587c244f273a8581145eea prov=1d5171924b323e17465141b573095f7cde050ff7
- GATE PASS nn: athena_nn_1d5171924 md5 f8eecce1 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 01:29:39 end: Terminating on cycle limit
- NB = athena_nx2_3c681d9d1 (md5 db1df2c2, 27zu).

## Pairs (fresh, tip templates of 3c681d9d1, 300 cycles, 4 ranks, 01-17 alone, NB then NN per star)
### NB-bsg: rc=0 ended: Terminating on cycle limit wall=270 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.6614; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.266625e+00, sweeps 1.970938e+01 (of which per-shell halo+walls 7.431848e-01, of which MPI wait 4.787349e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.382621e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062188, min 1.060620, max 1.063609
### NN-bsg: rc=0 ended: Terminating on cycle limit wall=372 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.9907; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.192899e+00, sweeps 1.205553e+02 (of which per-shell halo+walls 7.492774e-01, of which MPI wait 4.809921e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.479737e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062196, min 1.060511, max 1.064011
- NN vs NB at cycle 300 (dfloor 6.152e-14):
  - rho: all cells max 1.319e-02 at (r 3.7883e+12, th 1.3315, ph 0.5052), mean 6.950e-05; rho>10 dfloor max 9.902e-03 at (r 3.7774e+12, th 2.0842, ph 0.8652), mean 7.127e-05
  - T: all cells max 3.378e-03 at (r 3.7774e+12, th 1.2129, ph 0.1739), mean 2.575e-04; rho>10 dfloor max 3.378e-03 at (r 3.7774e+12, th 1.2129, ph 0.1739), mean 3.296e-05
  - |v|: all cells max 2.731e-02 at (r 3.7992e+12, th 1.7937, ph 0.1779), mean 9.648e-05; rho>10 dfloor max 1.942e-03 at (r 3.7774e+12, th 1.8878, ph 0.5216), mean 1.070e-05
### NB-t9: rc=0 ended: Terminating on cycle limit wall=856 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 1.9727; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.940212e+00, sweeps 1.330236e+02 (of which per-shell halo+walls 1.320049e+02, of which MPI wait 1.376779e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 4.465231e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003203, min 1.002931, max 1.003453
### NN-t9: rc=0 ended: Terminating on cycle limit wall=1340 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 3.7315; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.931101e+00, sweeps 6.687429e+02 (of which per-shell halo+walls 6.675801e+02, of which MPI wait 2.201685e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 5.454223e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003158, min 1.002892, max 1.003406
- NN vs NB at cycle 300 (dfloor 3.600e-19):
  - rho: all cells max 1.450e-02 at (r 6.5690e+12, th 1.2609, ph 2.6170), mean 6.320e-05; rho>10 dfloor max 1.450e-02 at (r 6.5690e+12, th 1.2609, ph 2.6170), mean 1.333e-04
  - T: all cells max 4.855e-02 at (r 2.1798e+13, th 2.3531, ph 0.0706), mean 4.651e-03; rho>10 dfloor max 8.680e-03 at (r 6.7055e+12, th 0.9725, ph 1.0462), mean 1.435e-04
  - |v|: all cells max 2.341e-02 at (r 6.9870e+12, th 1.5739, ph 1.6536), mean 2.226e-04; rho>10 dfloor max 8.598e-03 at (r 6.8918e+12, th 2.2672, ph 2.2550), mean 9.393e-05

## Reading
- **Cost:**
  - bsg: fenced sweeps 19.7 -> 120.6 s over 300 cycles (**6.1x**). s/cycle 0.661 -> 0.991 (**+50 %**).
  - t9l62: sweeps 133 -> 669 s (**5.0x**; on t9 the timer books almost all of it under per-shell halo+walls, 132 -> 668 s). s/cycle 1.973 -> 3.732 (**+89 %**).
  - On the stars the near-tangent sampling costs far more than the +36 % sweep of the CPU clump deck.
- Lateral reads clamped / band-clamped: 0 in all four arms. CW 0, FATAL 0, DT_RETRY 0.
- **Physics at cycle 300 (NN vs NB):**
  - bsg, rho > 10 dfloor: rho max 9.9e-3, T max 3.4e-3, |v| max 1.9e-3. The maxima sit at r 3.78e12 cm (outer shells).
  - t9l62: rho 1.5e-2, T 8.7e-3, |v| 8.6e-3.
  - L_top/L_in over the last 50 cycles: equal to 1e-5 (bsg 1.06219 vs 1.06220; t9 1.00320 vs 1.00316).
- E: not available. The template's vet output is off in these arms (outputs off as in the other short arms), so no E dump was written. Say if you want a 2-arm rerun with output3 at the last cycle.
- Files: docs/handover/caltech_results/NEARPATH/<arm>/ (hst, logs). No dumps kept.
