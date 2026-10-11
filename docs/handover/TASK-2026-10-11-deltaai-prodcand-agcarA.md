# TASK raven -> DeltaAI: AG Car A, production-candidate gate (1 short ghx4-interactive job)

From the Raven session. Raven gpudev estimates a start of 07:32 for this job (>20 min wait), so we are asking you as well.
The same job is queued on Raven as 31068627. **First result wins**: post your NOTE as soon as you have it.
If ghx4-interactive will not start within ~2 h, reply SKIP in the NOTE.

## Code (fork github.com/jing-ze-ma/athenak), PROBLEM=he_star_m1, `build_inc_deltaai.sh hes_gpu` per commit
- CAND: branch `prodcand-1011` **e66ee669**. This is the production candidate. It is deepatm-1011 + vgdrot, merged with
  icrad-1011, cbfix-1011, moveq-1011 and sptests-1011. New defaults:
  - `vet_bc_diff=face`;
  - `ic_rad=fs_relax` on fresh starts;
  - `vet_source_noesrc_h2be=true`;
  - `implicit_hr_face=sf`;
  - `implicit_f0_init=auto`, `implicit_rad_inertia=true`, `implicit_nofb_work=true`.
- REF: branch `rt-integration` **d28dacd8**. This is the old production code.

## Files
Use `docs/handover/agcarA_files/` on this branch (bsg-files-1009). It holds the IC, the Rosseland and Planck tables, and
agcar_A_va.athinput. Point he_ic_file, he_opac_table and he_planck_table at your copies.
Remember `export MPICH_GPU_SUPPORT_ENABLED=1`.

## Runs: one job, 1 node x 4 GH200, 1 rank per GPU (4 blocks 480x64x64), --time 0:40, plain srun
Common keys for all runs:
`time/nlim=60 time/ndiag=1 output5/dt=1e30 rad_m1/implicit_cfl=1251.3 output2/dcycle=60 output3/dcycle=60 output6/dcycle=60 output7/dcycle=60`
This gives fixed dt = 1239.996 s, hst every cycle, and dumps at cycles 0 and 60.
1. cand1, CAND binary:
   `rad_m1/time_scheme=hesdirk2 rad_m1/time2_vet_col=predict rad_m1/implicit_vimp=true rad_m1/time2_vstage=true rad_m1/implicit_blend_xthin=0 rad_m1/implicit_thin_switch=fse rad_m1/implicit_thin_corr=hc`.
   All other settings are defaults (face, fs_relax and the vgdrot speed-up).
2. cand2: the same as cand1, repeated (run-to-run check).
3. ref, REF binary: the input's production keys only (be, hr with `implicit_blend_xthin=30`).

## Report: NOTE-2026-10-11-deltaai-prodcand-agcarA.md
- Every run:
  - rc, FATAL count and NON-CONVERGED;
  - Picard mean and max;
  - s/cycle as steady (excluding rotation cycles 10, 20, ... and the last cycle), rotation and average;
  - cost per simulated time (s/cycle / dt);
  - the wall time of the ic_rad start for cand.
- L_top (hst) per cycle for cycles 0-60 for cand1 and ref. For cand it should be flat from cycle 0.
- Is cand2 bitwise identical to cand1? Compare the bin data and the hst rows.
- E/J_fs and E/J_tw in the envelope at cycle 60 (median, p5-p95) in the bands 1.0-1.05, 1.05-1.5, 1.5-2 and 2-3 R_ph,
  with R_ph = 388.3 Rsun. The vetacc ana_va.py `ej` mode computes these.
- cand1 vs ref at cycle 60: the shell-mean relative difference (max over shells in the band) of dens, velx, E and F1.
  Use the bands < 0.8, 0.95-1.05 (photosphere) and > 1.05 (envelope) R_ph.
- GPU-h used.
