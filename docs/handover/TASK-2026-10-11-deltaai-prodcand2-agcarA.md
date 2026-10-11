# TASK raven -> DeltaAI: AG Car A, FINAL production candidate prodcand2 (1 short ghx4-interactive job)

From the Raven session. Raven gpudev estimates this job (31071238) will start at 10:40, a wait of more than 20 min, so we are asking you as well.
**First result wins**: post your NOTE as soon as you have it.
If ghx4-interactive will not start within about 2 h, reply SKIP in the NOTE.

## Code (fork github.com/jing-ze-ma/athenak), PROBLEM=he_star_m1, `build_inc_deltaai.sh hes_gpu`
Use branch `prodcand2-1011` at **b89c137e**. Its code commit is **767d3313**; b89c137e only adds docs.
- It contains prodcand-1011 (e66ee669 + 558be1c2, which you ran in TASK ...-prodcand-agcarA).
- It also merges these branches:
  - predfix-1011 (`time2_pred_cap = 0.5` default; `vet_x1_periodic_inflow = source`; no fatal on an already-inadmissible hesdirk2 stage);
  - segfix-1011 (`vet_seg = mean | split | emit`, default mean);
  - icgen-1011 (opt-in `problem/he_ic_rad = sph1d`);
  - sptests-1011 (opt-in `implicit_hr_excess = iso`).
- No REF binary is needed.

Remember `export MPICH_GPU_SUPPORT_ENABLED=1`.

## Files
Use `docs/handover/agcarA_files/` on bsg-files-1009: the IC, the Rosseland and Planck tables, and agcar_A_va.athinput.
Point he_ic_file, he_opac_table and he_planck_table at your copies.
**AthenaK rejects a command-line key that the deck does not name** ("Parameter ... on command line not found"; this
happened on Raven). Make one deck copy for all arms, with these keys added at their default values:
- in `<rad_m1>`: `vet_seg = mean`, `implicit_hr_excess = none`, `ic_rad = fs_relax`;
- in `<problem>`: `he_ic_rad = file`, `he_ic_rad_gas = column`, `he_ic_rad_hse = false`.
With every key at its default, this deck runs exactly like agcar_A_va.athinput.

## Runs: one job, 1 node x 4 GH200, 1 rank per GPU (4 blocks 480x64x64), --time 0:40, plain srun
Common keys for all runs:
`time/nlim=60 time/ndiag=1 output5/dt=1e30 rad_m1/implicit_cfl=1251.3 output2/dcycle=60 output3/dcycle=60 output6/dcycle=60 output7/dcycle=60`
This gives fixed dt = 1239.996 s, hst every cycle, and dumps at cycles 0 and 60.

CAND = `rad_m1/time_scheme=hesdirk2 rad_m1/time2_vet_col=predict rad_m1/implicit_vimp=true rad_m1/time2_vstage=true rad_m1/implicit_blend_xthin=0 rad_m1/implicit_thin_switch=fse rad_m1/implicit_thin_corr=hc`
All other settings stay at their defaults. In a he_star_m1 build with no ic_rad key, that includes `ic_rad = fs_relax`.

| arm | keys |
|---|---|
| cand1 | CAND |
| cand2 | CAND (repeat; run-to-run check) |
| emit | CAND `rad_m1/vet_seg=emit` |
| ic | CAND `problem/he_ic_rad=sph1d problem/he_ic_rad_gas=re_atm problem/he_ic_rad_hse=true rad_m1/ic_rad=fs_relax` |
| iso | CAND `rad_m1/implicit_hr_excess=iso` |

## Report: NOTE-2026-10-11-deltaai-prodcand2-agcarA.md
- For every run:
  - rc, FATAL count and NON-CONVERGED;
  - Picard mean and max;
  - s/cycle as steady (excluding rotation cycles 10, 20, ... and the last cycle), rotation and average;
  - cost per simulated time (s/cycle / dt);
  - the init wall time. For ic, also the setup lines `he_ic_rad = sph1d ... N s` and `ic_rad = fs_relax ...`.
- L_top (hst) per cycle, cycles 0-60, for every arm, normalised by cand1's L_top at cycle 0.
- Is cand2 bitwise identical to cand1? Compare the bin data and the hst rows.
- cand1 vs prodcand e66ee669 cand1 (your previous NOTE): are they bitwise identical? Expected: yes, or differences only
  where `time2_pred_cap` limits. Report the run.log line with the count of limited cell-builds.
- E/J_fs and E/J_tw in the envelope at cycle 60 (median, p5-p95) in the bands 1.0-1.05, 1.05-1.5, 1.5-2 and 2-3 R_ph,
  with R_ph = 388.3 Rsun (vetacc ana_va.py `ej` mode), for every arm.
- emit, ic and iso vs cand1 at cycle 60: the shell-mean relative difference (max over shells in the band) of dens, velx,
  E and F1, in the bands < 0.8, 0.95-1.05 (photosphere) and > 1.05 (envelope) R_ph.
- GPU-h used.
