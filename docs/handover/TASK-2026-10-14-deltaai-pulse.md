# TASK raven -> DeltaAI: E1 pulsation diagnosis (worker PU, pulse-1014): 4 REDUCED 3-D E1 arms, 1 node x 4 GH200, <= 3 links of <= 2 h

From Raven worker PU (coordinator-delegated investigation; no production). SHORT TEST class: one ghx4-interactive job at a time (2 h),
at most 3 links in total. If it cannot start within ~2 h, or the user says skip, post SKIP in the NOTE and stop.

## Why
Caltech's full E1 diverged at 1.24e7 s after a growing outer-envelope pulsation (inflation to 1.2 R_ph, infall -30..-80 km/s).
Question: physical (strange-mode / radiation-dominated envelope) or a setup artefact (seed, MLT scaffold, resolution)?
Raven runs the same arms at nx1 = 268 (gpudev, 15-min links); DeltaAI gives the nx1 = 536 (production radial) partner.

## 1. Build (report md5 + HEAD)
- Fork branch `fin4-dtretry` @ **5b7e1bc5** (053bd996 is the same physics, only DT_RETRY log wording; either is fine, say which).
  `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1), as your earlier he_star_m1 builds.

## 2. Files (`docs/handover/pulse-1014/`)
- `agcA_E1_16b0bba7.athinput.in`: the E1 production deck with @FILES@ paths. Resolve @FILES@ to a dir holding
  `ic_agcar_A_ge_deep70.txt` (agcmodel-1012/deltaai/, md5 e81deb51), `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5) and
  `planck_ext2orion_gs98_x0.36_z0.02.txt` (493f9d01) (both in agcarA_files/). Change no other key.
- `pulse_dai.sh`: one link = 4 independent single-rank arms, 1 GH200 each (536 x 32 x 32 = one MeshBlock):
  A536 (prod), S0_536 (`he_seed=0`), S3_536 (`he_seed=1e-3`), F536 (`mlt_ramp_start=1e6 mlt_ramp_time=2e6`).
  Overrides are on the command line (lateral 32x32, dt_retry on, tlim 1.25e7, m1 bins + rst every 1e6 s). Adapt LAUNCH to your srun
  (one process per GPU; 4 concurrent). `sbatch --export=ALL,BIN=..,DECK=..,S=<scratch>,WALL=01:52:00 <flags> pulse_dai.sh`.
- Next link: the same command; every arm restarts from its newest rst. Stop after link 3 or when all arms reach 1.25e7 s.
- STOP an arm only on rc != 0 / FATAL. `DT_RETRY` lines are expected (dt halving on a failed Picard solve) -- not a stop.

## 3. Report: NOTE-2026-10-14-deltaai-pulse.md (push to bsg-files-1009)
- Per link and arm: job id, t reached, cycle, s/cycle, number of DT_RETRY lines, the last Picard line, FATAL/NaN counts.
- After each link (so Raven can analyse early): push per arm `agcar3d.user.hst`, `agcar3d.hydro.hst` and the output of
  `python3 pulse_ana.py run <arm dir> <arm>_pulse.npz > <arm>_pulse.txt` (needs bin_convert.py beside it) into
  `docs/handover/pulse-1014/deltaai_out/`. Keep the scratch dirs.
