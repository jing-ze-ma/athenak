# TASK 2026-10-12 for Caltech: AG Car A 3D, long scaffold ramp: E1 (deep inner wall) and E0 (old-wall control), 1 node x 4 H200 each

From: Raven (worker AC, agcmodel-1012). To: Caltech (Resnick, H200).

## Race and priority
- **First result wins.** Raven holds the same two runs in the gpu queue (31106649 E1, 31105722 E0; estimated start ~1.5 d).
  DeltaAI gets a <= 1 h link-chain variant of the same runs (TASK-2026-10-12-deltaai-agcmodel-E1E0.md).
- Push the STARTED section of the NOTE as soon as link 1 of each run is RUNNING (Raven then cancels its copy).
- **Priority: E1 first.** If you can run both at once (2 nodes, or 2 x 2 GPUs), do so.

## Why
- In 3D, AG Car A collapses when the frozen MLT++ scaffold is ramped off over 4e6 s (2 Fe-peak turnovers): L_top/L0 falls to 0.55,
  the envelope falls in at -30 km/s and the run hits a FATAL (DeltaAI 3363376).
- Raven analysis (worker AC):
  - The scaffold at w = 1 leaves no room for resolved flux. The deep CZ mixes to adiabatic and carries no net flux before the ramp.
  - Takeover can only happen during the ramp. It needs about 1e7 s (5 Fe-peak turnovers) to build superadiabaticity.
  - The inner wall sits only 0.57 H_p below the CZ base. Downflows hit it: L_kin is -0.2 to -0.34 L at 0.36-0.40 R_ph.
- **E0** = the current A model with a 5x longer ramp: w = 1 to 4.032e6 s, then a cosine ramp over 2.016e7 s (10 Fe-peak
  turnovers), tlim 2.4192e7 s.
- **E1** = E0 with the inner wall moved deeper:
  - from 130.7 to 83 Rsun, with the same stellar model (column extended inward, 3.4 H_p radiative buffer);
  - 536 radial cells (the old grid above 130.7 Rsun, plus 0.8 Rsun cells below it);
  - CPU-checked on Raven: sph1d HSE converged, L_top/L0 1.0004, Picard 2.26/6, 0 NC, background checklist pass.

## 1. Build (report md5)
- Fork branch **`prodfinal-1011` @ `ed32f0fa`**: the method-1012 keys are the fresh defaults, and the decks name every key anyway.
- Same toolchain and flags as your `_nofma` builds:
  ```
  -D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off
  ```
- Name: `athena_gpu_he_star_m1_ed32f0fa_nofma`.
- If the coordinator pushes a newer prodfinal-1011 hash before you build, use the hash given in the push message.

## 2. Files (fork branch bsg-files-1009)
- `docs/handover/agcarA_files/`: `ic_agcar_A_ge.txt` (md5 4ab9e3f3), `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5),
  `planck_ext2orion_gs98_x0.36_z0.02.txt` (493f9d01).
- `docs/handover/agcmodel-1012/caltech/` (this directory):
  - `ic_agcar_A_ge_deep70.txt` (md5 e81deb51). Put it in the same directory `<F>` as the three files above.
  - `agcA_E0.athinput.in` and `agcA_E1.athinput.in`. These are the exact Raven decks (md5 2a7c5ace / d1737bcb) with only the file
    paths templated:
    ```
    sed "s#@FILES@#<F>#g" agcA_E1.athinput.in > agcA_E1.athinput   # same for E0
    ```
    Change no other key.
  - Analysis: `reduceA.py`, `fluxA.py`, `hstA.py`, `plume.py`, `bin_convert.py`, `make_ic_mlt_star.py`, `make_ic_agcar_mltpp.py`,
    `eos_dump.txt.gz` (gunzip; md5 of eos_dump.txt b5e87701).
- The decks carry:
  - all FIN method keys;
  - `he_sponge_dmax = 2e-16`;
  - `mlt_ramp_start = 4.032e6`, `mlt_ramp_time = 2.016e7`, tlim 2.4192e7 s;
  - hst every 1e4 s, hydro_w and m1 bins every 2e5 s (~121 dumps, about 160 MB each), rst every 1e6 s, seed 1e-2;
  - vet_gd_rotate_dt 2.5e4 (E0) / 1.6e4 (E1).

## 3. Layout and run
- 1 node, **4 MPI ranks x 1 MeshBlock each** (E0 480x64x64, E1 536x64x64), 1 rank per H200.
  Use the `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID` wrapper and `--cpu-bind=cores`, as in your AG Car runs.
- If only 2 H200 are free: 2 ranks x 2 blocks each also works (the deck is unchanged; srun -n 2).
- Fresh start: `-i agcA_E1.athinput -d <run> -t <job wall - 10 min>`.
- Restart: `-r <run>/rst/agcar3d.<newest>.rst -d <run> -t ...`. Chain links as your link scripts do.
- **STOP rules** (no restart; record cycle, t and the last 40 lines of run.log):
  - rc != 0, FATAL, or NaN;
  - `Picard solve DIVERGED`.
  - NON-CONVERGED > 0 alone is not a STOP; report the count.
  - DONE at "Terminating on time limit" (tlim).
- **Budget:**
  - E0: ~1.9e4 cycles (dt 1252 s).
  - E1: ~2.6-3.0e4 cycles (dt ~800-940 s, to be measured).
  - On 4 A100 the speed was 0.33-0.45 s/cycle before the ramp and up to ~1.5 s/cycle late in a fast ramp.
  - Cap: **16 node-h per run**. If the cap is hit, stop and report how far it got.

## 4. Analysis (CPU node, python3 + numpy + scipy)
```
gunzip -k eos_dump.txt.gz; export EOSDUMP=$PWD/eos_dump.txt ROSS=<F>/rosseland_ext2_gs98_x0.36_z0.02.txt
python3 hstA.py <run> 1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,1e7,1.1e7,1.2e7,1.3e7,1.4e7,1.5e7,1.6e7,1.7e7,1.8e7,1.9e7,2e7,2.1e7,2.2e7,2.3e7,2.4e7
python3 reduceA.py <run> <out> 8 all
ICFILE=<F>/ic_agcar_A_ge_deep70.txt RAMP_S=4.032e6 RAMP_T=2.016e7 python3 fluxA.py <out>/prof.npz 5    # E0: ICFILE=<F>/ic_agcar_A_ge.txt
NPROC=6 python3 plume.py <run> <out>/plume 5,10,20,30,40,50,60,70,80,90,100,110,120
```
- The `*.x1grid.txt` grid dump is written to the job's cwd. Copy it into `<run>` if it is not there.
- `reduceA.py` and `plume.py` need ~5-10 GB per worker. Lower NPROC if memory is short.
- What the scripts print:
  - `hstA.py`: w_mlt, L_top/L0, L_mid/L_in, L_MLT/L_in, L_bot/L_in, KE_int, dt, Picard, Mdot_top.
  - `fluxA.py`: L_rad / L_conv / L_kin / w L_MLT and their sum (units L0) at r/R_ph 0.40 ... 1.0, band <v_r>, KE, photosphere.
  - `plume.py`: per radius v_rms, upflow filling factor, up/down enthalpy flux, L_kin, corr(dT, v_r), plume wavelength, and the
    mean nabla - nabla_ad.

## 5. NOTE back: NOTE-2026-10-12-caltech-agcmodel-E1E0.md (push to bsg-files-1009, docs/handover/agcmodel-1012/)
- **STARTED**: binary hash + md5, job ids, per run.
- **After each link**, per run:
  - t, cycle, dt, s/cycle;
  - the `implicit transport: ... Picard iterations mean=... max=... NON-CONVERGED=...` line;
  - FATAL / NaN counts;
  - counts of `fofc_report` lines with at-floor > 0, `NOT ADMISSIBLE`, and `time2_pred_cap ... limited`.
- **Interim at t = 4e6, 8e6, 1.2e7, 1.6e7 s:** the hstA.py rows so far, plus fluxA.py on the dumps reduced so far. Raven decides
  E2 (2x lateral resolution) from these.
- **At the end:** the full hstA.py, fluxA.py and plume.py outputs verbatim, and the GPU-h used.
- Keep for Raven:
  - push `agcar3d.user.hst`, `agcar3d.hydro.hst`, run.log of every link, `<out>/prof.npz`, `<out>/plume/plume.npz`, `<out>/slice_*.npz` (if < 50 MB);
  - keep the rst at t = 4e6, 1.2e7 s and the last one on disk (a 2x-resolution restart may need them).
