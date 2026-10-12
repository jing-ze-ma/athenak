# TASK raven -> DeltaAI: AG Car A full-size 3D run through the frozen-MLT++ scaffold ramp (fresh, 1 node x 4 GH200)

From: Raven (worker P, rampA-1011).

## Why
- Worker J (vimpbc-1011) found that AG Car B in full-size 3D collapses when the frozen MLT scaffold ramps off:
  - L_top/L0 falls 0.93 -> 0.47;
  - the envelope falls in at -63 km/s;
  - convective KE grows 5x but carries no flux.
- AG Car A in quasi-1D collapses the same way at its own ramp (L_top/L0 0.41 at w = 0).
- A's ramp has never been run in 3D. A production starts on Raven later today.

## Race and timing
- **First result wins.** Raven runs the same test as a gpudev 13-min chain (/raven/ptmp/jinma/rampA_1011/run). That copy
  is at t = 2.1e6 s and its chunks are interleaved with other users, so it finishes in roughly 4-7 h.
- Post **STARTED** (with the job id) in the NOTE as soon as your job runs.
- If you cannot start within ~2 h: SKIP and say so in the NOTE.

## Environment
- Same as your earlier AG Car / BSG tasks: `export MPICH_GPU_SUPPORT_ENABLED=1`, plain srun, 1 rank per GPU.
- 1 node x 4 GH200, 4 ranks. The deck has 4 MeshBlocks of 480x64x64, so 1 block per GPU.

## 1. Build (report the md5)
- Fork branch **`prodfinal-1011` @ `0197f020`**, built with `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1, GH200).
- On Raven, 0197f020 is bitwise equal to 1fbe8621 on this deck (A100).
- If the tip has moved, build 0197f020 exactly and do not use the tip.

## 2. Files
- From `bsg-files-1009`:
  - `docs/handover/agcarA_files/`: `ic_agcar_A_ge.txt` (md5 4ab9e3f3), `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5),
    `planck_ext2orion_gs98_x0.36_z0.02.txt` (493f9d01). Put all three in one dir `<F>`.
  - `docs/handover/rampA-1011/` (this TASK's dir):
    - `agcar_A_ramp.athinput.in`: the CURRENT Raven A production deck (prodfinal_1011/run/agcA/agcar_A_fresh.athinput v4,
      md5 18068a70). Differences from it: only the three file paths are templated (`@FILES@`), and the hydro_w and m1
      bins are written every 5e5 s instead of 1e5 s.
    - Key settings it carries:
      - `implicit_gas_src = matched`, `implicit_vimp_bc = true`, `he_sponge_dmax = 2e-16`, `vet_gd_halo_ipc = true`;
      - tlim 8.064e6 s;
      - hst every 1e4 s, rst every 5e5 s, log every 1e5 s;
      - `mlt_ramp_start = 4.032e6`, `mlt_ramp_time = 4.032e6` (cosine w: 1 -> 0 over 4.032e6 -> 8.064e6 s).
    - Analysis: `reduceA.py`, `fluxA.py`, `hstA.py`, `bin_convert.py`, `eos_dump.txt.gz` (gunzip it; md5 of
      eos_dump.txt b5e87701).
- Fill the deck:
  ```
  sed "s#@FILES@#<F>#g" agcar_A_ramp.athinput.in > agcar_A_ramp.athinput
  ```
  Do not change any other key.

## 3. Run (fresh from t = 0)
- Expected speed: A100 0.45 s/cycle, dt 1252 s (stable to 2.1e6 s on Raven). About 6450 cycles to tlim, so about 48 min
  on 4 A100. GH200 is probably faster, but dt may shrink during the ramp.
- **Preferred:** one job to tlim 8.064e6 s, if your wall limit allows ~1.5 h:
  ```
  srun -n 4 <bin> -i agcar_A_ramp.athinput -d <run> -t <job wall - 5 min>
  ```
  `-t` makes the code stop cleanly and write a final rst.
- **Otherwise:** at most 2 chained jobs of <= 1 h each, reaching at least t = 6.5e6 s (tlim if possible). Restart with:
  ```
  srun -n 4 <bin> -r <run>/rst/agcar3d.<newest>.rst -d <run> -t <...>
  ```
  The code writes a final rst at the `-t` wall limit.
- Stop at the first FATAL or NaN and report it.

## 4. Analysis (on the node after the run, or on a CPU node; python3 with numpy + scipy)
```
gunzip -k eos_dump.txt.gz
python3 hstA.py <run> 5e5,1e6,1.5e6,2e6,2.5e6,3e6,3.5e6,4e6,4.25e6,4.5e6,4.75e6,5e6,5.25e6,5.5e6,5.75e6,6e6,6.25e6,6.5e6,6.75e6,7e6,7.25e6,7.5e6,7.75e6,8e6,8.064e6
EOSDUMP=$PWD/eos_dump.txt ROSS=<F>/rosseland_ext2_gs98_x0.36_z0.02.txt python3 reduceA.py <run> <out> 8 all
ICFILE=<F>/ic_agcar_A_ge.txt python3 fluxA.py <out>/prof.npz 1
```
- `reduceA.py` reads hydro_w + m1 bins pairwise. Each worker needs ~5-10 GB; lower NPROC if memory is short.
- What the scripts print:
  - `hstA.py`: w_mlt, L_top/L0, L_mid/L_in, L_MLT/L_in, L_bot/L_in, KE_int, dt, Picard, Mdot_top at the listed times.
  - `fluxA.py`, per bin dump:
    - L_rad (comoving M1), L_conv (resolved fluctuation enthalpy flux, gas + radiation), L_kin, L_MLT = w x frozen IC
      MLT flux, and their sum, all in units of L0 at r/R_ph = 0.55 ... 1.0;
    - band max <v_r> [km/s];
    - KE below R_ph;
    - photosphere radius (per-column tau = 2/3, mean / max).

## 5. Report: NOTE-2026-10-12-deltaai-agcarA-ramp.md (push to bsg-files-1009, docs/handover/rampA-1011/)
- A **STARTED** line with the job id.
- Build: hash and md5. Per job: rc, FATAL / NaN counts, wall time, s/cycle, zone-cycles/s.
- The run.log line `implicit transport: solves=... Picard iterations mean=... max=... NON-CONVERGED=...`. Also count:
  - `fofc_report` lines with at-floor > 0;
  - `NOT ADMISSIBLE` / stage-fallback lines;
  - `time2_pred_cap ... limited` lines.
- The full `hstA.py` table and the full `fluxA.py` output. Paste them verbatim.
- Keep for Raven: `agcar3d.user.hst`, `agcar3d.hydro.hst`, `run.log`, `<out>/prof.npz`, and `<out>/slice_*.npz`.
  Push them if the total is < 50 MB; otherwise push prof.npz + hst + log only.
- GPU-h used.
