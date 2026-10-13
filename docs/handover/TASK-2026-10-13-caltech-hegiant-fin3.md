# TASK raven -> Caltech: He giant production on FIN3 (dic-prod-1013 @ dad00933) -- smoke + restart check first, chain HELD/gated

From the Raven coordinator (worker HG, hegiant-1013). USER asked: run the He giant (MESA r15140 `M2pt754_Porb100_profile29`, the
model of S. Wu; IC column `ic_giant_own.txt` md5 ddc72692) as production on the final FIN3 method, on Raven, Viper AND Caltech,
same deck physics and same binary hash everywhere. This REPLACES the 10-09 He giant fresh N897 run
(TASK-2026-10-09-caltech-hegiant-fresh.md, chain 4288960-65, stopped at 1.30 d on the old xthin-30 scheme): do NOT restart that
run dir; start a NEW dir from t = 0.
**Do not touch the BSG chain 4434961-63** or any other job. This run is extra; it may wait for free H200s.

## 1. Build (report md5)
- Fork branch **`dic-prod-1013` @ `dad00933`** = rotavg-rst2-1013 9972a56d (your BSG production hash, md5 4235b7dc build) + DI's
  he_ic_equil (17a22820, dad00933). Key-off it is BITWISE to 9972a56d (Raven/Viper checked). The He giant deck names
  he_ic_equil = false (OFF, see section 2), but all three sites must run the SAME hash, and AG Car E1 production on Raven/Viper
  runs dad00933 -> build dad00933.
- Same toolchain and flags as your 9972a56d build: `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On
  -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`, default `Athena_M1_VGD_NMAX`.
- Name `athena_gpu_he_star_m1_dad00933_nofma`. An incremental build from your 9972a56d tree is fine (2 commits, he_star_m1.cpp +
  driver.cpp).

## 2. Files and deck
- IC + tables: the SAME files as 10-09 (`docs/handover/hegiant-fresh-1009/`, MD5SUMS): ic_giant_own.txt ddc72692,
  rosseland_tops_hegiant_blend.txt c1ef063c, planck_tops_hegiant_blend.txt b699c60d. Reuse your 10-09 files dir.
- Deck: `docs/handover/hegiant-1013/hegiant_fin3.athinput.in` -> `sed "s#@FILES@#<your files dir>#g"`. Change no key.
  (= Viper's deck byte-for-byte apart from the 3 paths; Raven's adds only vet_gd_twin_lowmem = true, a bitwise memory key.)
- What changed vs the 10-09 deck (full list in `docs/handover/hegiant-1013/CHANGES.md`):
  - <rad_m1> = the AG Car E1 FINAL block (FIN3: hesdirk2, xthin 0 + fse/hc, vimp + vimp_bc, conservative2, vet_rotate_avg 4,
    vet_bc_inner auto, opac_newton true + slope_max 3, every key named), with the He giant's implicit_flux_x1min and
    vet_gd_rotate_dt 4.0e2 (~20 dt).
  - IC: ic_rad = fs_relax (E1), he_ic_rad = file (E1's sph1d atmosphere does NOT apply to this star: FATAL in its HSE loop, the
    frozen MLT scaffold carries 78-99 % of F up to R_ph), he_ic_equil = false (DI's equilibrium IC was tested: its re-marched
    atmosphere falls back, 1.03-1.10 R_ph density change 0.41 vs 0.037 without it; the coordinator may still switch it on).
  - floor tfloor 4000 -> 3000 K (below T_rad at x1max, the E1 rule) with eos_logt_min 3.5 -> 3.4; sponge r0 199 -> 180 Rsun
    (0.9 x1max, E1 rule); fofc_report diagnostics.
  - Unchanged: grid N897, EOS, IC column, tables, gravity (he_gm_column), wall, he_base_heat, scaffold ramp 5 -> 10 d, seed 1e-2,
    outflow top, cfl 0.3, tlim 2.592e6, outputs (hst 25 s, bins 0.25 d, rst 0.5 d).

## 3. Layout
1 node, **2 MPI ranks x 2 MeshBlocks (897x64x64), 1 rank per H200** (`--ntasks-per-node=2 --gres=gpu:h200:2`, your
CUDA_VISIBLE_DEVICES=$SLURM_LOCALID wrapper, `--cpu-bind=cores`). Do NOT use smaller lateral blocks: 64x32 FATALs in vet_gd
("clamped lateral reads"). Expected memory: Raven A100 1 block/GPU: >40 GB (OOM) without lowmem, 31.3 GB with vet_gd_twin_lowmem; old method was 36.6 GB. On H200 with 2 blocks/GPU and lowmem false expect ~75-95 GB (ESTIMATE; 141 GB fits). If it does not fit, add `vet_gd_twin_lowmem = true` (bitwise, memory only) and say so.

## 4. Smoke + restart check (required before any production link)
- A: fresh, `time/nlim=24 time/ndiag=1 output1/dt=1e-6 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output6/dt=1e30
  output5/dt=240` (rst 00001 at t = 240 s, ~cycle 11). Check: rc 0, no FATAL / NaN / DIVERGED / NON-CONVERGED; the log lines
  "he_ic_equil ... done", "vet_rotate_avg = 4"; dt, Picard mean, s/cycle (cycles 3-24), GPU peak.
  Raven/Viper reference: Raven gpudev 31132998 (4 A100, 1 block/GPU, lowmem): rc 0, NC 0, Picard mean 5.16 max 10, dt 20.0 s every cycle, 1.41 s/cycle, init 72 s
- B: NEW srun `-r <A>/rst/hegiant.00001.rst time/nlim=24` into a fresh dir. The cycle-24 bin dumps (every variable, byte for
  byte after the header) and the user.hst rows after the restart point must be BITWISE vs A. Raven A100: BITWISE (31132998: 12/12 user.hst + hydro.hst rows, all final hydro_w + m1 dumps).
- If not bitwise or not clean: STOP and report; do not start production.

## 5. Production (only after section 4 passes; queue the chain HELD or afterok on the smoke)
- Run dir e.g. `/resnick/groups/carnegie_poc/jingze/hegiant_1013/fin3`; fresh t = 0 to tlim 2.592e6 (30 d).
- Your BSG link script pattern (run.cfg read at job start): FRESH if no rst, else the newest VALID rst with the rst_info.py
  last_time keys, XKEYS empty. STOP on rc != 0 / FATAL / NaN / DIVERGED / NON-CONVERGED > 0; DONE at rst t >= 2.592e6.
- rst pruning keep 3 (section 4b of TASK-2026-10-13-caltech-bsg-fin3.md, your installed prune_rst.sh, same EXIT trap).
  A He giant rst is about 10 GB (cycle-0 file 7 GB); keep 3 = ~30 GB.
- Links: 12 h (as your BSG), afterany. Estimated 5-6 links (~1.4 (A100 x4; your old-method H200 x2 number was 1.11) s/cycle x ~1.3e5 cycles).
- Start rule: first link afterok on the smoke/restart job, or submit held (`--hold`) and release after the NOTE says PASS.

## 6. Report: NOTE-2026-10-13-caltech-hegiant-fin3.md (push to bsg-files-1009)
- build md5; smoke numbers; restart verdict; chain job ids; then per link t, cycle, dt, s/cycle, Picard, counts, newest rst.
