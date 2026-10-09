# TASK viper -> Caltech: test battery for the half-range VET xthin fix (Fix A), CPU + 2 H200 (user 10-09)

Uses the 2 H200 freed by stopping the He giant fresh run, plus CPU nodes. Do NOT touch BSG.

## What is being tested

`implicit_blend_xthin = 30` (the transparency override of the half-range implicit VET face flux) forces upwind on
every face where X = (c dt/dx)/(1 + c dt chi) >> 30. viper's ORDER study found it breaks convergence: grey atmosphere
vs Hopf error grows 1.5e-3 -> 4.2e-3 from 64 to 512 cells; thin pulses lose space order (0.69 / -0.82). Without xthin
(xthin 0) the scheme is 2nd order and equals central on smooth problems, but BiCGStab stagnated on beam tests.
Fix A (branch `xthinfix-1009`, commit `f1d355e8d`): new key `implicit_blend_xthin_mode = all | beam` (default `all` =
old behaviour, bitwise) and `implicit_blend_xthin_wmin` (default 0). `beam` applies the override only where the
beam-aware weight w > wmin, so smooth/diffusive faces stay central. NOTE: a command-line key that is absent from a
fresh input file is FATAL, so the scripts write the new keys into the inputs.

## Steps

1. Build from fork branch `xthinfix-1009` commit `f1d355e8d` (= vimphr-1009 6671987e + the fix; src/rad_m1 only):
   - CPU, no PROBLEM (built-in pgens), MPI on, Release: `ATHENA_CPU`.
   - H200 he_star_m1 exactly as your AG Car / BSG binaries (CUDA 12.9, HOPPER90, MPI, host -ffp-contract=off),
     incremental from your 98835d99 he build dir if possible: `ATHENA_GPU`. Report both md5.
2. CPU battery (correctness only, no timings): from a checkout of this branch (`bsg-files-1009`), with `REPO` = a
   checkout of `xthinfix-1009` (needs `tests_m1/` and `vis/python/`):

       B=<checkout>/docs/handover/xthinfix-1009
       REPO=<xthinfix-1009 checkout> ATHENA_CPU=<cpu binary> XTF_RUN=<run dir> P=<cores> bash $B/battery_cpu.sh all

   Arms: `cen` (central), `hr` (production, xthin 30 mode all), `hrb` (mode beam = the fix), `hrx0` (xthin 0).
   - order: pulses kappa 0.128 ... 12800 (tau_cell 2.5e-4 ... 400), grey atmosphere (steady, vs Hopf) 32-512,
     coupled radiation-acoustic wave P_rad/P_gas 100, tau_lambda 10 and 1e3; be and hesdirk2; X / T / C modes.
     ~1000 one-core runs, the finest 512 runs dominate. Split `order`, `gates`, `beams`, `ana` across jobs if
     useful (`bash battery_cpu.sh order` etc.; `ana` last).
   - gates: G1 (atm2d, cfl 1e2 / 1e4), G5 (thick pulse kappa 12.8-12800), G3 (Marshak cfl 1 / 10), arms hr / hrb.
   - beams: shd3b, cyl, xb20, ba0, ba20 (c dt/dx ~300, 100 cycles, np 1), arms hr / hrb / cen.
   - Output: `$XTF_RUN/RESULTS/{order_eval,atm_hopf,arm_minus_cen,gates,beams}.txt`.
3. GPU (2 H200, one job, ~15-25 min): AG Car A production input from your `agcarb_1009/files`
   (`agcar_rcxA_ge_accel_st_mg_hr.athinput`), t = 0, 30 cycles, arms all / beam / all / beam on 2 GPUs, then
   beam / all on 1 GPU (scaling point; 4 blocks on one rank, ~80 GB as your link_agcarb_1g.sh):

       ATHENA_GPU=... IN_A=.../agcar_rcxA_ge_accel_st_mg_hr.athinput XTF_RUN=... REPO=... BUNDLE=$B \
         sbatch --export=ALL [your 2-H200 header] $B/gpu_agcar.sh

   Adapt only the #SBATCH header and the srun / CUDA_VISIBLE_DEVICES line to your link scripts. Smoke first
   (`time/nlim=3` by editing O= for one arm) per the smoke rule. Output `$XTF_RUN/gpu.<job>/RESULTS_gpu.txt`.

## Expected (viper / earlier numbers; what counts as PASS)

| test | PASS for `hrb` |
| --- | --- |
| atm (G1v) Hopf L1 32 ... 512 | decreasing ~2nd order, 512 near hrx0's 3.7e-4 (hr: 5.5e-3 1.5e-3 1.9e-3 3.1e-3 4.2e-3) |
| pulses X / C orders | equal to cen / hrx0 (cen hesdirk2 X: k0.128 1.75 1.94 2.04; k12.8 1.90 2.25 2.24); hr breaks at k0.128 (0.69) and k12.8 (1.07) |
| pulse T orders | as cen (hesdirk2 + vet_sc is 1st order in time where the tensor changes; known, not this fix) |
| rw_t10 | E converges (hr: E X -1.67 -1.18, no convergence) |
| `arm_minus_cen` for hrb | ~0 on the pulses (hrx0 - cen = 0.00 exactly) |
| G1 / G5 / G3 | G1 max dE/E <= 2e-2; G5 within 2 % of 1 at kappa 1280 / 12800; G3 L1 <= 0.02 |
| beams | no FATAL, NON-CONVERGED 0 in 100 cycles; inner BiCGStab mean comparable to hr (hr: xb20 ~37); cyl direction median <= 5 deg; shd3b umbra depth ~0, edges ~0.23 / 0.35 |
| AG Car A GPU | rc 0, NON-CONVERGED 0; all vs beam differ only above ~1 R_ph (interior bands ~0 or round-off); s/cycle beam vs all within noise; repeats bitwise |

## Report

Push `docs/handover/NOTE-2026-10-09-caltech-xthinfix-tests.md` to this branch with: binary md5s, and the RAW
contents of the five RESULTS/*.txt files and RESULTS_gpu.txt (no summarising beyond a 5-line verdict), plus any
FATAL / NON-CONVERGED. Run trees stay on Caltech.
