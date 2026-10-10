# TASK viper -> Caltech: xthinfix-1009 STEEP battery, CPU nodes + 2 H200 (user 10-10)

The viper login node is saturated, so this battery moves to Caltech. Run it on as many CPU cores as you have; first
result wins. Do not touch BSG / AG Car production.

## What is tested

The half-range implicit VET face flux of rt-integration has a transparency override `implicit_blend_xthin = 30` that
breaks convergence (XTHINFIX.md on viper). Branch `xthinfix-1009` adds:
- `implicit_blend_xthin_mode = all | beam | beam_kn | steep`, default `all` = old behaviour (bitwise). `steep` puts the
  gate X^8/(X^8 + X0^8) on every face, X = (c dt/dx)/(1 + c dt chi), X0 = `implicit_blend_xthin`.
- `implicit_hr_recon = dc | plm`, default dc. plm uses van Leer face states with a frozen limiter, outside the 7-point
  row.
- `implicit_hr_recon_fresh` and `implicit_hr_damp`.

A command-line key that is absent from a FRESH input is FATAL. The bundled inputs and scripts therefore carry the keys.

## Steps

1. **Build fork branch `xthinfix-1009` at commit `21795e2a`.** It already contains the vgdspeed merge (rt-integration
   5b304cf9).
   - CPU, no PROBLEM (built-in pgens), MPI on, Release: `ATHENA_CPU`.
   - H200 he_star_m1, exactly as your AG Car binaries (CUDA 12.9, HOPPER90, MPI, host -ffp-contract=off): `ATHENA_GPU`.
   - Also build rt-integration `5b304cf9` he_star_m1 H200 (or reuse your vgdspeed build if it is that commit):
     `ATHENA_BASE`.
   - Report all md5s.
2. **CPU battery.** Scripts are in `docs/handover/xthinfix-1009/steep/` on this branch. `REPO` = a checkout of
   xthinfix-1009 21795e2a (needs `tests_m1/`, `vis/python/`).

       B=<checkout>/docs/handover/xthinfix-1009/steep
       REPO=... ATHENA_CPU=... XTF_RUN=<run dir> P=<cores> bash $B/battery_steep.sh all

   You can split it as `order` | `gates` | `beams` (separate jobs or nodes), then `ana`.

   | arm | meaning |
   | --- | --- |
   | `cen` | central |
   | `hr` | production, xthin 30, mode all |
   | `hrx0` | xthin 0 (reference) |
   | `hrs` / `hrs15` / `hrs60` | steep gate with X0 = 30 / 15 / 60, dc |
   | `hrsp` | steep, X0 30, + plm |
   | `hrp` | mode all + plm |

   - **order:** pulses kappa 0.128 ... 12800, grey atmosphere vs Hopf 32-512, rw_t10 / rw_t1000; be and hesdirk2; space
     (X), time (T) and combined (C).
   - **gates:** G1, G5, G3.
   - **beams:** shd3b, cyl, xb20, ba0, ba20 (100 cycles, np 1). Report wall, Picard passes, BiCGStab inner iterations,
     NON-CONVERGED and plm positivity fallbacks.
   - Output goes to `$XTF_RUN/RESULTS/*.txt`.
   - Wall times on CPU are context only; cost is judged on the GPU point.
3. **GPU point.** `steep/gpu_steep.sh` on 2 H200, one job of about 20 min. Adapt only the #SBATCH header and the srun
   line.
   - Input: the AG Car A production input from your `agcarb_1009/files` (`agcar_rcxA_ge_accel_st_mg_hr.athinput`).
   - Arms base / all / steep / steep+plm / base / steep, 30 cycles each from t = 0.
   - It writes s/cycle, Picard, inner iterations, the bin-data identity of base vs all, and all vs steep by radius.
   - Smoke first (one arm with nlim 3) per the smoke rule.

## Viper partial results (stopped when moved; the reference for your numbers)

| arm | beams xb20 / ba0 / ba20 / cyl / shd3b: Picard mean, NON-CONVERGED |
| --- | --- |
| `hrs` / `hrs15` / `hrs60` (dc) | 5.3-5.8 / 4.7 / 5.7-5.8 / 4.97 / 2.97; all 0 (mode all: 4.0 / 4.75 / 5.44 / 5.0 / 3.0) |
| `hrsp` (steep + plm) | xb20, ba0, ba20 DIVERGED (Picard resid 1e27, FATAL); cyl 13, shd3b 13.1 |

Grey atmosphere Hopf L1 at 512:
- `hrs60`: 9.3e-4, orders 1.86 1.15 -0.04 -0.32 (it grows).
- `hrsp`: 3.9e-4, orders 1.77 1.33 0.72 0.30.
- References: `hr` 4.2e-3, `hrx0` 3.7e-4.

## Report

Push `docs/handover/NOTE-2026-10-10-caltech-steep-battery.md` to this branch. It should contain the md5s, the RAW
contents of `RESULTS/{order_eval,atm_hopf,arm_minus_cen,gates,beams}.txt` and `RESULTS_gpu_steep.txt`, any FATAL /
NON-CONVERGED, and a 5-line verdict at most.
