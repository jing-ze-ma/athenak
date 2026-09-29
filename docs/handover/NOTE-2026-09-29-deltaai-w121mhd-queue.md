# NOTE for viper: check your queue for the w121 MHD max_eta scan (DeltaAI queue is jammed)

**From DeltaAI, 2026-09-29 ~16:15 CDT.** Refers to TASK-2026-09-29-deltaai-w121-mhd-etamax.md.

## Ask
Please check how soon viper would start **1-node jobs** for this scan and reply in a NOTE:
- the estimated start for a 1-node, 4 h job and a 1-node, 8 h job (`sbatch --test-only` or equivalent), in the
  partition you would use;
- the max wall time per job there;
- ms/cycle for MHD on 1 node, if you have it (the smoke gave 30.7 ms/cycle on 2 MI300A, including startup).

The user decides where the scan runs (all DeltaAI, all viper, or split) from your answer.

## Why
DeltaAI `ghx4` (the normal partition, 2-day limit): **1697 jobs pending; `sbatch --test-only` for a 1-node job
estimates a start on 10-05**. `ghx4-interactive` starts within ~1 h, but it caps jobs at **2 h** and charges
**2x** (8 SU per node-hour).

Rough scan cost on DeltaAI (from the TASK dt estimates; the measured ms/cycle comes from the smoke below):
about 25 node-hours for the 16 arms + noise twin. Most of it is 4 arms:
- b3_e14_f13 and b10_e14_f13 (Ohmic dt 0.185 s, no STS), ~5 h each;
- b10_e13_f16 (v_A-limited dt ~0.14 s), ~6 h, and b3_e13_f16 (0.46 s), ~2 h.
The other arms are ~0.1-0.5 h each. On `ghx4-interactive` the scan costs ~200 SU, and the long arms need chains
of 2 h restart links.

A natural split: DeltaAI runs the short arms (CFL- or cheap-Ohmic-limited) on ghx4-interactive; viper runs the
4 long arms above.

## DeltaAI state (nothing is running yet)
- Package unpacked at /work/nvme/bivj/jma20/w121_mhd_0929/deltaai_pkg (slim tgz md5 7f02a2c1 OK). The file
  `121_mhd_0929 scan (hydro production 1.0e-16)#` in deltaai_pkg/ came in the tarball (stray; harmless).
- 16 arm inputs made with scan.sh (NROT 2), plus the noise twin b3_e12_f13_r2 (to run at 2 ranks).
- Binary: rt-integration 350eee1c, dhj, CUDA HOPPER90 + ARMV9_GRACE, md5 722b75a461aef6d1bde13769e5f89c77.
  9e0b204a / 354a88b7 (bbot_gauss, units audit) are not in it; the arms use problem/bbot in code units,
  and the old key is bitwise.
- **User confirmed (09-29):** 3 G and 10 G are the polar field at the INNER wall (x1min = 1.127e10 cm):
  bbot 0.846 / 2.821 code. At x1max that is 1.0 / 3.3 G.
- run_arm.sub was adapted for DeltaAI as run_arm_deltaai.sub. It adds MPICH_GPU_SUPPORT_ENABLED=1, module
  reset, srun -c 16, --mem=0, a log-stall hang guard, and an athena -t limit of 10 min below the slurm limit.
- Smoke 3271618 (ghx4-interactive, pending): 300 fresh-start cycles each of b3_e13_f13, b10_e14_sts_f13 and
  b10_e13_f16 on 4 GH200 -> stability, dt/limiter, ms/cycle. Results will follow in a later NOTE.
