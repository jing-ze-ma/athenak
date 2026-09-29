# NOTE 2026-09-29 (viper): WASP-121b max_eta scan SPLIT between DeltaAI and viper (user 09-29 ~23:20)

The user chose to split the max_eta scan of TASK-2026-09-29-deltaai-w121-mhd-etamax.md:
- **DeltaAI** runs the 12 short arms plus the noise twin on ghx4-interactive
  (NOTE-2026-09-29-deltaai-w121mhd-queue.md).
- **viper** runs the 4 long arms below.
- Both sites use the same arm definitions (deltaai_pkg/make_arm.sh and scan.sh) and NROT = 2 (tlim =
  rot 302), so the results are comparable and will be combined. viper reports when its arms finish.

## viper arms
Each arm is one apu job: 1 node x 2 MI300A, 2 ranks, no chain. The jobs were submitted 09-29 ~23:30 and
are pending (Resources); Slurm gives no start estimate yet (apu is full).

| arm | bbot (code) | max_eta | STS | dfloor | job | --time | athena -t |
|---|---|---|---|---|---|---|---|
| b3_e14_f13  | 3 G (0.8463)  | 1e14 | false | 1e-13 | 12028522 | 24:00:00 | 23:45:00 |
| b10_e14_f13 | 10 G (2.8209) | 1e14 | false | 1e-13 | 12028523 | 24:00:00 | 23:45:00 |
| b10_e13_f16 | 10 G (2.8209) | 1e13 | false | 1e-16 | 12028524 | 24:00:00 | 23:45:00 |
| b3_e13_f16  | 3 G (0.8463)  | 1e13 | false | 1e-16 | 12028525 | 06:00:00 | 05:47:00 |

- Arm directories: `/viper/ptmp2/jinma/w121_mhd_0929/viper_arms/<arm>`. The job script is
  `/viper/ptmp2/jinma/w121_mhd_0929/run_viper.sub` (ROCm 7.2 modules, HSA_XNACK=1,
  HSA_NO_SCRATCH_RECLAIM=1). The outputs (hst, rst every 0.5 rot, bins 4 per rotation) are as in the package.
- The binary is `/viper/ptmp2/jinma/builds/bin/athena_dhj_gpu72_ae767d20`, md5
  b91ad6e0dfd597df7b0ef39fc4544de7. Its source is rt-integration ae767d20, the same source as DeltaAI's
  350eee1c build. The old bbot key is bitwise with the new one.
- Expected wall time on viper, at ~21 ms/cycle on 2 MI300A (smoke 12027359):
  - e14 arms: Ohmic dt ~0.185 s gives ~25 h for 2 rotations. They may stop at rot ~301.9 on the 24 h limit.
  - b10_e13_f16: v_A-limited dt ~0.14 s gives ~33 h. Expect rot ~301.4 at 24 h.
  - b3_e13_f16: dt ~0.5 s gives ~2.6 h.
- An arm that does not reach rot 302 is reported as it stands and is not chained. Its rst files allow a
  continuation if the user wants one.
