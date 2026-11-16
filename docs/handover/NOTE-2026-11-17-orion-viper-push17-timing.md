# NOTE 2026-11-17 Orion -> viper: push-17 CPU timing on Raven (the timing part of TASK-2026-11-17-viper-orion-push17-regression)

## Short answer
- **Push 17 e5c54306a recovers the push-16 CPU loss on cs_closed.**
  - cs_closed is -7 % vs push 16 and +3 % vs push 15.
  - Its "clean stage steps" wall is -7 % vs push 16 and +2 % vs push 15.
  - flz_q is +4 % vs push 15, inside the 5 % rule.
- ec96434ff (CleanLineRange) on its own, applied to dd468ab99 ("ddec"), gives:
  - flz_q +1 %, cs_closed +4 %, clean stage +4 % vs push 15;
  - cs_closed -6 % vs push 16.
- The remaining +3 % on flz_q in push 17 (2.613 vs ddec 2.536, about 3 spreads) comes from other push-17 content. A candidate is the new
  hydro and MHD flux-sweep barriers (a482258b5). We did not bisect it.
- **Bitwise:** on both decks, push 17 and ddec give byte-identical .bin and .hst files to push 16 dd468ab99.
- **Raven QOSGrpCpuLimit:** other workers' `small` jobs are pending on it. This `general` job started within about 2 min, so Raven was
  not capped for us.

## Setup (as in the bisect, NOTE 64e442cf3)
- **Machine:** Raven job 32444809, one dedicated node, ravc4053 (Xeon 8360Y, 72 cores, `general --exclusive`). It took 19 min, builds
  included. Compare Raven only with Raven.
- **Keys:** all keys off: clean_gas_exact, clean_fv_joint, clean_pc_factor, clean_yres = false and clean_rot_ring = off, written into
  a deck copy and on the command line. Other push-17 keys are at their defaults.
- **Decks:** flz_q (40 cycles, OMP 8) and cs_closed (4 cycles, OMP 6).
- **Order:** 2 repeats, the second round in reverse order.
- **Binaries:** built on Raven (gcc/13, openmpi/4.1, gpufix_reg recipe):

  | build | md5 |
  |---|---|
  | push 15 dfe21a746 | a80bdfce |
  | push 16 dd468ab99 | 666c9533 |
  | push 17 e5c54306a | 77b3563b |
  | ddec = dd468ab99 + ec96434ff | 85cd276b |

  To build ddec, we added only CleanLineRange to rad_m1_clean.hpp. The cherry-pick's other conflict was only context from push-17 code.

## Results (mean of 2 runs ± half-spread; % vs push 15, [%] vs push 16)

| build | flz_q s/cycle | cs_closed s/cycle | cs_closed clean stage steps (s) |
|---|---|---|---|
| push 15 dfe21a746 | 2.508 ±0.8% | 2.112 ±0.3% | 3.257 ±0.3% |
| push 16 dd468ab99 | 2.635 ±0.0% **+5 %** | 2.344 ±0.3% **+11 %** | 3.590 ±0.5% **+10 %** |
| push 17 e5c54306a | 2.613 ±1.0% **+4 %** [-1 %] | 2.177 ±0.5% **+3 %** [-7 %] | 3.321 ±0.4% **+2 %** [-7 %] |
| ddec | 2.536 ±0.6% **+1 %** [-4 %] | 2.202 ±0.0% **+4 %** [-6 %] | 3.385 ±0.1% **+4 %** [-6 %] |

- For comparison, P3 on 0678599eb (Raven job 32443460) gave +0 / +1 / +2 %. ddec is 1-2 % above that on cs_closed, which comes from
  the 0678599eb..dd468ab99 changes.

## Files (Raven, /raven/ptmp/jinma/stag2/p16bis)
- job17.sbatch, node17.log and p17tim.32444809.out;
- the tables come from `python3 ana.py Z dfe21a746,dd468ab99,e5c54306a,ddec`;
- the binaries are athena_dd468ab99, athena_e5c54306a and athena_ddec.

No restart dumps were left.
