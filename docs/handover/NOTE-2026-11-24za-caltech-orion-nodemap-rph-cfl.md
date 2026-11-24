# NOTE 2026-11-24za Caltech -> orion: node map, RSG dumps for R_ph, cfl check (re 24x / 24y / 24z)

## Node map (11-24 ~01:50 PST)
| node | holder (ends) | run |
|---|---|---|
| hpc-sm-01-17, 02-03, 02-09 | 7333072 (11-26 05:21) | RSG rc3 (12 ranks), step .89 |
| hpc-sm-02-16 | 7333035 (11-26 23:16) | bsg rc3 production (FOFC + mg4), moved here from 02-09 |
| hpc-sm-02-11 | 7333073 (11-26 17:54) | m35mms rc3 production (FOFC + mg4) |
| hpc-sm-02-17 | 7333073 | t9l62: production PAUSED for the 24v layout arms (T0 running now), moved here from 01-17 |
| hpc-sm-01-15 | 7333036 (11-26 23:27) | M35ZAMS: production PAUSED for the 24v layout arms (Z3 running now) |
| (pending) | 7333037 | est. start 11-27 22:00 |

- Nothing was displaced. bsg and t9l62 moved to other holders' nodes at their FOFC/mg4 restarts, so that RSG gets 3 nodes in ONE holder (srun can't span holders).
- t9l62 and M35ZAMS production resume after their layout arms, on the winning layout if 24v's rule says switch.
- The layout results follow from the agent.

## RSG R_ph (24y)
The old deck's bin cadence is 1 d = ~86 cycles at dt 1000, so dumps already exist (rsg_rc3_1124/bin/):
- hydro_w + m1: rsg.*.08362.bin (restart, cycle 771891), 08363, 08364, 08365 (~cycle 772150).
- m1_vet: rsg.m1_vet.01673.bin.
Tell us which to copy where, or which reduction to run on them here.

## cfl (24z)
The latest rst of every production stores cfl_number = 0.3: bsg3d.00007, m35mms.00005, t9l62.00005, m35zams.00003, rsg.01708. None needs a restart.
The ZD0-ZD2 arms will name time/cfl_number=0.3 hydro/fofc_retest=true.

## R136
- r136.00149.rst is on Caltech: r136_rc4_1124/rst/, md5 0bca5931aa3788d86b58090be042f3e7, equal to viper's.
- Its 3 tables (ic + 2 opac) are copied too, md5s equal.
- It waits for rc4 + a free node.

-- Caltech
