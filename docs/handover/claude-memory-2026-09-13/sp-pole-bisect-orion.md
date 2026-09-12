---
name: sp-pole-bisect-orion
description: The sp polar-row blow-up bisection restarted on orion CPU (jobs 193408-11); how the production grid was RECONSTRUCTED and verified, the one number that still does not match, and why the controls decide
metadata:
  type: project
---

Round 2 of [[sp-pole-bottom-radial-blowup]], rerun on orion CPU because viper is down.
`/orion/ptmp/jinma/Athenak/sp_pole_bisect`. Four detached worktrees + CPU MPI builds
(gcc/13, openmpi/4.1, `Kokkos_ARCH_SPR`), kokkos symlinked to the main repo (the submodule
pointer is d8e9af03 at ALL four commits, so one checkout serves them all):

| arm | commit | role |
| --- | --- | --- |
| 18f5dd21 | the dirty binary | CONTROL, must blow up |
| 2a64e2c7 | the clean binary | CONTROL, must stay flat |
| 4990eb41 | 18f5dd21 minus fefe6a17 (the wall mirror) | candidate |
| 863e8337 | the x3-face shift | candidate |

## The input had to be RECONSTRUCTED, and my first attempt was WRONG

Viper's `bench/` input is not in git. My first attempt took the repo input as of 18f5dd21
verbatim -- nx1 = 64 UNIFORM, `f_stretch_theta = 2.0`. The user challenged it, correctly.
The notes actually pin the production grid: [[sp-mhd-energy-excess]] says "nx1 128 refit
stretch, 2.8 deg" and "the theta stretch (f = 3, only on sp)".

**The check that settles the theta grid, and it is exact.** [[sp-pole-bottom-radial-blowup]]
records one geometric number, "8.1 deg" polar cells. `StretchTheta` in
`src/coordinates/grid_stretch.hpp` is
`theta = pi/2 (1 + sinh(a(2 xi - 1))/sinh(a))`, so the polar cell width is a pure function
of (nx2, a). nx2 = 64: a = 2 gives 5.66 deg, **a = 3 gives 8.10 deg** -- and its uniform
equivalent is 2.81 deg, the note's "2.8 deg". So nx2 = 64, f_stretch_theta = 3.

Final reconstruction: nx1 128 + refit poly stretch (-0.068392, -2.191487, 2.464818,
-1.366698), x1max 2.0556e10, nx2 64 f_stretch_theta 3, nx3 128, meshblock 128 x 16 x 16
(32 blocks), point-mass gravity, bbot 3 G, pfloor 1e-3, dfloor 5e-13, grey RT
(rt_ck = false, the era default), IDEAL MHD. Resistivity is off the way the code wants it:
the module is constructed iff `<mhd>/ohmic_resistivity` EXISTS, so the key is COMMENTED OUT,
not set to a null value. Round 1 showed the ideal arm still blows up, so this is the right
cheap arm.

**The meshblock came from the gate script, not from a note.** `polerow.py` probes j = 8
INSIDE a block, which is out of bounds for a 8-cell theta block -- so viper's theta block
was >= 9. Combined with [[meshblock-decomposition-gpu]] (32 blocks is the measured GPU
optimum) that gives 128 x 16 x 16 = 32 blocks. Block layout is not cosmetic for a POLAR
halo bug, so this was worth restarting the arms for.

## STILL NOT MATCHED: dt is 19.2 s here, 2-3 s there

The 200-cycle timing run gives dt 19.2 s and 1.37e7 zone-cycles/s on one node (16x7);
the notes record "~1 rot/h (dt 2-3 s)" for the arms being bisected. A 7x dt gap means
something in the config is still not viper's -- most likely the RT mode or the resistivity
settings. **Therefore the two CONTROLS decide whether any of this is meaningful**: if
18f5dd21 does not separate from 2a64e2c7 on this config, the reconstruction is wrong and
the candidate arms must be DISCARDED, not believed. Do not read the candidates first.

Gate: history column 11 (`1-ME`, radial magnetic energy) at rot 4.7, and `polerow.py` on
the rot-4 dump. Initial ME1 = 4.77e30 here (bbot 3 G; viper's baseline 8e31 was at its own
field and grid, so compare arms against EACH OTHER, never against viper's absolute number).

Cost: ~1.7 h per arm at dt 19 s, 8 h slots, rst every rotation, `arm.sh` resumes from the
newest rst. Jobs 193408 (18f5dd21), 193409 (2a64e2c7), 193410 (4990eb41), 193411 (863e8337).

## Status at 2026-09-08 01:55 (orion)

Jobs 193408-11 running, second slots 193425-28 queued behind them. Gate is ME1 at
rot 4.7; nowhere near it yet, and the three arms that must SEPARATE have not:

    18f5dd21 (dirty control)  rot 2.15  ME1 5.204e+31
    2a64e2c7 (clean control)  rot 1.86  ME1 4.669e+31
    4990eb41 (candidate)      rot 2.14  ME1 5.206e+31
    863e8337 (candidate)      rot 2.14  ME1 5.206e+31

Both candidates are tracking the DIRTY control to four digits, and the clean control is
behind in rotation so its lower ME1 is not yet evidence of anything -- compare at EQUAL
rot, not at equal wall time. Nothing is decidable before rot 4.7. Just let them run and
re-read this table.
