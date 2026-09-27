---
name: dhj-viper-dt-discrepancy
description: "LARGELY DISSOLVED 2026-08-18 evening. The user re-checked: viper's dt bottoms at ~2.6, NOT 'never below 3', and V_old (a build of the viper commit 2af153a3) reproduces it. Old and new binaries track within 7% to t=3e5. ONE question left: the newer binary dips to dt~1.44 at t=1e6 -- does V_old do the same?"
metadata:
  node_type: memory
  type: project
---

Opened and mostly closed 2026-08-18.

## RESOLUTION SO FAR: there is little or no discrepancy up to t = 3e5

The original report was "viper holds dt > 3 out to t = 1e7". **The user re-checked and
corrected it: viper's dt actually goes down to ~2.6.** With that number the alarm mostly
disappears -- `V_old`, built from the viper commit, sits at **dt = 2.599 at t = 3.03e5 with a
minimum of 2.491**, which is exactly "2.6 ish". The user confirms V_old matches their viper run.

**Old binary vs new, same input, same decomposition, only the binary differs:**

| t | V_old (2af153a3) | V_exact (d1c289dd) | ratio |
|---|---|---|---|
| 1e5 | 3.920 | 3.917 | 1.001 |
| 1.5e5 | 3.101 | 2.957 | 1.048 |
| 2e5 | 2.796 | 2.820 | 0.991 |
| 2.5e5 | 2.682 | 2.851 | 0.941 |
| 3e5 | 2.749 | 2.781 | 0.989 |

Within 7%, sign wandering = chaotic decorrelation, not a systematic code effect. **So the
entire general-EOS refactor (aecba9f3, d2dae200 and the HLLD/HLLE rewrite) does NOT change this
run's dt** over the range compared so far.

## THE ONE QUESTION LEFT

dt is **NOT monotonic**. `V_perna_12` (current binary) reads 2.561 at t=5e5, **1.585 at 8e5,
1.436 at 1e6**, then recovers to **1.522 at 1.5e6**. So there is a DIP around t ~ 1e6, not a
collapse. Whether the viper commit shows the same dip is unknown -- `V_old` had only reached
t = 3.03e5 at session end.

**NEXT SESSION: read `idflr/V_old/stdout.txt` and compare its dt over t = 5e5 .. 2e6 against
`V_perna_12` / `V_exact`.** If V_old dips to ~1.4 too, the matter is fully closed and the whole
thing was a misremembered number plus chaotic scatter. If V_old stays near 2.6 through the dip,
the code-version question reopens -- but only for t > 5e5, and the bisect list is
`aecba9f3, d2dae200, eb39f8db, 1c3cf1e7, de4f5180`.

A waiter was armed at t = 6e5 (task b1jgcrtvk) but the session ended; just read the file.

## What was ELIMINATED along the way (all by experiment)


1. **MeshBlock decomposition.** `V_exact` runs the original byte-for-byte (2 MeshBlocks, 2
   ranks x 56 threads) against `V_perna_12` (16 blocks, 16 ranks x 7). dt agrees to 4-5
   significant figures: 10.145/10.145 at t=500, 8.131/8.129 at 1e4, 5.806/5.799 at 2e4.
   Roundoff only. **My 64x32x16 change is not the cause.**
2. **The outer-x1 Maxwell term.** `V_bc_off` (`bc_outer_maxwell=false`) vs `V_bc_on`: ratio
   1.004 / 1.078 / 1.073 / 1.087 / 0.943 at t = 1e5..1e6. **Under 10%, and it changes sign.**
   Since removing the term ENTIRELY is worth <10%, the `bc6b5774` -> `767fff43` reformulation
   of that term (raw ghost-energy subtraction -> clamped effective gravity) is bounded by the
   same 10% and cannot explain 2.2x. The pgen's own note agrees: the term is 3% of the
   hydrostatic term at 3 G.
3. `c0edbe49` (eta_b out-of-bounds) -- `constant` type and 1D/2D only, the commit says 3D is
   in bounds. `d52f3610` -- GPU host pointers, inert on CPU. `2a154e7d`/`d1c289dd` -- need the
   well-balanced scheme (commented out in this file) or the general EOS (off).

## What the physics is doing

First hst ever obtained for this configuration (`V_bc_on`, output1/dt=1e4): **the field
amplifies by ~4 orders of magnitude.** 1-ME goes 4.77e10 (t=0) -> 5.74e14 (2e5) -> 7.40e14
(4e5) -> 6.42e14 (8e5); 3-ME 0 -> 1.03e14. Mass conserved to 0.14%, tot-E to 0.5%. dt decays
as the field grows -- consistent with [[dhj-dt-limited-by-alfven-floor]] (radial Alfven CFL in
the floored outer atmosphere). **On viper the field must amplify far less, or the top stays
denser.**

## THE VIPER BINARY IS `2af153a3` -- CONFIRMED BY THE USER + GIT

The user built it on **2026-08-11**, and the branch was **REBASED**: `aecba9f3` and `d2dae200`
carry author date 08-07 but **commit date 2026-08-12 06:00**. So on 2026-08-11 the branch HEAD
was `2af153a3` "small bug fix" (commit date 2026-08-07 21:35) -- **before the entire general-EOS
series**. ALWAYS check %cd, not %ad, on this branch.

`2af153a3 .. d1c289dd` over the files an ideal+perna run touches is 5,761 insertions / 309
deletions across 35 files, including **`rsolvers/hlld_mhd.hpp` (+32)**, `hlle_mhd.hpp` (+84),
`mhd_fluxes.cpp` (+190), `mhd_newdt.cpp` (+23) and the pgen (+458). An ideal-gas MHD run calls
HLLD on every face, so the Stage 1/2 general-EOS refactor is the prime suspect.

## IN PROGRESS: building 2af153a3

`git worktree add --detach /orion/u/jinma/ATHENAK/wt_2af153a3 2af153a3`, then the SAME flags the
current binary used (recovered from `build_dhjrt/CMakeCache.txt`):
`-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release
-DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_SERIAL=ON -DAthena_ENABLE_MPI=ON -DKokkos_ARCH_SPR=ON
-DPROBLEM=deep_hot_jupiter_rt`, modules `gcc/13 openmpi/4.1 cmake`. Binary lands at
`idflr/athena_2af153a3`. The pgen DOES exist at that commit.

Then: run `idflr/base_viper_2blk.athinput` (byte-identical to the original) with it to t ~ 3e5
and compare dt against `V_exact`. If dt >= 3, viper is reproduced -> bisect forward over
`aecba9f3, d2dae200, eb39f8db, 1c3cf1e7, de4f5180`.

## OLD NEXT STEP (superseded, kept for the candidate list)

The input is now proven irrelevant. Ask the user **which commit the viper binary was built
from** -- the Aug-11 file date was only ever an assumption. If unknown, bisect: build in a
`git worktree` (do NOT check out in the main tree) with `-DPROBLEM=deep_hot_jupiter_rt` and run
the original file to t = 3e5, comparing dt. Candidate commits, newest first, all touching MHD
or this pgen: `767fff43`, `bc6b5774`, `5d1c3436`, `c0edbe49`, `1c3cf1e7`, `eb39f8db`,
`de4f5180`, **`d2dae200` (refactored HLLD for the general EOS -- an ideal-gas run uses HLLD)**,
`aecba9f3` (Stage 1 hydro/MHD). d2dae200/aecba9f3 are the strongest remaining suspects: they
touched the solver an ideal-gas run actually calls.

Also still unchecked: whether the viper build was single precision, or a GPU/Kokkos backend.

Related: [[dhj-dt-limited-by-alfven-floor]], [[dhj-ideal-vs-general-cost]], [[general-eos-project]].
