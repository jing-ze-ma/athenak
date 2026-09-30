---
name: session-state-2026-08-18
description: "START HERE. End of 2026-08-18: origin/general-eos still d1c289dd, NO commits today, working tree clean under src/. All cost questions answered; the one open thread is the viper dt discrepancy. Six runs left in flight."
metadata:
  node_type: memory
  type: project
---

End of 2026-08-18. Supersedes [[session-state-2026-08-17]] as the entry point.

## NO CODE CHANGED TODAY

`origin/general-eos` is still `d1c289dd`. Nothing committed, nothing edited under
`src/ inputs/ tst/ docs/`. Today was measurement only. Everything new lives in
`/orion/u/jinma/ATHENAK/idflr/`.

## RESUME HERE -> [[dhj-viper-dt-discrepancy]] -- but it LARGELY DISSOLVED at session end

The alarm was "viper holds dt > 3 to t = 1e7 while the current binary decays to 1.34". The user
**re-checked and corrected it: viper's dt actually bottoms at ~2.6**, and `V_old` -- a build of
the viper commit `2af153a3` in a worktree -- reproduces that (dt = 2.599 at t = 3.03e5, min
2.491) and matches the old and new binaries to within 7%. So the general-EOS refactor does not
change this run.

**All that is left:** dt is NOT monotonic -- the current binary dips to **1.44 at t = 1e6** and
recovers to **1.52 at 1.5e6**. `V_old` had only reached t = 3.03e5 when the session ended.
Just read `idflr/V_old/stdout.txt` and compare dt over t = 5e5..2e6.

Binary built and kept: `idflr/athena_2af153a3`, worktree `/orion/u/jinma/ATHENAK/wt_2af153a3`
(remove with `git worktree remove` when done). **The branch is REBASED -- use `%cd`, not `%ad`:**
on 2026-08-11 the HEAD was `2af153a3`, and `aecba9f3`/`d2dae200` only landed 2026-08-12 06:00.

## Answered today (all measured, none acted on in code)

- **General EOS is 3.09x ideal per SIMULATED SECOND**, not 1.77x. The dt advantage that made it
  look cheap is transient and gone by t > 1e5. [[dhj-ideal-vs-general-cost]]
- **dhj's dt is the radial Alfven CFL at r/Rp ~ 1.28**, in a cell sitting exactly on dfloor.
  1.27% of the volume sets dt for the whole run. Raising dfloor x10 buys 2.68x.
  [[dhj-dt-limited-by-alfven-floor]]
- **`max_eta = 1e13` beats 1e14 for both EOSs**; at 3 G the margin is 2.9x, at 10 G 1.12x.
- **On the viper base at 3 G the sign REVERSES: 1e13 beats 1e12 by 1.38x**, because a lower cap
  suppresses diffusion, the field builds up, v_A rises and the Alfven CFL bites earlier. So
  max_eta acts on dt through FIELD STRENGTH here, not through dt_diff.
- **`pfloor=1`/`tfloor=10 K` is not worth it** on the tracked input -- the 10 K floor fires
  exactly zero times and pfloor merely heats the thinnest gas 10x for no dt gain.
- **perna vs the EOS x_e table cost the same** (within 1-3%) once max_eta is low enough that
  diffusion never binds.

## Still not done (unchanged from yesterday)

- `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` still has `max_eta = 1e14` + `use_rkg_sts =
  true`. Today's evidence says 1e13 + STS off. **NOT CHANGED -- the user has not decided.**
- General EOS Stage 4 (rho,e) table: still not started ([[general-eos-stage4-rho-e-table]]).
- solar_convection atmospheric runaway: untouched ([[solar-convection-general-eos]]).

## LEFT RUNNING -- do not cancel

Six jobs, all 3 G, all from the viper base, 23:30 wall from ~12:46 and ~15:30:

| job | dir | what |
|---|---|---|
| 190502 | `V_perna_12` | perna, max_eta 1e12 |
| 190503 | `V_perna_13` | perna, 1e13 |
| 190504 | `V_eos_12` | eos x_e, 1e12 |
| 190527 | `V_exact` | **the original file byte-for-byte**, 2 blocks, 2 ranks x 56 thr |
| 190528 | `V_bc_on` | 1e12 + hst/dumps, Maxwell term ON |
| 190546 | `V_old` | **the viper commit `2af153a3`**, original file, 2 blocks, 2 ranks x 56 thr |

`V_eos_13` (190505) and `V_bc_off` (190529) were cancelled to free slots (the 6-job assoc cap);
`V_bc_off` had already answered its question. Resubmit either if wanted.
The four 10 G runs and the four earlier 3 G runs were all cancelled today.

## Files in `idflr/` (all new today)

`base_viper.athinput` (the original + MeshBlock 64x32x16 -- a TWO-LINE diff),
`base_viper_eos.athinput` (+ `<units>` + `eos_*` keys + `ohmic_resistivity=eos`),
`base_viper_bc.athinput` (+ an explicit `bc_outer_maxwell` line),
`base_viper_2blk.athinput` (byte-identical to the original),
`run_long_2rank.sh` (2 ranks x 56 threads).

## Three operational traps learned today

1. **AthenaK command-line overrides can only MODIFY a parameter that already exists in the
   input file.** `problem/bc_outer_maxwell=false` killed two jobs instantly because the key was
   absent. Add the line to the file first.
2. **The assoc limit is 6 concurrent jobs** (`sacctmgr show assoc user=jinma` -> MaxJobs=6);
   a 7th sits in `AssocMaxJobsLimit`.
3. **`athena -c` says `OpenMP parallelism: OFF` and that is NOT the Kokkos backend.** Kokkos
   OpenMP is on; verify occupancy with `scontrol show node <n> | grep CPULoad` -- it read
   112.00 on every node. 16 ranks x 7 threads = the full 112 physical cores.

## MeshBlock sizing, measured

64x8x8 -> 2.39x ghost overhead; **64x32x16 -> 1.49x**. Predicted 1.60x, **delivered 1.24x per
real cell**: the ghost work fell 29% but the per-cell-update rate got 29% WORSE (32,768-cell
working set vs 4,096). 2 blocks on 2 ranks x 56 threads is slower still (0.0191 vs 0.0177
s/cycle) -- OpenMP scales less well than ranks here.
