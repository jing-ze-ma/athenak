# NOTE 2026-11-24v orion -> Caltech: layout batch for t9l62 + M35ZAMS (USER: stop runs for quick tests OK; bsg / m35mms keep)


USER (11-24 ~04:00 PST): "try find the optimal block layout for those 4 stars. they don't need to be the same." USER via COORD: "you can stop those 4 runs if you need to test. but make it quick."

## Summary
- **bsg and m35mms: KEEP. No arms, no stop.** Each already runs 1 block / GPU (4 blocks, 256x128x128 / 336x128x128) with tile 128, which is the largest tile that fits 4 GPUs. Every other layout (2-4 blocks / GPU, or 1 / GPU strips 256x64) needs tile 64. That gives a new answer and a weaker tile-local preconditioner, so it has no expected gain.
- **t9l62** (8 blocks 512x128x128, 2 / GPU) and **M35ZAMS** (16 blocks 384x32x32, 4 / GPU) are the two candidates. Each gets a fresh 200-cycle layout A/B.
- A restart can't change the layout (FATAL). A star that switches therefore starts fresh. It switches only if the measured gain beats the decision rule in the last section.

## Common setup for all arms
- Binary: `/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3`, md5 09f50abe. 4 ranks = 4 H200, 1 node. Use the production launcher settings (prodg.sh: romio321, pmix).
- Decks: prod_rc3_1124/decks/tile/IA123_t9l62.athinput and IA123_m35zams.athinput (= ia123rc3_decks/IA123_m35zams_wallsym.athinput, md5 f29145b5).
- Name these keys on EVERY arm (current production keys; rst-stored-key rule): `hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4`
- Run length and outputs: `time/nlim=250` for fresh arms. Outputs off except hst: `output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`, plus any further bin/rst blocks of the deck. Use a scratch run dir under prod_rc3_1124/layopt/<arm>/. Delete it after posting, and keep no rst.
- Return per arm:
  - s/cyc over cycles 50-250, from the 50-cycle elapsed prints;
  - setup time (launch to cycle 0);
  - peak GPU memory per GPU: run `nvidia-smi --query-gpu=index,memory.used --format=csv -l 10 > mem.csv` in the background during the run, and report the max per GPU;
  - mean BiCGStab its / stage and solves / cycle from the CLEAN lines;
  - CW / FATAL / fallback counts;
  - md5 of the first 201 non-# rows of user.hst and hydro.hst, and dt at cycle 200.

## Arms
### t9l62 (mesh 512 x 256 x 512; x1 stays whole, as the deck requires ("vet_gd: one x1 block"))
| arm | blocks | meshblock (nx1 x nx2 x nx3) | per GPU | layout_tile | answer vs T0 | command-line delta |
|---|---|---|---|---|---|---|
| T0 | 8 (2 theta x 4 phi) | 512x128x128 | 2 | 128 | baseline (current production) | none |
| T1 | 4 (1 theta x 4 phi) | 512x256x128 | 1 | 128 | **identical** (tile rule; x2-only change, LAYMOM-safe; the GPU PCR t[2] fix 3efc8cb24 is in rc3) | `meshblock/nx2=256` |
| T2 | 4 (2 theta x 2 phi) | 512x128x256 | 1 | 128 | **identical** (tile rule; x3-only change, 2 tiles / block along phi -> LAYMOM fwb_ring per-tile fix is in rc3) | `meshblock/nx3=256` |

### M35ZAMS (mesh 384 x 128 x 128; x1 whole)
| arm | blocks | meshblock | per GPU | layout_tile | answer vs Z0 | command-line delta |
|---|---|---|---|---|---|---|
| Z0 | 16 (4x4) | 384x32x32 | 4 | 32 | baseline (current production) | none |
| Z1 | 4 (2x2) | 384x64x64 | 1 | 32 | **identical** (tile rule) | `meshblock/nx2=64 meshblock/nx3=64` |
| Z2 | 4 (1 theta x 4 phi) | 384x128x32 | 1 | 32 | **identical** (tile rule; theta whole; theta-periodic wedge, so the mirror-wall "1 block along x2" constraint is satisfied anyway) | `meshblock/nx2=128` |
| Z3 | 4 (2x2) | 384x64x64 | 1 | **64** | **NEW run** (round-off different: the tile-local MG preconditioner gets 64-cell tiles) | `meshblock/nx2=64 meshblock/nx3=64 mesh/layout_tile=64` |
| Z4 (optional, only if a slot is free) | 8 | 384x64x32 | 2 | 32 | identical (x3-only) | `meshblock/nx2=64` |

### ZAMSFAST's M35ZAMS dt arms (ZAMSFAST owns them; RESTARTS from the latest m35zams production rst, same 16x(384x32x32) layout, 300 cycles, no bin/rst outputs)
- Name explicitly on each: `rad_vet/clean_stag_mg_levels=4 hydro/fofc=true rad_vet/floorcap_v=3e7`.
- **ZD0** DTREP base: add `hydro/dt_report=10`. Return:
  - the DTREP / DTREG lines;
  - s/cyc over cycles 50-300 (this is also the late-state s/cyc for the layout projection);
  - CLEAN solves and BiCGStab its per cycle;
  - the floorcap_report lines (25).
- **ZD1** (diagnostic, not for adoption): ZD0 + `rad_vet/floorcap_v=1e7`. Return dt, the DTREP binder, the floorcap clipped-KE booking, and s/cyc.
- **ZD2** (only if a slot is free): ZD0 with `rad_vet/clean_stag_mg_levels=2` named.

## Node plan (COORD 04:xx: batch goes FIRST, then RSG takes 02-03 / 02-16 / 02-17; all arms start at once)
1. Stop t9l62 (hpc-sm-01-17) and m35zams (hpc-sm-01-15) gracefully: USR1 = rst + exit, as for bsg at 00:42. ZD0 / ZD1 / ZD2 restart from the m35zams rst written here.
2. bsg (02-09) and m35mms (02-11) **keep running**: they need no arm. Stopping them would only cost progress.
3. Start all at once:
   | node | arm(s) | est. wall | note |
   |---|---|---|---|
   | hpc-sm-02-17 | T0 | ~12 min | timing |
   | hpc-sm-02-03 | T1 | ~12 min | timing |
   | hpc-sm-01-17 | T2 | ~12 min | timing (prod t9l62 node, freed in step 1) |
   | hpc-sm-01-15 | Z0, Z1, Z2, Z3, Z0 repeat, in series (Z4 if time) | ~5 x 4 = 20 min | all Z arms on ONE node, so the ratios are internal; the Z0 repeat brackets drift |
   | hpc-sm-02-16 | ZD0, ZD1 (ZD2 if time) | ~3 x 5 min | timing-suspect node: ZAMSFAST's arms are about dt / binder, not s/cyc. ZD0's s/cyc is cross-checked against Z0 on 01-15 |
4. The test nodes are free for RSG after ~15 min (02-17 / 02-03 / 02-16). 01-15 is done after ~20-25 min.
5. Then each star follows "Stop / restart steps" below:
   - **keep** = resume from the step-1 rst with the production command line + `rad_vet/clean_stag_mg_levels=4` (planned for the next rst anyway), so nothing is lost beyond the ~25 min pause;
   - **switch** = fresh start with layout X.

## Predictions (from what we measured; honest ranges)
- **What does NOT change with blocks / GPU on one node:**
  - The PCR line solve: x1 lines, 32768 lines / GPU for t9l62 in every layout; the team size is set by nx1, so the 25 % occupancy is the same.
  - The Krylov dot / update kernels: AthenaK packs all of a rank's blocks into one launch.
  - Halo bandwidth: intra-node NVLink, where MPILAT measured 0 % lever on single-node H200.
  - In tile mode, the preconditioner strength: it is set by the tile, not the block.
- **What does change:**
  - per-block host loops and buffer packing, mainly the vet_gd per-shell ray halo (251 exchanges / cyc on m35mms; ~1 per shell on t9l62's 512 shells, per block);
  - intra-rank face copies;
  - per-block launches outside the packed kernels.
  - This is small on big blocks and large on tiny ones.
- **Per-arm predictions:**

  | arm | prediction | basis |
  |---|---|---|
  | T1 / T2 vs T0 | -5 to -15 % s/cyc (central ~-8 %) | HALFDAY-T9's 1.35x is an upper bound (it compares A100 2-node vs H200, not the same hardware). T1 has fewer off-rank faces than T2 (2 x 256 vs 2 x 256 + 2 x 128 lines / block), within noise. |
  | Z1 / Z2 vs Z0 | -10 to -25 % | 1.57 M cells / GPU (~1/10 of t9l62's), so latency / per-block bound: 4 blocks of 32^2 x 384 per GPU. |
  | Z3 vs Z1 | a further -3 to -10 % from fewer BiCGStab its (stronger 64-tile PC); new answer | |

- **Memory:**
  - t9l62 1 / GPU holds the same 16.8 M cells / GPU as today. Expect today's value or slightly less (fewer ghost cells), well under 141 GB. Scaled from m35mms 36 GB for 5.5 M cells: ~100-110 GB. Please confirm with mem.csv.
  - M35ZAMS: ~10-12 GB / GPU in every layout.

## Decision rule and expected time to tlim (fill in with the arm numbers)
Let s0 = baseline s/cyc, sX = best arm, N = total cycles, Nd = cycles done by the current production at decision time. A fresh start with X beats keeping if sX / s0 < (N - Nd) / N.
| star | started | N (est.) | Nd at the step-1 stop (est.) | break-even sX/s0 | keep: time to tlim | switch: time to tlim at the predicted gain |
|---|---|---|---|---|---|---|
| bsg | 00:25 | ~29k (dt ~170) | — | — | keep (~4-5 h left) | — |
| m35mms | 00:25 | ~72k (dt 12.0) | — | — | keep (~16-18 h left) | — |
| t9l62 | 00:25 | ~51k (HALFDAY-T9) | ~6k (2.35 s/cyc) | **< 0.88** | ~45k x 2.4-3.0 s = 30-38 h | -8 %: 51k x 2.2-2.8 = 31-39 h (a wash); -15 %: 29-36 h (saves ~1-2 h) |
| M35ZAMS | 00:25 | ~3.5e5 at dt 2.49 (2.5e5 at the 3.5 s ceiling, ZAMSFAST) | ~4e4 (0.37 s/cyc) | **< 0.89** | ~3.1e5 x 0.37 = ~32 h | -15 %: 3.5e5 x 0.31 = ~30 h; -25 %: ~27 h (saves ~5 h) |

- Use ZD0's late-state s/cyc for M35ZAMS's s0, not the fresh value, and multiply the fresh-arm ratio sX / s0 by it (ZAMSFAST: dt 2.49 s late, ~3.5 s ceiling).
- If the gain is under the break-even, keep.
- A switch to T1 / T2 / Z1 / Z2 keeps the CURRENT answer bitwise (check: the hst md5 equals T0 / Z0). The fresh restart then reproduces the current run's history exactly and simply runs faster. A switch to Z3 is a new run.
- If T1 / T2 / Z1 / Z2's hst md5 does NOT equal the baseline, post the max rel diff of L_top and do not switch: that would be a layout bug to report, not a layout choice.

## Stop / restart steps per star after the arms
- bsg, m35mms: nothing. At their next rst, add rad_vet/clean_stag_mg_levels=4 (NOTE 24s).
- t9l62 / M35ZAMS, if KEEP: restart from the step-1 rst with the production command line + `rad_vet/clean_stag_mg_levels=4`.
- t9l62 / M35ZAMS, if SWITCH to X:
  1. Keep the step-1 rst until the new run passes its watch.
  2. Fresh-start in a new run dir prod_rc3_1124/<star>_L<X>/ with `-i <deck>`, the production command line (FOFC keys), `rad_vet/clean_stag_mg_levels=4` and the arm's meshblock / layout_tile delta. Use the same holder and node.
  3. Same 300-cycle watch as for FOFC; also check that the hst md5 of rows 0-200 equals that arm's. Post s/cyc.

-- orion (coordinator)
