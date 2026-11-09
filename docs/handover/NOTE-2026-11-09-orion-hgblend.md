# NOTE 2026-11-09 Orion -> viper: hgblend CPU gates for `<rad_m1>/implicit_blend_t_upwind` (answer to TASK-2026-11-09-orion-hgblend)

Short answer: all criteria PASS.
- The key is active on the small fresh HG state: O4 differs from N4 in 99 % of the cells.
  - Final state: max |dE|/E = 3.3e-6 and max |dT|/T = 9.2e-7 over all cells; 3.2e-7 and 8.1e-8 at tau > 1.
  - Energy: total-energy change at round-off (-2.7e-16).
- Keys absent: N4 and R4 are bitwise identical (whole files).
- Rank layout: O1 and O4 have a bitwise identical state. The hst differs only at the MPI_Reduce round-off level, 2.4e-13.
- Restart: bitwise.
- No stage fault occurred (min E from the solve 0.536), so no XSPLIT_BLEND line was printed.
- Cost on CPU: the first two steps take ~3.3x the Picard passes; wall time is +51 % over 6 cycles.

## Commits, builds, jobs
- NEW = `1f314023a` (origin/hgpos-1109b): detached worktree `/orion/ptmp/jinma/hgblend_1109/wt_new`, kokkos d8e9af03.
  Built with the TASK recipe (gcc/13 openmpi/4.1 cmake/3.30, Release, MPI, he_star_m1, `-ffp-contract=off`), `make -j32`:
  120 s, 0 warnings.
- REF = `b798f1a30`: the hgpos NEW binary, reused (`/orion/ptmp/jinma/hgpos_1109/wt_new/build/src/athena`, same recipe).
- Deck: `docs/handover/hgpos-1109/hg_cs_small.athinput` + `files/` from bsg-files-1009. All 4 updated MD5SUMS OK.
  It now carries `layout_invariant = true`.
  - `implicit_blend_t_upwind` is not in the deck, and a fresh start cannot add an absent key from the command line.
    So `on.athinput` = deck + `implicit_blend_t_upwind = true` right after `<rad_m1>` (sed).
    O4, O1, OA and D4 use it; R4 and N4 use the deck as is.
  - `rad_m1/implicit_blend_t_upwind=true` was also passed on the command line in every ON run, OB included.
- Keys for all runs: `time/nlim=6 output1/dt=1e-30 output2/dt=1e30 output3/dt=1e30 output6/dt=1e30` (OA: nlim=3).
  - OB = `-r OA/rst/hegiant.00001.rst` (cycle 3) + the same keys + nlim=6 + ON.
  - D4 = ON + `rad_m1/dbg_t2_admiss=50`. The deck already has `dbg_t2_admiss = 20` in every arm.
- Slurm (all on p.shared, `--qos=serial`):
  - **227953** = O1 alone, from the start (1 rank, `--mem=120G`, orion001). MaxRSS 25.2 GB.
  - **227954** = R4 N4 O4 D4 OA->OB (`-n 20 --mem=150G`, each step `srun -n 4 --exact --mem=24G --cpu-bind=cores`, orion002).
    It pended ~5 min on AssocMaxJobsLimit.
    - R4 and N4 finished (MaxRSS 5.78 / 5.80 GB per task).
    - **O4, D4 and OA were OOM-killed at the 24G step limit.** O4 and D4 died after cycle 6, OA after cycle 3, each while
      writing the final restart. sacct MaxRSS at the kill was 6.29 GB per task, so 4 x 6.3 GB > 24G.
    - OB had then started from OA's t = 0 restart; I cancelled that step (227954.5).
    - The killed runs' logs are kept in `runs_oom24G/`.
  - **227959** = O4 D4 OA->OB rerun with `--mem=48G` per step (orion002). MaxRSS 5.82 GB per task.
    The numbers below are from this job and 227954 (R4, N4) / 227953 (O1).
  - Note: the ON key needs a little more than 6 GB per rank at peak on this deck, against ~5.8 GB without it.
    24G per 4 ranks was enough for N4 and R4 (and for the hgpos arms), but not with the key.
- Run dirs: `/orion/ptmp/jinma/hgblend_1109/runs/<id>` (cmd.txt, log, hst, wall.txt).
  All rst and bin were deleted after the comparisons.
- Scripts, in `/orion/ptmp/jinma/hgblend_1109/scripts/`:
  - `bincmp.py` and `rstcmp.py` compare the data payload after the header (the header embeds the deck).
  - `dev_cs.py`: max |dE|/E and |dT|/T per MeshBlock from the final restarts.
    - E = m1 u0 and T = wtemp x 1.202724e-8.
    - tau = radial sum from the top of rho kappa_R(rho, T) dr, on N4's state, with the deck's Rosseland table. dr comes
      from the version-1.2 x1f coordinate arrays of the bin dump.
    - It uses `rr_cs.py`, the consfix reader with the cubed-sphere restart block `M1ROTAV1` added.
    - Check: T/T_rad = 1.0000 at tau > 3. 445,440 of 1,081,344 cells have tau > 1 (tau_max 2.7e6).
- Every run: 6 cycles (Terminating on cycle limit, t = 857.14 s, dt = 142.86 s): cycle 1 backward Euler, cycles 2-6 hesdirk2.

## Criteria
1. **Keys absent, N4 vs R4: PASS.** `cmp`:
   - hegiant.hydro.hst and hegiant.user.hst identical;
   - all bin files (hydro_w, m1, m1_vet at 00000 and 00001) identical as whole files;
   - rst 00000 and 00001 identical as whole files.
2. **Rank layout, O1 vs O4: PASS.**
   - The data payloads of hydro_w, m1 and m1_vet are identical at 00000 and 00001.
   - The rst payload is identical at 00000 (538,235,662 B) and at the final 00001 (2,595,577,318 B, all double-precision state).
   - hst: time and dt are identical in all 7 rows. The rest differs at round-off: user.hst max rel 2.44e-13 (col 9, e_gas),
     hydro.hst 1.92e-13. This is the same MPI_Reduce effect as in hgpos (outputs/history.cpp:415): it is already there at
     t = 0, where the state is bitwise identical.
   - D4 vs O4: hst and the final rst payload are bitwise identical, so the dbg key does not change the result.
3. **Restart, OA+OB vs O4: PASS.**
   - The last 3 rows of both hst files are identical.
   - The final hydro_w, m1 and m1_vet payloads are identical.
   - OB's final rst (00002) payload is identical to O4's final rst (00001).
   - OB: stage steps 3, backward-Euler steps 0.
4. **Active path: PASS (O4 differs from N4).**
   - About half of the float32 words of the final bin dumps differ. 1,071,904 of 1,081,344 cells have a different E.

   | final state, O4 vs N4 | max \|dE\|/E | max \|dT\|/T |
   |---|---|---|
   | all cells | 3.326e-06 | 9.150e-07 |
   | tau > 1 | 3.240e-07 | 8.091e-08 |
   | tau < 1 | 3.326e-06 | 9.150e-07 |

   - The largest dE/E sits at the outer boundary cell: MeshBlock 12, i = 703, rho = 1e-14 (the floor), T = 3072 K,
     tau = 3e-7. The deep interior moves by < 3.3e-7.
   - Counters (O4; D4 identical):
     ```
     <rad_m1> implicit_blend_t_upwind (rank 0): x2/x3 face evaluations=9.400000e+01 faces with |q| capped at c=0.000000e+00
     <rad_m1> time_scheme=hesdirk2: stage steps=5.000000e+00 backward-Euler steps=1.000000e+00 stage fallbacks=0.000000e+00 ...
     <rad_m1> implicit transport: solves=1.100000e+01 Picard iterations mean=8.545455e+00 max=2.400000e+01 NON-CONVERGED=0
     <rad_m1> bicgstab: outer passes=9.400000e+01 linear solves=9.400000e+01 inner iterations mean=2.595745e+00 max=4 total=244
     min E from the solve=5.364873e-01 ; floor clips: solved E <= e_floor cell-solves=0, written-back gas eint <= 0 cell-solves=0
     ```
   - N4: stage steps 5, backward-Euler steps 1, Picard mean 4.18, max 8; bicgstab outer passes 46, inner total 114;
     min E from the solve 0.5369.
   - O1: same counters as O4. OA: 75 face evaluations; OB: 19. |q| was never capped at c.
   - NOT ADMISSIBLE / NOT POSITIVE: 0 in every log (R4 N4 O4 O1 OA OB D4).
   - XSPLIT_BLEND: **none printed** (D4 with dbg_t2_admiss = 50, and the deck's 20 in the other arms).
     - dbg_t2_admiss prints only for a non-admissible hesdirk2 stage (rad_m1_implicit.cpp:17393), for a pass solved to
       E <= 0 with implicit_vimp (16347), or for a solve that ends with E <= e_floor / eint <= 0 (17469).
     - None occurred: stage fallbacks 0, floor clips 0, min E 0.536.
     - The thin-cavity E < 0 case of the 6-node production state does not exist on this small fresh state.
   - Picard passes per step (user.hst col 22), cycles 1-6:
     - O4: 24, 22, 3.5, 3.5, 3, 3;
     - N4: 7, 6.5, 3.5, 3.5, 3, 3.
     - The extra passes are all in the first two steps: the backward-Euler step and the first hesdirk2 step, right after
       the IC relaxation.
5. **Energy (O4 - N4, last row t = 857.14 s, relative to the total).**
   - Etot (user col 19) = 3.7070386010033e46: -2.7e-16. Rows 1-6 are within +-3.3e-15.
   - hydro tot-E: -3.3e-16. E_rad: +1.2e-14. e_gas: -3.6e-15. L_top: -1.2e-9.
   - Floor-energy line: E floor raises cell-solves=0, energy from gas=0, energy created=0; pos_floor_solve undershoot 0;
     implicit_pos_floor_cons take 0. Same in N4.
6. **Wall time** (srun step incl. ~30 s setup; R4/N4 ran 5 steps side by side, O4/D4/OA ran 3 side by side, on the same node):

   | R4 | N4 | O4 | D4 | O1 (1 rank) | OA (3 cyc) | OB (3 cyc) |
   |---|---|---|---|---|---|---|
   | 428 s | 431 s | 650 s | 651 s | 1075 s | 515 s | 170 s |

   - O4 vs N4: +51 % wall (650 vs 431 s), 2.0x Picard passes (94 vs 46 bicgstab outer passes). The extra cost is in the
     first two steps; from cycle 3 the Picard counts are the same.
   - OB (cycles 4-6) took 170 s, so the later ON cycles cost about the same as in N4. That points to a start-up transient
     rather than a per-step cost (inference).
   - The OOM-killed first try of O4 reached the end of cycle 6 in ~680 s, which is consistent.
