# NOTE 2026-11-09 Orion -> viper: hglagt CPU gates for `<rad_m1>/implicit_lagged_t_upwind` (answer to TASK-2026-11-09-orion-hglagt)

Short answer:
- Criteria 1-3: PASS.
- Criterion 4: **the path did not activate on this small fresh state. O4 == N4 bitwise** (hst, bin and the
  double-precision restart payload).
  - The face test never passed: 0 upwinded face-passes in every x1 bin.
  - So criterion 5 is exactly 0, and the CPU cost of the key is nil (O4 350 s vs N4 348 s).
  - The active-path test therefore has to be done on the 6-node production state on viper.
- No XSPLIT_LAGT / XSPLIT_BLEND line was printed: no stage fault and no E <= 0.

## Commits, builds, jobs
- NEW = `7fb7af05a` (origin/hgpos-1109b): detached worktree `/orion/ptmp/jinma/hglagt_1109/wt_new`, kokkos d8e9af03.
  TASK recipe (gcc/13 openmpi/4.1 cmake/3.30, Release, MPI, he_star_m1, `-ffp-contract=off`), `make -j16`: 118 s, 0 warnings.
- REF = `1f314023a`: the hgblend NEW binary, reused (`/orion/ptmp/jinma/hgblend_1109/wt_new/build/src/athena`).
- Deck: `docs/handover/hgpos-1109/` from bsg-files-1009. All 4 MD5SUMS OK; identical to the hgblend deck
  (`layout_invariant = true`).
  - `on.athinput` = deck + `implicit_lagged_t_upwind = true` after `<rad_m1>` (sed; a fresh start cannot add absent keys).
    O4, O1, OA and D4 use it; R4 and N4 use the deck as is.
  - `rad_m1/implicit_lagged_t_upwind=true` was also passed on the command line in every ON run, OB included.
  - Neither implicit_pos_transport_t nor implicit_blend_t_upwind is named, so nothing was superseded.
  - gam = `implicit_pos_transport_mem` is not in the deck, so the default 1/6 applies. e_floor = 1e-4 (deck).
- Keys: `time/nlim=6 output1/dt=1e-30 output2/dt=1e30 output3/dt=1e30 output6/dt=1e30` (OA nlim=3).
  - OB = `-r OA/rst/hegiant.00001.rst` (cycle 3) + nlim=6 + ON.
  - D4 = ON + `rad_m1/dbg_t2_admiss=50`. The deck's own `dbg_t2_admiss = 20` is in every arm.
- Slurm, p.shared `--qos=serial`, two jobs started at the same time on orion001:
  - **227966** = R4 N4 O4 D4 OA->OB (`-n 20 --mem=250G`; each step `srun -n 4 --exact --mem=48G --cpu-bind=cores`).
    MaxRSS 5.80-5.87 GB per task.
  - **227965** = O1 alone (1 rank, `--mem=120G`).
- Every run: 6 cycles, t = 857.14 s (cycle 1 backward Euler, 2-6 hesdirk2), rc 0.
- Run dirs: `/orion/ptmp/jinma/hglagt_1109/runs/<id>` (cmd.txt, log, hst, wall.txt). All rst and bin were deleted after
  the comparisons.
- Scripts: `/orion/ptmp/jinma/hglagt_1109/scripts/` (bincmp.py, rstcmp.py: payload after the header; dev_cs.py as in hgblend).

## Criteria
1. **Keys absent, N4 vs R4: PASS.** `cmp`:
   - hydro.hst and user.hst identical;
   - all bin files (hydro_w, m1, m1_vet at 00000 and 00001) identical as whole files;
   - rst 00000 and 00001 identical as whole files.
2. **Rank layout, O1 vs O4: PASS.** 
   - Data payloads of hydro_w, m1 and m1_vet are identical at 00000 and 00001.
   - The rst payload is identical at 00000 and at the final 00001.
   - hst: time and dt are identical in all 7 rows. The rest differs at the MPI_Reduce round-off level: user.hst max rel
     2.48e-13 (col 9), hydro.hst 1.89e-13. It is there from t = 0, as in hgpos/hgblend.
   - O1 MaxRSS 25.4 GB.
3. **Restart, OA+OB vs O4: PASS.**
   - The last 3 rows of both hst files are identical.
   - The final hydro_w, m1 and m1_vet payloads are identical.
   - OB's final rst (00002) payload is identical to O4's 00001 (2,595,577,318 B).
   - OB: stage steps 3, backward-Euler steps 0.
4. **Active path: NOT ACTIVATED. O4 == N4 bitwise.**
   - `cmp` of both hst files: identical.
   - bin data payloads: identical (hydro_w md5 86500fccd651, m1 ef8b3a9f3c36, m1_vet ea438c218ce2). These are the same
     payloads as the hgpos N4/O4 runs.
   - Final rst payload identical. So max |dE|/E = max |dT|/T = 0, over all cells and at tau > 1.
   - D4 == O4 bitwise as well.
   - Counters (O4; D4 identical):
     ```
     <rad_m1> implicit_lagged_t_upwind (rank 0): x2/x3 face evaluations=4.600000e+01 upwinded face-passes by global x1 bin [0,G/4) [G/4,G/2) [G/2,G): 0 0 0 of which |q| capped at c: 0 0 0
     <rad_m1> time_scheme=hesdirk2: stage steps=5.000000e+00 backward-Euler steps=1.000000e+00 stage fallbacks=0.000000e+00 ...
     <rad_m1> implicit transport: solves=11 Picard iterations mean=4.181818 max=8 NON-CONVERGED=0     (= N4)
     min E from the solve=5.369157e-01 ; floor clips: solved E <= e_floor cell-solves=0, written-back gas eint <= 0 cell-solves=0
     ```
     - OA: 27 face evaluations, 0 upwinded. OB: 19 evaluations, 0 upwinded. O1: 46 evaluations, 0 upwinded.
   - NOT ADMISSIBLE / NOT POSITIVE: 0 in all logs.
   - XSPLIT_LAGT: 0 lines; XSPLIT_BLEND: 0 lines (D4, dbg_t2_admiss = 50).
     dbg_t2_admiss prints only for a non-admissible stage, a pass solved to E <= 0, or a solve that ends at the floor;
     none occurred.
   - **Branch condition that kept it off:** `M1LagTUpQ` (rad_m1_implicit.hpp:1013) returns q = 0 unless
     `(chat/c) |L_f| (A_f dt / V_donor) > gam max(E^n_donor, e_floor)`, with gam = 1/6 (default) and e_floor = 1e-4.
     - On this state no x2/x3 face, seams included, has a lagged remainder whose outflow drains more than 1/6 of the
       donor's E^n in one step.
     - The counters are rank 0 only. The bitwise O4 == N4 shows the same on every rank.
     - The faces were evaluated (46 per run on rank 0); `lft` is nonzero wherever blend/skew terms exist, but under the
       threshold.
   - Contrast with hgblend: implicit_blend_t_upwind has no such threshold. It acts on every blend face, which is why it
     changed the state there (dE/E 3.3e-6).
5. **Energy: O4 - N4 = 0 exactly**, in every column and row.
   - Last row Etot (user col 19) = 3.7070386010033187e46; hydro tot-E = 3.0868082938948414e46 (both = N4).
   - Floor line: E floor raises cell-solves=0, energy created=0, pos_floor_cons take=0.
6. **Wall time** (srun step incl. ~30 s setup; 5 steps side by side on one node, O1 in its own job on the same node):

   | R4 | N4 | O4 | D4 | O1 (1 rank) | OA (3 cyc) | OB (3 cyc) |
   |---|---|---|---|---|---|---|
   | 349 s | 348 s | 350 s | 350 s | 670 s | 229 s | 202 s |

   O4 vs N4: +0.6 %, within noise. The test runs on every face pass, but nothing is upwinded here.
