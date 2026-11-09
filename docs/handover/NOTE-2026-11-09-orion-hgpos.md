# NOTE 2026-11-09 Orion -> viper: hgpos-1109 CPU gates (answer to TASK-2026-11-09-orion-hgpos)

Short answer: gates 1-4 PASS, with one exception under gate 2. Gate 5 is a diagnostic; the difference is exactly 0.
- Gate 2 exception: the O1 and O4 hst files are not bitwise identical (max rel 2.5e-13).
  The final dumps, including the double-precision restart payload, are bitwise identical.
- On this small fresh state both new paths stay inactive. Every new counter is 0: admitted stages 0, upwinded x2/x3 faces 0.
- O4, E4 and T4 are bitwise identical to N4 (hst, bin, rst payload). So the active-path test is still to be done on the
  6-node GPU state on viper, as the TASK foresaw.

## Commits, builds, deck change, jobs
- NEW = `b798f1a30` (origin/hgpos-1109b tip), REF = `f34ada236` (origin/prod-1112-cand).
  Detached worktrees `/orion/ptmp/jinma/hgpos_1109/wt_{new,ref}`, kokkos d8e9af03.
- Build: the TASK recipe exactly (gcc/13 openmpi/4.1 cmake/3.30, `-ffp-contract=off`, Release, MPI, he_star_m1).
  `make -j16`, both built side by side: about 8 min each, 0 warnings.
- Inputs: `docs/handover/hgpos-1109/*`, all 4 MD5SUMS OK.
- **Deck deviation (no code change).** On a fresh start, AthenaK's `ModifyFromCmdline` refuses a key that is not in the deck:
  `Parameter 'layout_invariant' in block 'mesh' on command line not found` (parameter_input.cpp:430).
  Only a restart adds absent keys. The first job, 227905, failed this way in all arms within seconds; its outputs were deleted.
  None of `mesh/layout_invariant`, `rad_m1/time2_adm_eflr` or `rad_m1/implicit_pos_transport_t` is in `hg_cs_small.athinput`.
  So the keys were inserted into deck copies with sed:
  - `base.athinput` = deck + `layout_invariant = true` in `<mesh>`. Used by R4 and N4, so the new keys are truly absent there.
  - `on.athinput` = base + both rad_m1 keys. Used by O4, O1 and OA. OB reads its input from OA's restart.
  - `e.athinput` = base + `time2_adm_eflr = true`.
  - `t.athinput` = base + `implicit_pos_transport_t = upwind`.

  The TASK's command-line keys were still all passed. They now override existing lines with the same values.
  Without the mesh key, recipe none gives layout_invariant = false (mesh.cpp:197-210). Suggestion: add the line to the
  deck on the branch.
- Command keys for every run: `mesh/layout_invariant=true time/nlim=6 output1/dt=1e-30 output2/dt=1e30 output3/dt=1e30 output6/dt=1e30` (+ arm keys).
  OA uses `time/nlim=3`.
  OB = `-r OA/rst/hegiant.00001.rst` (OA's final restart, cycle 3) + the same keys + `time/nlim=6` + ON.
- Slurm:
  - Job 227907: p.shared, `--qos=serial -n 32 --mem=200G`, one node (orion002). Each arm ran as
    `srun -n N --exact --mem=24G --cpu-bind=cores`, all in parallel, with OB after OA.
  - O1 was OOM-killed at 24G in 227907. The ragged vet_gd band is 4.45 GB on one rank; MaxRSS was 25 GB at the kill.
  - O1 was rerun alone in job 227909 (p.shared, `--mem=120G`, orion003). MaxRSS 27.1 GB.
- Run dirs: `/orion/ptmp/jinma/hgpos_1109/runs/<id>` (cmd.txt, log, hst, wall.txt). All rst/bin were deleted after the comparisons.
- Compare scripts: `/orion/ptmp/jinma/hgpos_1109/scripts/{bincmp.py,rstcmp.py}`. They compare the data payload after the header.
  The headers embed the PAR_DUMP of the deck, which differs between the deck copies by design.
- Cycles: every arm ran 6 cycles (Terminating on cycle limit, t = 857.14 s, dt = 142.86 s): cycle 1 backward Euler,
  cycles 2-6 hesdirk2. The hst has 7 rows (cycles 0-6). No nlim reduction was needed.

## Gates
1. **Keys absent bitwise, N4 vs R4: PASS.**
   - `cmp`: hegiant.hydro.hst and hegiant.user.hst identical.
   - All bin dumps (hydro_w, m1, m1_vet at 00000 and 00001) are byte-identical as whole files. N4 and R4 used the same deck.
   - The final rst (00000 and 00001) are byte-identical as whole files too.
2. **Rank layout, O1 vs O4: PASS for the state, hst not bitwise.**
   - Final dumps: the hydro_w, m1 and m1_vet data payloads are identical at 00000 and 00001. The final restart payload
     (2,595,577,318 bytes, all double-precision state incl. hesdirk2 / vet_gd blocks) is identical, and so is rst 00000.
   - hst: time and dt are identical in all 7 rows. The other columns differ at round-off:
     - user.hst max rel 2.48e-13 (col 9, e_gas);
     - hydro.hst max rel 1.89e-13.
     - They differ already at row 0 (t = 0), where the restart-0 state is bitwise identical.
   - Inference: the hst sums come from `MPI_Reduce` of per-rank partial sums (outputs/history.cpp:415). That reduction is
     not in the layout_invariant (repro) machinery, so its summation order depends on the rank count. This is a diagnostic
     output difference, not a state difference.
3. **Restart, OA+OB vs O4: PASS.**
   - The last 3 rows of hydro.hst and of user.hst are identical (`diff` of `tail -n 3`).
   - The final hydro_w, m1 and m1_vet payloads are identical.
   - OB's final rst (00002) payload is identical to O4's final rst (00001).
   - OB counters: stage steps 3, backward-Euler steps 0. The restart carries the hesdirk2 slope.
4. **Runs: PASS.** O4, E4, T4 (and O1, OA, OB): rc 0, no FATAL, no NaN, 0 NOT ADMISSIBLE / NOT POSITIVE lines.
   ```
   O4: <rad_m1> time_scheme=hesdirk2: stage steps=5.000000e+00 backward-Euler steps=1.000000e+00 stage fallbacks=0.000000e+00 ... accepted-with-fallback stages=0.000000e+00 time2_adm_eflr admitted stages=0.000000e+00
   O4: <rad_m1> implicit_pos_transport_t = upwind (rank 0): memory-upwinded x2/x3 faces (first passes)=0.000000e+00 in solves=1.100000e+01
   E4: <rad_m1> time_scheme=hesdirk2: stage steps=5.000000e+00 backward-Euler steps=1.000000e+00 ... accepted-with-fallback stages=0.000000e+00 time2_adm_eflr admitted stages=0.000000e+00
   T4: <rad_m1> time_scheme=hesdirk2: stage steps=5.000000e+00 backward-Euler steps=1.000000e+00 ... accepted-with-fallback stages=0.000000e+00
   T4: <rad_m1> implicit_pos_transport_t = upwind (rank 0): memory-upwinded x2/x3 faces (first passes)=0.000000e+00 in solves=1.100000e+01
   O1: stage steps=5 backward-Euler steps=1 admitted stages=0; x2/x3 faces=0 in solves=11
   OA: stage steps=2 backward-Euler steps=1 admitted stages=0; x2/x3 faces=0 in solves=5
   OB: stage steps=3 backward-Euler steps=0 admitted stages=0; x2/x3 faces=0 in solves=6
   ```
   - The x1 `memory-upwinded x1 face passes` line is not printed. The deck has no `<rad_m1>/implicit_pos_transport`, so the
     x1 path is off, and the transport line counts only transverse faces.
   - The `time2_adm_eflr admitted stages=` field appears only when the key is named (it is absent in the N4, R4 and T4 lines).
   - No stage faults: stage fallbacks 0, time2_auto_be skipped-stage BE 0, positivity fallbacks 0, min E from the solve 0.537.
   - **So both new paths were inactive on this state.**
5. **Energy (diagnostic).** O4 - N4 = 0 exactly, in every column and every row of both hst files (bitwise identical).
   - Last row (t = 857.14 s): Etot (user col 19) = 3.7070386010033187e46, E_rad = 6.2023e45, e_gas = 3.4886e45;
     hydro tot-E = 3.0868082938948414e46.
   - Floor lines (O4 = N4 = E4 = T4):
     - `E floor raises cell-solves=0 energy from gas=0 energy created=0 ... solve undershoot (pos_floor_solve)=0 taken from the column's radiation (implicit_pos_floor_cons)=0`;
     - `floor clips (all ranks, every solve incl. BE): solved E <= e_floor cell-solves=0 written-back gas eint <= 0 cell-solves=0`;
     - g0 clipped cell-passes = 1.578161e6 in every arm.

## Wall time per run (srun step, incl. start-up; cycle 0 setup ~29 s)
| R4 | N4 | O4 | O1 | OA (3 cyc) | OB (3 cyc) | E4 | T4 |
|---|---|---|---|---|---|---|---|
| 341 s | 344 s | 345 s | 662 s (1 rank) | 227 s | 159 s | 344 s | 346 s |

The new keys cost no measurable time here (O4 345 s vs N4 344 s); they did no work.
