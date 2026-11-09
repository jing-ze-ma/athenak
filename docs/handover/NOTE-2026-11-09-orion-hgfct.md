# NOTE 2026-11-09 Orion -> viper: hgfct CPU gates for `<rad_m1>/implicit_lagged_t_fct` (answer to TASK-2026-11-09-orion-hgfct)

Short answer:
- Criteria 1-3: PASS (keys absent bitwise; O1/O4 state bitwise; restart bitwise).
- Criterion 4: the limiter ACTS.
  - O4 differs from N4 in 91 % of the cells. FCT counters > 0. E > 0 everywhere (min E 0.556, min E from the solve 0.513).
  - No NOT ADMISSIBLE / NOT POSITIVE, no fallbacks.
- **Energy is not at round-off:** Etot(O4) - Etot(N4) = -1.43e-10 of the total at the last row.
  - The difference is the extra radiation that leaves through the top boundary. O4's L_top is larger by 3.5e-5.
  - In step 1 (backward Euler) dEtot equals -dL_top dt to 1.0e-4. Over 6 steps -5.29e36 erg vs int dL_top dt =
    4.57e36 (trapezoid) .. 4.96e36 (right sum) erg, within the hst quadrature.
  - So the limiter itself does not appear to create or destroy energy (inference; a `dbg_cons_ledger` run would settle it).
  - The hydro tot-E differs by -6.4e-13.
- Cost on CPU: nil (O4 371 s vs N4 366 s; Picard counts identical).

## Commits, builds, jobs
- NEW = `a4bac617e` (origin/hgpos-1109b): detached worktree `/orion/ptmp/jinma/hgfct_1109/wt_new`, kokkos d8e9af03.
  TASK recipe (gcc/13 openmpi/4.1 cmake/3.30, Release, MPI, he_star_m1, `-ffp-contract=off`), `make -j16`: 122 s, 0 warnings.
- REF = `7fb7af05a`: the hglagt NEW binary, reused read-only (`/orion/ptmp/jinma/hglagt_1109/wt_new/build/src/athena`).
- Deck: `docs/handover/hgpos-1109/` from bsg-files-1009. All 4 MD5SUMS OK; same as hglagt.
  - `on.athinput` = deck + `implicit_lagged_t_fct = 0.001` after `<rad_m1>` (sed; a fresh start cannot add absent keys).
    O4, O1, OA and D4 use it; R4 and N4 use the deck as is.
  - `rad_m1/implicit_lagged_t_fct=0.001` was also passed on the command line in every ON run, OB included.
- Keys: `time/nlim=6 output1/dt=1e-30 output2/dt=1e30 output3/dt=1e30 output6/dt=1e30` (OA nlim=3).
  - OB = `-r OA/rst/hegiant.00001.rst` (cycle 3) + nlim=6 + ON.
  - D4 = ON + `rad_m1/dbg_t2_admiss=50`. The deck's own `dbg_t2_admiss = 20` is in every arm.
- Slurm, p.shared `--qos=serial`, two jobs started at the same time on orion001:
  - **228046** = R4 N4 O4 D4 OA->OB (`-n 20 --mem=250G`; each step `srun -n 4 --exact --mem=48G --cpu-bind=cores`).
    MaxRSS 5.6-6.3 GB per task.
  - **228045** = O1 alone (1 rank, `--mem=120G`). MaxRSS 25.4 GB.
- Every run: 6 cycles, t = 857.14 s (cycle 1 backward Euler, 2-6 hesdirk2), rc 0, no FATAL / NaN / OOM.
- Run dirs: `/orion/ptmp/jinma/hgfct_1109/runs/<id>` (cmd.txt, log, hst, wall.txt). All rst and bin were deleted after
  the comparisons.
- Scripts: `/orion/ptmp/jinma/hgfct_1109/scripts/`.
  - bincmp.py and rstcmp.py compare the payload after the header.
  - dev_cs.py: E = m1 u0, T = wtemp; tau from the deck's Rosseland table on N4's state; T/T_rad = 1.0000 at tau > 3;
    445,440 of 1,081,344 cells at tau > 1.

## Criteria
1. **Keys absent, N4 vs R4: PASS.** `cmp`:
   - hydro.hst and user.hst identical;
   - all bin files (hydro_w, m1, m1_vet at 00000 and 00001) identical as whole files;
   - rst 00000 and 00001 identical as whole files.
2. **Rank layout, O1 vs O4: PASS.**
   - The data payloads of hydro_w, m1 and m1_vet are identical at 00000 and 00001.
   - The rst payload is identical at 00000 and at the final 00001 (2,595,577,318 B).
   - hst: time and dt are identical in all 7 rows. The rest differs at the MPI_Reduce round-off level: user.hst max rel
     2.64e-13 (col 9), hydro.hst 1.91e-13. It is there from t = 0, as before.
   - D4 vs O4: hst and the final rst payload are bitwise identical.
3. **Restart, OA+OB vs O4: PASS.**
   - The last 3 rows of both hst files are identical.
   - The final hydro_w, m1 and m1_vet payloads are identical.
   - OB's final rst (00002) payload is identical to O4's 00001.
   - OB: stage steps 3, backward-Euler steps 0.
4. **Active path: PASS (it acts), with the energy caveat in 5.**
   - O4 != N4: the hst differs from row 2 (cycle 1) on, and the final rst payload differs. 984,174 of 1,081,344 cells
     have a different E.

   | final state, O4 vs N4 | max \|dE\|/E | max \|dT\|/T |
   |---|---|---|
   | all cells | 1.980e-01 | 4.621e-02 |
   | tau > 1 | 7.020e-05 | 1.755e-05 |
   | tau < 1 | 1.980e-01 | 4.621e-02 |

   - The largest change is at the outer boundary cell: MeshBlock 12, i = 703, rho = 1e-14 (floor), T = 3072 K, tau = 3e-7.
     alpha = 1e-3 forces the limiter to bind hard in the outer layers.
   - E > 0 everywhere: min E of the final state is 0.556 (N4 0.596); min E from the solve is 0.513 (N4 0.537).
     Floor clips 0; E floor raises 0.
   - FCT counter lines (O4, rank 0; D4 identical):
     ```
     <rad_m1> implicit_lagged_t_fct=1.000000e-03 (rank 0): solves=1.100000e+01 x2/x3 face evaluations=4.600000e+01 | limited face-passes (R_donor<1) by global x1 bin [0,G/4) [G/4,G/2) [G/2,G): 0.000000e+00 2.240515e+06 1.366352e+07 removed drain (1-R)(chat/c)|L|A dt/V: 0.000000e+00 9.799095e+04 3.151015e+05 | cells R<1 (per solve, summed): 0.000000e+00 2.152180e+05 1.361718e+06 of which E^n<=0 (R=0): 0.000000e+00
     <rad_m1> implicit_lagged_t_fct RHS signs (rank 0, cell-passes=1.243546e+07): SRCR<0 0.000000e+00 E^n+SRCR<0 0.000000e+00 (1-alpha)max(E^n,0)+SRCR<0 0.000000e+00 JRHS<0 1.183407e+07 E^n+SRCR+JRHS<0 0.000000e+00
     ```
     - OA: limited face-passes 0 / 1.32e6 / 8.03e6; removed drain 0 / 6.3e4 / 2.0e5.
     - OB: 0 / 9.17e5 / 5.63e6; removed drain 0 / 3.5e4 / 1.1e5; cells R<1 0 / 1.17e5 / 7.41e5.
     - O1 (rank 0 = all 24 MeshBlocks): limited face-passes 0 / 8.97e6 / 5.47e7.
     - No limiting at all in the inner quarter of the radial grid.
   - Run counters:
     - `time_scheme=hesdirk2: stage steps=5 backward-Euler steps=1 stage fallbacks=0`, accepted-with-fallback stages 0.
     - `implicit transport: solves=11 Picard iterations mean=4.181818 max=8 NON-CONVERGED=0`. These are the same as N4,
       and Picard per step (user col 22) is identical: 7, 6.5, 3.5, 3.5, 3, 3.
     - No implicit_lagged_t_upwind line (the key is not named).
   - NOT ADMISSIBLE / NOT POSITIVE: 0 in all logs (R4 N4 O4 O1 OA OB D4).
   - D4: no XSPLIT / XSPLIT_LAGF / XSPLIT_BLEND line was printed. dbg_t2_admiss prints only for a non-admissible stage, a
     pass solved to E <= 0, or a solve that ends at the floor; none occurred.
5. **Energy: Etot not at round-off; the difference is top-boundary outflow (inference).**
   - user.hst, O4 - N4, relative, last row (t = 857.14 s):
     - Etot (col 19, 3.7070386e46): **-1.427e-10**;
     - E_rad -8.50e-10; e_gas -5.6e-12; L_top +3.466e-5; L_int +1.28e-7; L_in 0.
     - hydro tot-E: -6.41e-13.
   - Etot relative per row: 0, -3.6e-11, -6.9e-11, -8.8e-11, -9.9e-11, -1.18e-10, -1.43e-10.
   - Budget check, from the hst only:

     | step | d(Etot O4-N4) [erg] | -dL_top dt (end-of-step L_top) [erg] |
     |---|---|---|
     | 1 (BE) | -1.3215e36 | -1.3213e36 (ratio 1.0001) |
     | 2-6 | -1.24e36, -7.0e35, -3.9e35, -7.2e35, -9.1e35 | -1.04e36, -4.0e35, -4.7e35, -9.5e35, -7.9e35 |
     | sum 1-6 | -5.29e36 | -4.96e36 (right sum), -4.57e36 (trapezoid) |

     - For the single-stage BE step, the hst L_top is that step's boundary flux, and the budget closes to 1e-4.
     - For hesdirk2 steps the hst gives only the end-of-step L_top, not the stage fluxes. The per-step scatter of
       +-25-75 % is consistent with that.
     - Other boundary terms are tiny: dL_bot dt ~ 1e25 erg; Mdot_top changes by ~1e14 g/s.
   - Inference: the limiter keeps more radiation in the outer cells' lagged outflow balance, so more of it leaves through
     the top; it is not an internal source. A per-solve `dbg_cons_ledger` run (as in consfix-1109) would confirm it.
     I did not run it (not in the TASK).
   - Floor line: E floor raises cell-solves=0, energy created=0, pos_floor_cons take=0.
6. **Wall time** (srun step incl. ~30 s setup; 5 steps side by side on one node, O1 in its own job on the same node):

   | R4 | N4 | O4 | D4 | O1 (1 rank) | OA (3 cyc) | OB (3 cyc) |
   |---|---|---|---|---|---|---|
   | 366 s | 366 s | 371 s | 371 s | 688 s | 239 s | 167 s |

   - O4 vs N4: +1.4 % (371 vs 366 s). Picard and bicgstab counts are unchanged, so the cost is the extra first-pass
     Q_out sum and the E^n exchange.
   - MaxRSS per task: 6.32 GB (R4, N4) vs 5.89 GB (O4, D4).
