# TASK 2026-11-09 viper -> Orion: CPU correctness gates for two opt-in rad_m1 keys (hgpos-1109)

Self-contained. CPU only. Answer with `docs/handover/NOTE-2026-11-09-orion-hgpos.md` on fork/bsg-files-1009
(push it, then say so). Do not push to any other branch. Delete your run dirs' restart dumps after the
comparisons (quota).

## What is being gated
Branch **fork/hgpos-1109b** = fork/prod-1112-cand (f34ada236) + 3 commits:
- `15ec17cab` + `2f43a98d6`: `<rad_m1>/time2_adm_eflr` (bool, read only when named, default off). A hesdirk2 stage whose
  only fault is a solved E <= 0 is admitted. Those cells are then taken by the conservative write-back floor, the same
  as the backward-Euler redo. The key needs implicit_pos_floor, implicit_pos_floor_solve and implicit_pos_floor_cons =
  true, which are all on in the deck below.
- `b798f1a30`: `<rad_m1>/implicit_pos_transport_t = upwind` (read only when named, default none). The x2/x3 face memory
  outflow of a dim donor is taken implicitly upwind (the M1PtMemQ treatment of x1, now on transverse/cs faces).
  dbg_t2_admiss gets an XSPLIT line.
Gate commit: **b798f1a30** (the branch tip at writing). Reference commit: **f34ada236** (prod-1112-cand).

## Build (twice: b798f1a30 and f34ada236), as in NOTE-2026-11-08-orion-viper-consfix-ledger
Use detached worktrees (`git worktree add --detach <dir> <sha>`, then `git submodule update --init --recursive`):
```
module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build -j32
```
Call the binaries NEW (b798f1a30) and REF (f34ada236).

## Deck and files (on fork/bsg-files-1009)
`docs/handover/hgpos-1109/`:
- `hg_cs_small.athinput`: the HG production deck (He giant, cubed sphere, M1 implicit, hesdirk2, vet_col, general
  EOS) as embedded in the small HG restart, with these changes:
  - mesh nx2 = nx3 = 16, meshblock 704 x 8 x 8, so 24 MeshBlocks and 6 x 16^2 x 704 = 1.08e6 cells;
  - the table paths point at `files/`;
  - the output last_time / file_number lines are removed.
- `files/` has the IC column, the Rosseland table and the Planck table. `MD5SUMS` checks all four files.

Run from a directory that holds the deck and `files/` (the paths are relative). It is a fresh start: the pgen builds
the IC and the radiation relaxation passes. Expect the first step to be backward Euler and hesdirk2 stages from cycle 2.

Keys for all runs (command line): `mesh/layout_invariant=true time/nlim=6 output1/dt=1e-30 output2/dt=1e30 output3/dt=1e30 output6/dt=1e30`.
output1 (hst) then writes every cycle. ON = `rad_m1/time2_adm_eflr=true rad_m1/implicit_pos_transport_t=upwind`.

## Runs (srun -n <N> --exact --mem=... --cpu-bind=cores <bin> -i hg_cs_small.athinput <keys> -d <dir>)
| id | bin | ranks | extra keys |
|---|---|---|---|
| R4 | REF | 4 | - |
| N4 | NEW | 4 | - |
| O4 | NEW | 4 | ON |
| O1 | NEW | 1 | ON |
| OA | NEW | 4 | ON, `time/nlim=3`; then the final restart it writes in `<dir>/rst/` is used by OB |
| OB | NEW | 4 | `-r <OA rst> time/nlim=6` + ON |
| E4 | NEW | 4 | `rad_m1/time2_adm_eflr=true` only |
| T4 | NEW | 4 | `rad_m1/implicit_pos_transport_t=upwind` only |

If 6 cycles take more than ~20 min per run, use nlim 4 / OA nlim 2 and say so.

## Pass criteria
1. **Keys absent bitwise**: N4 vs R4. Identical `hegiant.hydro.hst` and `hegiant.user.hst` (`cmp`), and identical final
   `bin/*.bin` dumps (the hydro_w and m1 dumps written at the end).
2. **Rank layout**: O1 vs O4. Identical hst files and final dumps. (mesh/layout_invariant = true is on the command line, as in production); if not
   identical, report the max relative difference.
3. **Restart**: OA+OB vs O4. The last 3 hst rows and the final dumps are identical.
4. **Runs**: O4, E4 and T4 cycle with no FATAL and no NaN. Report the end-of-run summary lines:
   - `time_scheme=hesdirk2: stage steps=... backward-Euler steps=...`;
   - `time2_adm_eflr admitted stages=` (printed with the hesdirk2 counters);
   - the implicit_pos_transport line, which counts x1 and transverse faces upwinded.
   Also report any NOT ADMISSIBLE / NOT POSITIVE lines.
5. **Energy** (diagnostic): the total energy column(s) of hegiant.user.hst. Report O4 - N4 at the last row relative to
   the total. Also grep the log for the `implicit_pos_floor_cons` / floor-energy lines and report them.

This small fresh state may well have no stage faults, so the new paths may stay inactive. Say so if the counters are 0:
gates 1-3 still hold (keys absent / layout / restart), and the active-path test is done on the 6-node GPU state on viper.

## NOTE back
`docs/handover/NOTE-2026-11-09-orion-hgpos.md` on fork/bsg-files-1009 should contain:
- the commits, builds and job ids;
- per gate PASS/FAIL with the evidence (cmp output, max differences);
- the summary lines from criterion 4;
- wall time per run.
