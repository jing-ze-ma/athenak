# TASK 2026-11-09 viper -> Orion: CPU gates for `<rad_m1>/implicit_lagged_t_upwind` (hgpos-1109b)

Self-contained; CPU only; same deck, build recipe and run layout as `TASK-2026-11-09-orion-hgpos.md` (read that first;
your NOTE-2026-11-09-orion-hgpos 35331b94a shows it works). Answer with `docs/handover/NOTE-2026-11-09-orion-hgfct.md`
on fork/bsg-files-1009 (push, say so). Delete restart dumps after the comparisons.

## What is gated
fork/hgpos-1109b **a4bac617e** = 7fb7af05a (hglagt) + `<rad_m1>/implicit_lagged_t_fct = alpha` (Real, read only when
named, off when unset; cs + sp; FATAL with the older transverse keys): a Zalesak/FCT donor limiter on the lagged x2/x3
remainder L_f of the face flux (no matrix change): on the first Picard pass of each solve each cell's lagged outflow
Q_out is summed, R = min(1, alpha max(EN,0)/Q_out) is frozen for the solve, and every face's L_f is multiplied by the
donor's R (one number per face, conservative). ImplicitReport prints its counters (limited faces/cells per x1 bin,
removed outflow, RHS-sign counters).
Reference: **7fb7af05a** (your hglagt NEW binary).

## Deck
`docs/handover/hgpos-1109/hg_cs_small.athinput` + `files/` (now carries `layout_invariant = true`; MD5SUMS updated).
Keys for all runs: `time/nlim=6 output1/dt=1e-30 output2/dt=1e30 output3/dt=1e30 output6/dt=1e30`.
ON = `rad_m1/implicit_lagged_t_fct=0.001`.

## Runs (4 ranks unless stated; srun ... -i hg_cs_small.athinput <keys> -d <dir>)
| id | bin | ranks | extra |
|---|---|---|---|
| R4 | REF (7fb7af05a) | 4 | - |
| N4 | NEW (a4bac617e) | 4 | - |
| O4 | NEW | 4 | ON |
| O1 | NEW | 1 | ON |
| OA | NEW | 4 | ON, nlim=3 (final rst -> OB) |
| OB | NEW | 4 | -r <OA rst> nlim=6, ON |
| D4 | NEW | 4 | ON + `rad_m1/dbg_t2_admiss=50` (diagnostic only) |

## Pass criteria
1. Keys absent bitwise: N4 vs R4, identical hst files and final dumps.
2. Rank layout: O1 vs O4 identical state (final dumps); hst identical or round-off only (the MPI_Reduce you found).
3. Restart: OA+OB vs O4, last 3 hst rows + final dumps identical.
4. **Active path** (alpha = 1e-3 forces the limiter to bind; the gate is that it ACTS, not accuracy): O4 vs N4 must DIFFER, the FCT counters must be > 0, total energy O4 vs N4 at round-off (conservative limiter), E > 0 everywhere; (the key acts wherever a transverse blend outflow exists; the deck runs
   implicit_flux = blend with thin outer layers). Report max |dE|/E and max |dT|/T of the final dumps over all cells and
   over tau > 1 (as in your consfix NOTE), the end-of-run counters (`time_scheme=hesdirk2: stage steps / backward-Euler
   steps`, any implicit_lagged_t_upwind line), NOT ADMISSIBLE / NOT POSITIVE counts, and the XSPLIT_BLEND lines of D4 if
   any were printed. If O4 == N4 bitwise, say so: the path did not activate (then report which branch condition kept it off,
   e.g. no blend faces in x2/x3).
5. Energy: total-energy columns of hegiant.user.hst, O4 - N4 relative at the last row.
6. Wall time per run (O4 vs N4: the cost of the key on CPU, informational only; GPU cost is measured on viper).

## NOTE back
`docs/handover/NOTE-2026-11-09-orion-hgfct.md` on fork/bsg-files-1009: commits/builds/job ids, PASS/FAIL per criterion
with evidence, the numbers of criterion 4-6.
