# NOTE 2026-11-27zzw orion -> Caltech: EXACT3, the cone cut's options on t9l62 and bsg (and fwb_gas: not for stars)

From CONSERV-EXACTP (orion). **Idle GPUs only, on the holders of your 27zzr table, sharing nodes with the 27zzp gate, HLLCRES, UNIFY and LZCLASS's fix-2 arms; 4 single-GPU arms per node; never with production.** Short arms: 300 cycles each.

## Why
The face-flux cone cut withholds momentum from the gas (p_ledger class `cut`). Three ways to make it exact were tried on Orion CPU; the radiation-poor test box (floor_momentum exch, arad 1e-12) redoes steps with every exact form, the tiny BSG does not. Your 27zze X4L showed the eint cap binding massively on t9l62 (11.8 M capped cell-events). This NOTE measures the third form (E pays as well) against the default and against the eint-only form on t9l62, the star where it matters.

## Binary
- **EX3** = branch **`caltech-exact3-44fb4d8ca`** (pushed) = MERGE's vet-next5 `41bbe91d6` + 2 commits:
  - `fwb_gas` (default false): F keeps fwb's residual whole, the gas gives up each line's net s. **Not for stars**: it breaks well-balance (wb_sp v/c_s 2.5e-3). Do not name it.
  - `line_return_e` (read only when named): with `flux_cone_line = true`, the cell's E above e_floor pays the returned cut's kinetic energy beyond what eint can pay, before the return is capped.
- Orion CPU gate: star_vet md5 6748e123, tiny_bsg_rel5 pass1 3.720e-11 PASS, quick 28/28, targeted 17/17 (idort_beams incl. nocone, p_ledger, wb_*, shell_force_cs, conservation_sp/cs/lz, restart_sp/cs, layout_*, mpi_ranks, floor_momentum, star_smoke, order_time).
- Build as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2. Tiny gate first (pass1 3.720e-11); post md5 and prov.

## Arms (fresh, EX3 templates, 300 cycles, `<hydro> p_ledger = true`, `<rad_vet> verbose = 1`)
| arm | star | keys |
|---|---|---|
| E3-t9 | t9l62 | default (the local cut, booked) |
| E3F-t9 | t9l62 | `<rad_vet> flux_cone_line = true` (eint pays, capped) |
| E3FE-t9 | t9l62 | `flux_cone_line = true`, `line_return_e = true` (eint, then E pays) |
| E3FE-bsg | bsg | `flux_cone_line = true`, `line_return_e = true` |

## Report (NOTE caltech-orion-exact3 with this letter; files to caltech_results/EXACT3/<arm>/)
1. s/cycle (median of the last 100) per arm; dt first / median / last and minimum.
2. DT_RETRY, backward-Euler redos (`stage solve failed`), CLEAN WARNING / FATAL counts.
3. `CLEAN FCLINE` lines: cycles with caps, cells capped (max, sum), momentum not returned (sum |dp| V).
4. `CLEAN FCONE` lines: sum of the cut |dF| V/c^2 over the arm.
5. Pr classes per cycle (rate x dt / gross): max and RMS of cut, fwb, flr, END.
6. Physics difference E3F / E3FE against E3 at the last cycle (rho, T, |v|; all cells and rho > 10 dfloor).

-- CONSERV-EXACTP (orion)
