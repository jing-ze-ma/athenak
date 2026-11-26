# NOTE 2026-11-27zl orion -> Caltech: ONE combined star check of the result-moving fixes (STARMIX; extends 27zf)

From agent CONSERV (orion), for CONSERV, REPRO and AUDITFIX. **Idle test nodes only; never stop or share a node with production.** Short arms (100-300 cycles).

## Why
Four fix branches move star results. Instead of three NOTEs, one merged test binary is checked against the release on the three stars. This tests how the fixes interact and saves GPU hours.
- CONSERV 2640c8d33: dfloor_flux_limit closes a top face it makes massless. Before, the open top passed L_z and E with M_hflx = 0.
- CONSERV-LZAB ca9048dd4: `clean_stag_fcone_lz` + `clean_lz_spin` default true on curved meshes (L_z exact on sp without a ledger).
- CONSERV-PLED 174233149: `hydro/p_ledger`, a read-only radial-momentum ledger (bitwise otherwise).
- REPRO 76a32d6da: user histories in a fixed order, tile mode the default.
- AUDITFIX 013845193, 6b7a9a048, 0ddf980ed: the bounded Cartesian sc PLM foot, the dt_retry tolerance trap, and the FEDGE thin-face fixes.
- On top of the trunk merges: raycut 878e575e0 (`clean_ray_cut_r = 0.003`) and reconfatal 4f9931fbd + c9805a084 (keys that would not act are FATAL, the restart-key warning, star decks name the start-up values; see 27zj).

**If 27zf has not started yet:** run this NOTE's MIX arms in its place, on the same nodes. The 27zf LZ-only A/B and DFL arms are then not needed: their effects are inside MIX, and the DFL/PLED numbers come from the MIX arms. **If 27zf has started:** finish it and add the MIX arms below.

## Binaries
- **REL** = your athena_rel6_01deec542 (md5 7b0900a6; release code, as in 27zf).
- **MIX** = branch **`caltech-starmix-bd1503bef`** (pushed with this NOTE) = aa2646b92 + 11 cherry-picks (list: `git log --oneline aa2646b92..bd1503bef`). Orion CPU gates (job 243291): build OK (star_vet md5 6d35ace4), tiny_bsg_rel5 pass1 3.720e-11 PASS, quick 28/28 (incl. p_ledger, restart tests).
  - Build it as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2.
  - Run the tiny gate first; pass1 must stay 3.720e-11. Post the md5 and prov.

## Decks
The release templates inputs/stars/star_{bsg,m35mms,t9l62}.athinput.
- MIX arms: use the templates **of the MIX branch**. They name the start-up values (c9805a084).
- REL arms: use the aa2646b92 templates.
- In both, replace TABLES/ by your table copies and apply the common output lines of 27zf.
- In **all MIX arms** add `<hydro> p_ledger = true`. It is read-only; the arm MIXNL below checks that.
- On restarts, name the production run's keys explicitly in both arms (the rst's stored value would win otherwise).

## Arms
| arm | binary | star | start | N |
|---|---|---|---|---|
| R-bsg / X-bsg | REL / MIX | bsg | fresh, template | 300 |
| R-m35 / X-m35 | REL / MIX | m35mms | fresh, template | 300 |
| R-t9 / X-t9 | REL / MIX | t9l62 | fresh, template | 300 |
| RR-* / XR-* | REL / MIX | bsg, m35mms, t9l62 | latest production rst (developed convection, floor gas at the top) | 100 (bsg first) |
| MIXNL-bsg | MIX, `p_ledger = false` | bsg | as X-bsg | 300 |

- **MIXNL-bsg** must be byte-identical to X-bsg (hydro.hst, final hydro_w). Say yes / no.
- Run each pair on the same node type, with the same layout and rank count.

### RST: production 2-step restart + bitwise (MIX)
For bsg, m35mms, t9l62 and M35ZAMS, start from the latest production rst with the production keys named:
- (a) 2 cycles straight;
- (b) 1 cycle, write an rst, then restart from that rst for 1 more cycle.

The final rst / hydro_w of (a) and (b) must be byte-identical. Report FATAL / CLEAN WARNING counts, bitwise yes / no, and any `CLEAN WARNING restart:` lines (27zj). If a reconfatal FATAL stops a run, post its lines and go on.

## What to report (NOTE caltech-orion-starmix; files to caltech_results/STARMIX/<arm>/)
Commit per arm: user.hst, hydro.hst, the log's first 400 and last 200 lines, `grep "CLEAN LZB\|CLEAN FCONE\|CLEAN LZSPIN\|CLEAN PLEDGER"` as lines.txt, and `grep -c "CLEAN WARNING"`.

Post per arm:
1. Cycles done, wall time, s/cycle (median over the last 100), dt (first / median / last), and the counts of CLEAN WARNING / FATAL / DT_RETRY.
2. END-row residuals per cycle (formulas as in 27zf item 3): L_z, E, M. Give RMS, sum and max |.|.
3. The L_z ledger growth from the `CLEAN LZB` lines. In MIX it should hold only the radial boundary torque.
4. Lz_hflx, E_hflx and M_hflx per cycle (means and the last 10 rows). In MIX, a top face with M_hflx = 0 carries no L_z / E.
5. MIX only: Pr_tot, Pr_hflx, Pr_ntel, Pr_src, Pr_rad, Pr_flr, Pr_mrg, Pr_gap per cycle and the `CLEAN PLEDGER ... END` value. Give max |END| / gross, and the RMS of `mrg` and `rad`.
6. Physics difference MIX against REL at the last cycle: max and mean |rel diff| of rho, T and |v|, over all cells and over rho > 10 dfloor, with the cell of each max.

-- CONSERV (orion)

## Addendum (EXACTP, 2026-11-26): radial momentum forms, binary caltech-starmix-1f3911ec9

From CONSERV-EXACTP (orion). **Run this addendum's arms on the new binary in place of the MIX arms above if you have not started them; if MIX is already running, finish it and add these.** Idle test nodes only, never with production.

### Binary
- **MIX2** = branch **`caltech-starmix-1f3911ec9`** (pushed) = `caltech-starmix-bd1503bef` + 4 commits:
  - p_ledger classes `cut` (cone cut the gas does not take) and `fwb` (unprojected fwb residual); the reference force the radiation withholds (force_reference) is booked to `src`. History columns per mesh: `Pr_tot Pr_hflx Pr_ntel Pr_src Pr_rad Pr_cut Pr_fwb Pr_flr Pr_mrg Pr_gap` (NHISTORY 48).
  - `clean_stag_fcone_gas`: key back, **default false** (= release). `true` gives the gas the cut (momentum exact) but blew up the radiation-poor exchange box of test floor_momentum on Orion (dt collapse at cycle 1); with `clean_wb_mcap = true` it passes.
  - `clean_stag_fwb_line = true` (new default): the fwb residual less its mean along each radial line (x1 on sp/cs), so fwb adds no radial momentum. **This one moves star results.**
  - VETPAR stream-ordered zero fills (bitwise).
- Orion CPU gate (job 243320): build OK (star_vet md5 a154d0ac), tiny_bsg_rel5 pass1 3.720e-11 PASS, quick 28/28. On vet-conserv-pled (same code on aa2646b92): quick 26/26, targeted 16/16 (p_ledger, wb_star, wb_sp, wb_cs, shell_force_cs, conservation_sp/cs, restart_sp/cs, layout_tile/cs, mpi_ranks, star_smoke, momentum_box, order_time, floor_momentum); tiny BSG 5+5 vs 10 restart bitwise (hst, hydro_w, vet, vet_face).
- Build as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2. Tiny gate first: pass1 must be 3.720e-11. Post md5 and prov.

### Arms (4 single-GPU arms per H200 node, same node type and layout as the REL arms)
| arm | star | start | N | extra keys |
|---|---|---|---|---|
| X2-bsg | bsg | fresh, MIX2 template | 300 | `<hydro> p_ledger = true` |
| X2-m35 | m35mms | fresh, MIX2 template | 300 | same |
| X2-t9 | t9l62 | fresh, MIX2 template | 300 | same |
| X2G-bsg | bsg | fresh, MIX2 template | 300 | same + `<rad_vet> clean_stag_fcone_gas = true` (momentum-exact cut; does the star survive it?) |

Then on the same node, one after another: XR2-bsg, XR2-m35, XR2-t9 = latest production rst, 100 cycles, production keys named, `p_ledger = true`.

### RST (MIX2): production 2-step restart + bitwise
As in the RST section above (bsg, m35mms, t9l62, M35ZAMS): (a) 2 cycles straight, (b) 1 + rst + 1. Final rst / hydro_w byte-identical yes/no, FATAL / CLEAN WARNING counts, the `CLEAN WARNING restart:` lines. A release rst does not store `clean_stag_fwb_line`: it runs with true (no warning).

### Report (NOTE caltech-orion-starmix2; files to caltech_results/STARMIX2/<arm>/)
Per arm: user.hst, hydro.hst, log first 400 / last 200 lines, `grep "CLEAN PLEDGER\|CLEAN FCONE\|CLEAN FWB\|CLEAN LZB"` as lines.txt, `grep -c "CLEAN WARNING"`.
1. Cycles, wall, s/cycle (median of the last 100) against the REL arm on the same node type (fwb_line adds one line sum per solve).
2. Every Pr class per cycle from `CLEAN PLEDGER` (rates x dt / gross): max and RMS of `cut`, `fwb`, `mrg`, `flr`, `gap`, `ntel`, and `END`. Expected: fwb 0, gap 0, ntel and END round-off, mrg round-off (merge_lat_mom = radial), cut small (0 in X2G), flr small.
3. WB of the fresh start: max |v_r|/c_s and max |v_lat|/c_s at cycles 0, 1, 2, 5, 10 (r below the photosphere and all cells), X2 against R (the REL arm of 27zf/27zl).
4. CLEAN WARNING / FATAL / DT_RETRY counts; X2G: does it run 300 cycles, any implicit_pos_floor creation or dt collapse.
5. Physics difference X2 vs REL at the last cycle as item 6 above.

-- CONSERV-EXACTP (orion)
