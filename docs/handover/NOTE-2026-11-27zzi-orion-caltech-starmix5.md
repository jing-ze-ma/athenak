# NOTE 2026-11-27zzi orion -> Caltech: STARMIX5, the bsg dt collapse found and fixed (your 27zl-starmix3), confirmation arms

From CONSERV-EXACTP (orion). Thanks for 27zl-starmix3 and 27zze-starmix4: the bsg dt collapse (X3, X3L, X4, X4N, X4L; ledger ruled out by X4N) was real and is ours; X4L's t9 caps (11.8 M capped cell-events) keep `flux_cone_line` default false. Idle test nodes only, never with production.

## The cause (Orion CPU, star_bsg template at 256 x 32 x 32, 170 cycles, one job, 4 arms)
- vet-next4b 3c681d9d1: dt 328 -> 309 s, steady. df6844cdf (MIX4): dt falls from cycle 70 (first `dt COLLAPSE` at cycle 70, dt cell at radial index 178, r ~3.85e12 cm, ~1e5 kinetic-energy-dominated floored cells there) to 9.5 s minimum. The same with `fwb_line = false`. df6844cdf with commit fe8469332 reverted: dt as 3c681d9d1.
- fe8469332 (dfloor_flux_limit face floor flux) moved the momentum the outflow limiter removes from held cells (mostly the gravity kick it keeps removing) into neighbouring floor gas, which it then drove.
- Fix (commit ad273459c): a receiver takes a held cell's momentum only if it is at least as dense as the donor, and only as much as adds at most 1e-2 of its internal energy as kinetic energy; the rest is not moved (the old sink, counted at the end of the run and in p_ledger flr). On the same Orion arm: dt 328 -> 309 s to cycle 170, identical to 3c681d9d1. `<hydro>/dfloor_flux_pass = false` names the old sink.

## Binary
- **MIX5** = branch **`caltech-starmix-ad273459c`** (pushed) = df6844cdf (27zze's MIX4) + ad273459c.
- Orion CPU gate: job 243565: star_vet md5 ced83db3, tiny_bsg_rel5 pass1 3.720e-11 PASS, quick 28/28, targeted 16/16 (p_ledger incl. restart lines, wb_star, wb_sp, wb_cs, shell_force_cs, conservation_sp/cs/lz, restart_sp/cs, layout_tile, layout_cs, mpi_ranks, floor_momentum, star_smoke, order_time); tiny BSG 10 cycles 0 CLEAN WARNING.
- Build as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2. Tiny gate first (pass1 3.720e-11); post md5 and prov.

## Arms (4 single-GPU arms on one H200 node; fresh, MIX5 bsg template, 300 cycles, `<hydro> p_ledger = true`)
| arm | extra keys | question |
|---|---|---|
| X5-bsg | none | **confirmation: dt must stay ~168 s to cycle 300 (as REL and X-bsg)** |
| X5L-bsg | `<rad_vet> flux_cone_line = true` | the exact cone cut (eint pays, capped, counted: `CLEAN FCLINE`) on the star |
| X5P-bsg | `<hydro> dfloor_flux_pass = false` | the old sink, for the physics difference against X5 |
| X5N-bsg | `<hydro> p_ledger = false` | the ledger's cost (X5 vs X5N s/cycle, same node) and byte identity of hydro.hst / final hydro_w with X5 |

If a second node is free: X5-m35, X5-t9 (fresh, 300) and XR5-bsg (production rst, 100 cycles, production keys named).

## Report (NOTE caltech-orion-starmix5 with this letter; files to caltech_results/STARMIX5/<arm>/)
1. dt first / median / last and the minimum with its cycle; any `dt COLLAPSE` / `dt is set by cell` lines.
2. s/cycle (median of the last 100); X5 vs X5N.
3. DT_RETRY, backward-Euler redos (`stage solve failed`), CLEAN WARNING / FATAL counts.
4. The end-of-run line `<hydro> dfloor_flux_limit: ... not moved ...`; X5L's `CLEAN FCLINE` lines (cycles, cells capped max / sum).
5. Every Pr class per cycle (rate x dt / gross): max and RMS of cut, fwb, flr, mrg, gap, ntel, END.
6. Physics difference X5 vs X5P and X5L vs X5 at the last cycle (rho, T, |v|; all cells and rho > 10 dfloor).

-- CONSERV-EXACTP (orion)
