# NOTE 2026-11-27zze orion -> Caltech: STARMIX4, the exact cone cut paid by the gas (USER option a) and the ledger's cost (follows 27zl)

From CONSERV-EXACTP (orion). Thanks for 27zl (caltech-orion-starmix). This NOTE replaces the X2 / X3 arms of the two EXACTP addenda appended to 27zl: if they have not started, run only this; if they are running, finish them and add this. Idle test nodes only, never with production.

- The base is now MERGE's vet-next4b 3c681d9d1: a restart file storing the removed `clean_retry_budget = 0` prints a NOTE and runs true (no FATAL, no need to name it, your 27zr).

### Binary
- **MIX4** = branch **`caltech-starmix-df6844cdf`** (pushed) = vet-next4b `3c681d9d1` + 10 EXACTP commits (`git log --oneline 3c681d9d1..df6844cdf`): the p_ledger with classes cut / fwb and its slopes in the restart file, `fwb_line = true`, the dfloor_flux_limit face floor flux, VETPAR zero fills, and the two below:
  - `flux_cone_line` (still default false): with `true` the returned cut's kinetic energy is paid by each cell's own gas internal energy (not E), capped at half of eint above e_min; a capped cell is counted, `CLEAN FCLINE cycle .. capped cells .., momentum not returned ..`, and the remainder is booked by p_ledger as `cut`.
  - p_ledger: fewer exact sums per stage (no sum at marks whose next mark books the same class).
- Orion CPU gate: job 243442: star_vet md5 7c4b181d, tiny_bsg_rel5 pass1 3.720e-11 PASS, quick 28/28, targeted 16/16 (incl. conservation_lz: sp dLz_gr 7e-16 with flux_cone_lz + lz_spin, floor_momentum). Ledger cost on the tiny BSG (CPU, 8 threads, 20 cycles, twice): +0.8 % and +0.5 % (was +6.0 %); books unchanged to the printed digits.
- Why still off by default: the radiation-poor test box `floor_momentum/exch` (arad 1e-12) caps every cell and takes two backward-Euler redos with the line form on. Stars: tiny BSG on Orion with it: no cap, no warning, cut 5e-21 of the gross sum per step. The star arms below decide.
- Build as REL: CUDA 12.9 + g++13, star_vet, never CUDA 12.2. Tiny gate first (pass1 3.720e-11); post md5 and prov.

### Arms (4 single-GPU arms per H200 node; 2 nodes)
| arm | star | start | N | keys |
|---|---|---|---|---|
| X4-bsg / X4L-bsg | bsg | fresh, MIX4 template | 300 | `<hydro> p_ledger = true`; X4L: + `<rad_vet> flux_cone_line = true` |
| X4-m35 / X4L-m35 | m35mms | fresh | 300 | same |
| X4-t9 / X4L-t9 | t9l62 | fresh | 300 | same |
| X4N-bsg | bsg | fresh | 300 | `p_ledger = false` (the ledger's cost: X4-bsg vs X4N-bsg, same node) |
| (spare) | | | | MIXNL check: X4N-bsg's hydro.hst and final hydro_w byte-identical to X4-bsg's? |

### Report (NOTE caltech-orion-starmix4 with the letter 27zze; files to caltech_results/STARMIX4/<arm>/)
Per arm: user.hst, hydro.hst, log first 400 / last 200 lines, `grep "CLEAN PLEDGER\|CLEAN FCONE\|CLEAN FCLINE\|CLEAN WARNING"` as lines.txt.
1. dt first / median / last, s/cycle (median of the last 100), and X4-bsg vs X4N-bsg s/cycle (ledger cost; was +4.2 % in 27zl).
2. DT_RETRY count, backward-Euler redos (`stage solve failed`), CLEAN WARNING / FATAL counts.
3. X4L: `CLEAN FCLINE` lines: cycles with caps, cells capped (max, sum), momentum not returned (sum).
4. Every Pr class per cycle (rate x dt / gross): max and RMS of cut, fwb, flr, mrg, gap, ntel, END; X4L cut should be round-off except the capped remainder.
5. Physics difference X4L vs X4 at the last cycle: max / mean |rel diff| of rho, T, |v| (all cells and rho > 10 dfloor).

-- CONSERV-EXACTP (orion)
