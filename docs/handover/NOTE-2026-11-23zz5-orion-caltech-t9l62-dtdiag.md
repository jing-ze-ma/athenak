# NOTE 2026-11-23zz5 orion -> Caltech: t9l62 dt-limiter diagnostic arms (USER goal: t9l62 in ≤ 24 h on one node)

Run these short arms on the next free nodes, AFTER arm W (23zz4) is started. Priority: W > T0/T0d > T1/T1d/T2d > LAY/LAY0.


Posted by COORD. Goal: find what limits t9l62's dt in each phase, and whether
`hydro/rad_signal_speed_diff` (USER-adopted) or a 1-block/GPU layout shortens the run.

## Binary
- `caltech-m35fast-13679216a` (the m35mms D1/D2 binary: IA12.1 + rad_signal_speed_diff + max_stage2; DTREP prints the diff weight).
- Every arm: `IA121_t9l62_restart.args` (name every IA12.1 key), env `ATHENA_DT_REPORT=5`, rst/bin dt 1e30, `job/final_restart=false`, 4 ranks.
- Stop rule as 23zk (dt < 0.01, NaN/inf, FATAL -> scancel the step).

## Arms (each ~20-40 min on 4 H200; run 2 at a time if a node is free)
| arm | rst | extra keys | nlim |
|---|---|---|---|
| T0 | prod_fresh_1121b/t9l62/rst/t9l62.00031.rst (J4, t ~1.34e6; else 00037 in rst_p1_seg) | none | rst cycle + 30 |
| T0d | same | `hydro/rad_signal_speed_diff=true` | rst cycle + 30 |
| T1 | prod_ia11_1123/t9l62/rst_keep/t9l62.00010.rst (t ~4.3e5) | none | + 30 |
| T1d | same | `hydro/rad_signal_speed_diff=true` | + 30 |
| T2d | prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst | `hydro/rad_signal_speed_diff=true` | + 30 |
| LAY | t9l62.00010.rst | `meshblock/nx2=256` `meshblock/nx3=128` (4 blocks, 1 per GPU; tile 128 unchanged) | + 100 |
| LAY0 | t9l62.00010.rst | none (2 blocks/GPU, reference) | + 100 |

- T0/T0d: the J4 rst comes from an older binary. If the restart FATALs on a missing key/array, report the line; the
  cycle-0 DTREP line is still the answer (dt is computed from the read state), so even 1 cycle is useful.
- LAY: if vet_gd / tile FATALs on nx2=256, report the FATAL and skip it.
- If the W / SB arm is past cycle ~6600, add 100 cycles of `ATHENA_DT_REPORT=10` there too (no restart needed if you
  can run a continuation from its final state; otherwise skip).

## Report back
- `grep "DTREP\|DTREG" log` for every arm (first 30 cycles), + run.log tail (s/cycle, elapsed per 10 cycles).
- LAY vs LAY0: s/cycle over cycles 20-100, and the hst max relative difference (expect ~1e-12, tile-invariant).
- From your existing IA11 prod log (prod_ia11_1123/t9l62): the write times of one rst and one bin of each kind
  (hydro_w, m1, m1_vet, m1_face) — grep `written at` / output timing lines.

## How Orion reads it
- x = tauT_h * ceff / 2.998e10 at the binding cell: x >> 1 trapped (physical, key changes nothing);
  x < 1 the key lifts dt toward cfl*dt_gasonly.
- T0 vs T0d is the late-phase (70 % of the cycles to tlim) dt gain; T1/T2 the rise / peak.

-- orion (coordinator; drafted by HALFDAY-T9)
