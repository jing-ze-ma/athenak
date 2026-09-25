---
name: index-general-eos
description: Memory links for the general/tabulated equation of state project.
metadata:
  type: reference
---

## General EOS

- [EOS electron regression test](eos-electron-regression-test.md) — DONE: mhd_eos_electrons.py, plus how to dry-run a test without letting run_tests.py delete build/
- [EOS picked by run name](eos-selection-by-run-name-bug.md) — analysis scripts guessed the EOS from the directory name and silently used the wrong one; fixed, plus the audit that found the scope
- [EOS x_e inert in the resistivity](eos-xe-resistivity-capped.md) — RESOLVED: not the EOS at all; two out-of-bounds reads of eta_b, one of which made ohmic_resistivity=constant a silent no-op
- [eostest input inventory](eostest-input-inventory.md) — what the ideal/general comparison inputs cover and how to run them; .bin output is single precision
- [General EOS: MHD+polar difference](general-eos-mhd-polar-bug.md) — resolved, not a bug: the polar boundary forces HLLD→HLLE at the pole, and HLLE is the solver whose Roe average has no general-EOS analogue
- [General EOS project](general-eos-project.md) — non-ideal EOS for Newtonian hydro/MHD; branch `general-eos`; ALL stages including the analytic EOS are done and verified. START HERE.
- [General EOS: Stage 3 loose ends](general-eos-stage3-loose-ends.md) — the last three TODO(stage3)s closed; why the WB background hands back a TEMPERATURE and not a pressure
- [General EOS: Stage 3 perf cleanup](general-eos-stage3-perf-cleanup.md) — one T solve per cell instead of ~11-30; fully verified and committed as 7d1fed3f on 2026-08-12
- [General EOS: Stage 3 tabulated EOS](general-eos-stage3-table.md) — the analytic EOS itself: H2 + Saha + radiation via a (log rho, log T) table; done and verified, commit 03febafa
- [General EOS Stage 4: (rho,e) table](general-eos-stage4-rho-e-table.md) — SCOPED, NOT STARTED; re-baselined 2026-08-17 to ~1.5x (3.07x -> ~2.0x), proxy ceiling measured twice; task 1 is the (rho,e) domain shape
- [General EOS table cost](general-eos-table-cost.md) — RE-MEASURED 2026-08-17: 3.07x per cycle / 2.32x per simulated second, NOT 4.4x; profile (49% libm) still stands. Restrict any such comparison to the hydro-limited window
- [General EOS: table linear-wave tests](general-eos-table-linear-wave.md) — how general_eos=table is regression tested; the tabulated EOS cannot be called from host code
