# NOTE 2026-11-23zv orion -> Caltech: m35mms arms D1 / D2 / D1h on N2 node A (USER 24 h goal)

The new branch is `caltech-m35fast-13679216a` on GitHub = your IA12.1 fa59bffd6 + 2 commits, both key-gated:
- `hydro/rad_signal_speed_diff` (USER-adopted): the radiation term of the hydro signal speed counts only where radiation is trapped on the grid scale. Wave tests: stable and accurate.
- `rad_m1/clean_gas_secant=max_stage2`.
Gates: key off bitwise vs fa59bffd6; HG / RSG with the key on bitwise vs IA12.1.

## Build
CUDA 12.9 + g++ 13, then the tiny BSG gate (expect pass 1 = 3.686e-11).

## Arms (node A, 4 H200, in this order), from `prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst`, exactly as D0
`env ATHENA_DT_REPORT=25 <production env> $BIN -r .../m35mms.00010.rst -d m35fast_1123/<arm> $(cat <arm>_args.txt)`
- **D1:** key on, secant stage2, cfl 0.3, 2000 cycles. This is the key's effect vs D0.
- **D2:** key on + clean_gas_secant=max_stage2, cfl 0.3, 2000 cycles. Both levers.
- **D1h:** key on, cfl 0.15, 1000 cycles. The dt-halving accuracy yardstick.
- Args: `docs/handover/caltech_m35fast/{D1,D2,D1h}_args.txt` (md5 in MD5SUMS).
- Same returns as D0: user + hydro hst; `grep -E "DTREP|DTREG|CLEAN WARNING|CLEAN FALLBACK|FATAL|^CLEAN cycle|elapsed="` of run.log; s/cycle; GPU memory.

m35mms production, the t9l62 IA12.1 run and the node B arms: unchanged.

-- orion (coordinator)
