# NOTE 2026-11-16 Orion -> viper: CORRECTION to the stag2 progress NOTE 6d7b20a97

A full benchmark of the "recommended" staggered key set against the cell layout (d6b6ed4e2) found two regressions that the round tests
had missed. The report is at Orion /orion/ptmp/jinma/clean_1115/stag2/bench/BENCH.md.

1. **rw accuracy.** With clean_stag_fsolve, staggered is FIRST order at tau 0.1 and 10; orders 0.67-0.94.
   N256 L1(E): P 1/3 tau 0.1 2.9e-2 vs cell 3.1e-4; P 3.3 tau 10 2.8e-3 vs 6.6e-5.
   - split alone, or split + consistency: order 2.04 and better than cell (P 1/3 tau 10: 9.7e-4 / 2.4e-4 at N64/128).
   - Cause: the F solve's first-order LF jumps.
2. **Well-balancing.** On rhse_sp2, max v/c_s is 2.2e-5 vs cell 1.7e-7. split alone gives 6.4e-9; consistency raises it to 7.3e-6,
   and fsolve/fkin raise it further.

## The rest of 6d7b20a97 stands
- staggered is 1.2-4.6x cheaper per cycle on Cartesian and 1.25-1.6x on sp in this benchmark, a fairer timing than the round runs;
- flz q/Q and sp3 decomposition results stand;
- with the chain on, every staggered final state has E >= floor, while cell ends negative on cb/shd/shdg.

## Next on Orion
- A second-order F solve: lagged PLM jumps, as in the cell layout.
- Well-balanced F dissipation.
- A WB check for the consistency term on staggered.
The "recommended set" in 6d7b20a97 is withdrawn until these pass.
