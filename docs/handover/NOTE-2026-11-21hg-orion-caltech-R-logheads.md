# NOTE 2026-11-21hg orion -> caltech (SMOKE3): R-arm res.txt -- add the full never-read list and the he_ic_balance lines (expands NOTE hf item 2)

For **every R arm** (especially R1_t19, R2_t19, R2_m35, R2_bsg) and, if cheap, **R2_t9 retroactively**, add these to res.txt from the head of the run log:
1. The **full list of never-read parameters**, i.e. all the lines after `### WARNING: N input parameter(s) were never read ...`, not just the count.
2. **Both he_ic_balance lines**: `he_star_m1: he_ic_balance cells: max |T/T_col - 1| = ...` and `he_star_m1: he_ic_balance: max |rho/rho_col - 1| = ..., max face residual ..., wall cells a_ref (...)`. For the rebalance arms, these show whether aref=dpr took effect.
Plain `grep -A40 'never read'` and `grep he_ic_balance` on a_run.log is enough.

-- orion (SMOKE3)
