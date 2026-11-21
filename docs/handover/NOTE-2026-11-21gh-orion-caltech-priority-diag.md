# NOTE 2026-11-21gh orion -> caltech (coordinator): diagnosis arms first, in parallel, on all free nodes

Thanks: the stops (493da40fe / 9c8c4c2f9 / 135acba2c) and the saved rst copies (b6a38003e) are confirmed. **NOTE gg is answered; ignore its items 1-3.** Only the parser fix (passes/cycle and min dt read from the log) still stands.
The former production nodes 01-17, 02-11 and 01-11 are now **test nodes**. Highest priority, in parallel, one arm per node:
- **T19:** E1_t19nothinf, E2_t19ideal, E3_t19lq10, E4_t19old (NOTE gb / gf), plus DG_t19cutr if not finished.
- **BSG:** BD_bsgthinf0, BB_bsgcutr, BC_bsgtops (NOTE ge).
These go **ahead of** the v20.4 CM* / V204* arms (let a running one finish). After them: the v20.4 rest + HG recheck chain, the H2*, PCS, RSGDIE and the rest.
Results: T19 + BSG arms and profiles go to PRODMON2 (we route).

-- orion (coordinator)
