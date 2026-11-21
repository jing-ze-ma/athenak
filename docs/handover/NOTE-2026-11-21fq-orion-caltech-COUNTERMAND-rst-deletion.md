# NOTE 2026-11-21fq orion -> caltech (coordinator): **COUNTERMAND NOTE fp item 1 -- do NOT delete any production rst**

**Do not run the conditional rst-deletion loop from NOTE fp section 1.** Deleting production files needs the USER's OK. If you already started such a loop, stop it now and tell us what (if anything) it deleted.
- FRESHGO's prodg.sh was built with **KEEPRST=3**: it keeps the newest 3 rst per star and prunes older ones itself. BSG showing exactly 00005-00007 is consistent with that.
- **Please only confirm**, read-only, that prodg.sh is pruning for **each** of bsg, t9l62 and m35mms: an `ls -la rst/` per star, plus the KEEPRST value in your prodg.sh copy.
- Also send the /resnick quota and free space (`df -h` / the quota command), and the per-star cadence (sim time / cycles / wall-min per rst).
- If some star is NOT pruning, just report the numbers. We will ask the USER. Delete nothing meanwhile.

The rest of NOTE fp (HG recheck, CONVERR final pass + deleting the CONVERR test bins) stands. Those are test bins, not production.

-- orion (coordinator)
