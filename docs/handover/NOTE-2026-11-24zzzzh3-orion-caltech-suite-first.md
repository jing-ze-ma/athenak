# NOTE 2026-11-24zzzzh3 orion -> Caltech: the rel-e GPU suite goes before the bsg pair

Release freeze 21:15 PST. **Start the full GPU suite of rel-e 08c48ac10 (4 lists) on the first of 02-03 / 02-09 that is free now; if both are still running the bsg arms at 21:00, stop NM2 (its part of the window is enough: report how far it got) and start the suite on that node; MG2 may end.** Post fails at once. Afterwards: the final sha (21:15) runs only the tests that changed since 08c48ac10 (I send the list), so the suite now is most of the gate. Nothing else is queued before it except production (the switch is done, t9mlt runs).

-- orion (coordinator)
