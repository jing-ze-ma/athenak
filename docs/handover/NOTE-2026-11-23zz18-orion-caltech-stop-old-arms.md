# NOTE 2026-11-23zz18 orion -> Caltech: stop testing old combinations; one new candidate binary is coming

USER: we are running too many tests on old keys while newer updates exist.
- **Stop W now.** Its answer is already clear: the floor fix covers only ~10 % of the floor energy (the kedom cells), and the dominant skip(dfloor/vceil) path is untouched (your 22:05 evidence). T9FLOOR / VCOV are fixing that path.
- Let T0a / T0da / T0h finish their 1000 cycles (~30 min): that answers the diff key's late-state accuracy.
- The m35mms E arms (23zz13): run only E0 vs E4 (all three speed keys), on IA12.2 = your production binary. Skip E1 / E2 / E3 / E0n / E1p.
- LAY fresh-start speed: still useful, run when a node is free.
- **Coming in ~1 h:** IA12.3-rc1 = the cleaned release trunk + all gated updates, with a deck set naming every current best key. After that, all t9l62 / m35mms tests use rc1 only. Don't build anything older.
- m35mms production on IA12.2 B: unchanged.

-- orion (coordinator)
