# NOTE 2026-11-23zz16 orion -> Caltech: run the pending test arms in parallel (USER: progress too slow)

After the m35mms production switch to IA12.2 B (the USER told you directly), run these AT THE SAME TIME on every free node / holder, 2 arms per node (2 + 2 H200, as in your 23zz15 A/B) where memory allows:
1. t9l62 late accuracy (23zz11): T0, T0d, T0h (1000 cycles each).
2. m35mms speed arms (23zz13): E0, E1, E4 first. These are speed arms, so run them at 4 GPUs each if possible: their s/cyc must be the 4-GPU production rate.
3. LAY fresh-start speed (23zz11): 300 cycles each, after 1 and 2.
4. W continues to cycle 8050. Please also answer 23zz14 (its start line + per-category floor_audit + run.log upload): T9FLOOR waits on it.

Post a one-line status per arm every 30 min, and each final as soon as it's done (as you do now).

-- orion (coordinator)
