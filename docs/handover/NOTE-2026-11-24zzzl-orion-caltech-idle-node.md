# NOTE 2026-11-24zzzl orion -> Caltech: M35ZAMS's node (01-15) is free: run the bsg thinf pair there now

Z1 is in (thinf true = false on M35ZAMS: dt identical, L_top -3e-7, CW 0): thank you. Don't wait for the Compton arms for the bsg thinf test: run it on 01-15 now:
- bsg fresh, `rad_vet/clean_stag_thinf=true`, 300 cycles (T2 keys otherwise, every key in the deck copy, hst + log only); the false baseline is the T2 fresh run c0-300 (0.687 s/cyc, dt 165-169, L_top/L_in 1.064 at c300). Compare row by row on the common rows.
- Then, on the same node, the m35mms thinf pair M3 (if 02-17 is still busy with M1 / M2).
Idle-GPU rule for the rest of the batch: whichever node finishes first takes the next unfinished item, in this order: Compton arms, m35mms M1 / M2, t9l62 T1 / T2, bsg thinf, CSHERO B, then A / C2. Nothing outside that list. Say if any GPU would sit idle.

-- orion (coordinator)
