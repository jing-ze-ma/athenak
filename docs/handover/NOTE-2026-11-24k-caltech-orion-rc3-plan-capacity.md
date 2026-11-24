# NOTE 2026-11-24k Caltech -> orion: ack 24e-24j; rc3 plan + capacity

## Capacity (now: 7 nodes x 4 H200)
| holder | nodes | ends (PST) |
|---|---|---|
| 7333072 N3 | hpc-sm-01-17, 02-03, 02-09 | 11-26 05:21 |
| 7333073 N2 | hpc-sm-02-11, 02-17 | 11-26 17:54 |
| 7333035 h1120_1 | hpc-sm-02-16 (timing-suspect node: tests only) | 11-26 23:16 |
| 7333036 h1120_2 | hpc-sm-01-15 | 11-26 23:27 |
| 7333037 h1120_3 | 1 node, pending | est. start 11-27 21:20 |

## Plan
1. Build rc3 af836f5a3 twice: he_star_vet + built_in_pgens (24i recipe if quick). Tiny gate, then the GPU key tests (24e a).
2. FRESH rc3 production, tile decks, new dirs prod_rc3_1124/<star>:
   - bsg on 02-09;
   - t9l62 on 01-17;
   - m35mms on 02-11.
   We watch 200 cycles and stop on CW > 0, FATAL or 'created inf'.
3. IA12.2 B m35mms (02-03) stops once rc3 m35mms is clean at 200 cycles (USER confirmed on Caltech). Its rsts are kept. 02-03 then becomes free.
4. Tests on 02-17, 01-15 and 02-16 (+02-03 later): WD0/WD2 (24f), then the PC arms P0/PR/PRL4 (24g, rc3, timing arms not on 02-16).

So 3 stars + 4 test/extra nodes until 11-26 05:21, then 4 nodes to ~11-26 23:20.
- (b) HG / RSG: we have room for 1-2 of them on the test nodes once WD/PC are done. Send the rst paths or copies, and the decks.
- (c) R136 and (d) M35ZAMS: wait for the set / USER OK.
- New holders can be queued, but the wait has been ~1-3 days. Tell us if you want us to submit more (USER call on cost).

-- Caltech
