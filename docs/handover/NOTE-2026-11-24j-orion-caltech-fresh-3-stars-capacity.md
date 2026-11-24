# NOTE 2026-11-24j orion -> Caltech: start bsg, m35mms and t9l62 FRESH on rc3 (USER); how many more stars can you run?

**USER:** start all three stars fresh on the latest binary, i.e. rc3 = caltech-ia123rc3-af836f5a3 (NOTE 24h), decks in `docs/handover/ia123rc3_decks/{tile,notile}/IA123_{bsg,m35mms,t9l62}.athinput`.

## Order
1. Build rc3 (incremental; NOTE 24i recipe if convenient). Tiny gate with the seam keys: pass 1 = 3.686e-11.
2. The GPU key tests (24e (a), on the built_in_pgens build): post PASS / FAIL.
3. Start **bsg, m35mms, t9l62 fresh** (tile decks = production form), each on its own 4-H200 node, with the decks' tlim and output / rst cadence as in the IA12.2 production. Watch the first 200 cycles; stop and post if CW > 0, FATAL, or a "created inf" line appears.
4. The current m35mms IA12.2 B production keeps running until the fresh rc3 m35mms is past its first 200 cycles cleanly. Then stop it, keeping its rst as a reference.
5. t9l62: the outburst (t ~4.8e5) is days away. The fix from WD0 / WD2 (24f) will come in by a restart long before then.

## Capacity question
How many 4-H200 nodes can you hold at once now, and for how long (holder end times)? Besides the three stars we'd want, in order:
- (a) a test node for the WD0 / WD2 and PC arms (24f, 24g);
- (b) HG and RSG from their old-production rsts, if you have them or can take a copy (~30 GB / ~20 GB);
- (c) R136, once its mirror-symmetric rotating set lands (in progress);
- (d) M35ZAMS (held: needs a USER OK).
Reply with your node count and what fits.

-- orion (coordinator)
