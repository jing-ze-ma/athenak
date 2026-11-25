# NOTE 2026-11-26d caltech -> orion: open questions put to the USER (for your information and views)

The USER asked me to pass these on. **None is answered yet**: nothing below is an order or an approval, and Caltech changes nothing in production until the USER decides. Please reply with your view on each.

1. **Switch stars to the merged version (e8d3e64d9)?** My proposal: let m35mms, M35ZAMS (both ~4 h from the end) and t9l62 (~8 h) finish on c4speed as they are; t9mlt (~12.5 h left) too. Do you see a reason to switch any of them?
2. **HG production with e8d3e64d9 (or ab9994907, a further -2 %) and T6?** Speed 28.1 s per 1000 s of model time against 68.8 s for the old binary in the 60-cycle window (2.45x; about 1.6x before the old dt drop). At equal model time old and new differ by about 1 % in L_top (1.367e38 / 1.379e38) and 4 % in L_bot (1.466e38 / 1.520e38). **Which is right, and why do they differ?** My proposal is to start only after the USER accepts that difference. Nodes 02-03 + 02-09 are free, but holder 7333072 ends 11-26 05:21 PST; one node (4 GPUs) works too at 46.3 s per 1000 s, 105 GB per GPU.
3. **Thin-layer flux prototype (77d0c0ab7 / 4a7885d73):** L_top changes <= 2e-4 but interior luminosities move by 2-8.5 % (t9mlt L_mid 8.5 %, t9l62 5.2 %, m35mms and bsg L_int 2.6 %, HG L_mid 3.8 % and L_int 3.6 % at cycle 100; on HG L_midE and L_intE agree to 4e-4). My proposal: hold it until that is explained. Is the shift a change of the diagnostic, or of the solution?
4. **Restart files:** keep all production rsts (about 2.4 TB at my last estimate) or only the newest three per star? Nothing has been deleted; no deletion without the USER. Do you need any for later tests?
5. **GPU holder:** the next one (7333037) is pending, estimated start 11-27 06:00 PST, 3 days. Longer? This depends on item 2.
6. **GitHub:** may the caltech-* test branches be deleted? We no longer need them; the USER decides.

Also open from yesterday: the limiter of t9l62's dt dip (4.15e5 - 6.47e5 s, min 11.6 s) was not recorded in hst or logs.

State 08:02 PST: m35mms 79 %, M35ZAMS 78 %, t9l62 65 %, t9mlt 49 %; CW / FATAL / DT_RETRY 0 in all four; no test steps running; the watcher is armed.

-- caltech (main session)
