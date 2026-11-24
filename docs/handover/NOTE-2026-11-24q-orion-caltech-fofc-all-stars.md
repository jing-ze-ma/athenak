# NOTE 2026-11-24q orion -> Caltech: USER: FOFC ON for ALL stars (uniform settings); bsg restart

USER: turn FOFC on for every star, not only t9l62. All the decks meet its needs (nghost 3, eos general, plm). Add to every star's command line (the rsts store fofc = false, so name both):
`hydro/fofc=true hydro/fofc_retest=true`

1. **bsg** (stopped at c777 on "write-back closure inf"): restart from its newest rst BEFORE c777 with FOFC on, as production. Watch it past c777 and post CW / FATAL / "created inf" / the end-of-run "energy created x vol" line.
   - Our BSGINF agent is still finding the cause, and will send separate short test arms (one key toggled each) to run on a test node.
   - If bsg with FOFC fails at c777 again, stop it and post the log.
2. **t9l62**: FOFC on at the next rst (24p).
3. **m35mms, M35ZAMS**: FOFC on at their next rst. Post s/cyc before / after (FOFC's cost) for each.
4. **RSG**: add both keys to IA123_rsg_restart.args when it starts.
5. Ignore "FOFC first-order cell-stages" (broken counter, reads 0). The fix is coming.

-- orion (coordinator)
