# NOTE 2026-11-24zzzza5 orion -> Caltech: START bsg, m35mms, M35ZAMS (USER)

USER: start the three stars now that their preflight passes; t9l62 comes later (its steady init warns with thinf = true).
- **bsg and m35mms: start fresh now**, one node (4 H200) each. **M35ZAMS: start as soon as the Z3 preflight (24zzzza4) is clean.** If any fails its preflight, don't start that one; post.
- Binary **caltech-sdorder3-95813cf85 (md5 30d01765)**; decks from `ia123rc4_decks/candidate_1124` (MD5SUMS), unchanged except paths and job names; new run dirs (e.g. prod_c1_1124/<star>); do not touch the old production or T2 dirs, nor other people's runs. Final rst ON, rst cadence as in the deck, hst / outputs as the deck; `dt_report` off. No command-line keys beyond paths.
- First check at cycle 100, then every 500 cycles to cycle 3000 and then with your usual monitor: CW 0, FATAL / inf / DT_RETRY / fallback 0, clamped reads 0, s/cyc, dt, L_top / L_bot, floor_pay. A CLEAN WARNING, clamped read, FATAL or retry: post at once; stop only on FATAL / inf.
- Holder end times: state each node's allocation end against the projected finish (bsg ~6-7 h, M35ZAMS ~20-22 h, m35mms ~22-24 h) and chain restarts (from the newest rst) if a holder ends first. 
- The remaining tests (S5G, init arms, CSHERO, ...) keep their priority only on nodes the three don't use.
Post the start lines (node, dir, binary md5, deck md5, first cycle).

-- orion (coordinator)
