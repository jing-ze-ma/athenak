# NOTE 2026-11-20dk orion -> caltech (USER): when the beam arms finish, start 8 single-node arms at once

When the TG*/XCP* 100-cycle arms free their nodes, start all of these **together**, one node each (8 arms = 8 nodes):
1. R_LF and R_off (`arms_rsgdie.txt`).
2. LFa_rsg, L0a_rsg, LFb_rsg, L0b_rsg (`arms_lfseed.txt`).
3. G_on_rsg and G_off_rsg (`arms_rush.txt`).
Keep one node back only if a beam/thinf follow-up arm arrives urgently. If so, R_off (short) or the G RSG pair can wait for it.
Keep the 2-node reservation from 11-21 13:00 PST: the 3000-cycle seeds (~2.8 h) and R_LF (~2.6 h) finish well before it.

-- orion (for the USER)
