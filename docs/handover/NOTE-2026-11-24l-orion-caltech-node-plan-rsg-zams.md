# NOTE 2026-11-24l orion -> Caltech: USER node plan: 4 stars incl. M35ZAMS on 4 nodes, RSG on 3 nodes from its latest rst (rc3)

**USER:** 4 nodes for 4 fresh stars on rc3, including M35ZAMS (USER OK now); 3 nodes for RSG continued from its latest restart on rc3.

## 4 stars, fresh, rc3 (one node each)
- bsg on 02-09, t9l62 on 01-17, m35mms on 02-11 (as you planned);
- **M35ZAMS on 01-15**: deck `docs/handover/ia123rc3_decks/IA123_m35zams_wallsym.athinput`. It is the paper's 0.064 rad narrow wedge, S10 (vet_gd_ls 10), wall_interp + wall_sym + theta_seam_sym, 16 blocks of 384x32x32.
  - If it FATALs at start on the seam / set check, post the line; don't edit the deck.
- Same 200-cycle watch as the others (stop on CW > 0, FATAL, "created inf"). M35ZAMS also: post the max |F_phi|/|F_r| in the thin layer at cycle 200 (the earlier spurious F_phi should be gone).

## RSG on 3 nodes (02-17, 02-16, 02-03 once IA12.2 m35mms stops)
- The latest RSG rst we can reach is `rsg.01707.rst` (48.5 GB, 11-20 04:18) on viper: `/viper/ptmp2/jinma/prod_1107/rsg/rst/rsg.01707.rst`. A newer one may sit on Raven, which is down for maintenance.
- **Question:** can you pull a 48 GB file from viper (or Orion /orion/ptmp) yourselves? Or does the USER need to transfer it (Globus / scp)? Say which.
- Deck + keys: we'll send an rc3 RSG restart-args file that names every key (the rst stores old keys), once we know the transfer route.
- Until the rst arrives, use those 3 nodes for the tests: the GPU key tests (24e a), WD0 / WD2 (24f, ~25 cycles each), then the PC arms P0 / PR / PRL4 (24g, on rc3; timing arms not on 02-16).

-- orion (coordinator)
