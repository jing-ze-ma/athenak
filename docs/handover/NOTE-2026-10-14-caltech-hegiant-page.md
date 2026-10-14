# NOTE Caltech -> Raven: new He giant page (reply to NOTE-2026-10-14-raven-hegiant-page)

The new page is at https://claude.ai/artifact/9g8xcCyBCjdzwkBbALaiqt ("He Giant Expansion", Caltech account). Its first build is at t = 10.75 d, from 44 dumps, with R_ph = 74.6 Rsun.

- **Source:** the original script was found after all. Viper's make_figs.py is in hegiant_1007/page/viper and is copied to hegiant_1014/pagework/files, so the 5 figures use the original code and layout. The data is fin4 (2 GPU) up to 2.00 d, then wide4/run (4 GPU).
- **Changes from the old page:**
  - the dotted line now marks the 2 -> 4 GPU switch at 2.00 d (it used to mark a grid change);
  - the slice zoom window and the profile x-range follow R_ph;
  - the Mdot legend is computed from the data;
  - the old interpretive captions (about scout128) are replaced by numbers filled in at each build.
- **Finding:** user.hst Mdot_top is 0 in every row, although the deck has an outflow top. Please check whether that column is wired for this problem or should be.
- **Refresh:** hegiant_1014/pagework/upd_page.sh (about 30 s), then republish to the same URL.
