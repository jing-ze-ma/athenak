# NOTE 2026-11-21fi orion -> caltech (USER via coordinator): KEEP T19 + T19 corona running (watched); queue WGI_t9l62_corona (low priority)

1. **The USER decided to KEEP T19L6.4 and the T19 corona test running**, watched closely. Do not stop them.
   - The PROD_<star> monitor files from NOTE fh (full hydro + user hst every ~30 min) are what we track: KE_int (user col 11) and L_bot/L_in (cols 3/7) every 50 cycles.
   - Please make sure T19, the T19 corona, T9 and m35mms are always in each refresh (BSG too).
   - Also please give the **cycle numbers of the 2 DT_RETRY events** each run had in the first hour (cycle 0?).
2. **Low-priority test, free node only, never a production node:** WGI_t9l62_corona, the isothermal 2.1 MK corona.
   - The files are in 1104bbfcd: decks/t9l62_coronaiso_wngl_smoke.athinput (md5 8a395d6b) and t9l62/corona_test/ic_t9l62_corona_iso.txt.gz (md5 846d2dd1). Gunzip it into your FRESH copy.
   - `RT=1 armg.sh WGI_t9l62_corona $B_v203 decks/t9l62_coronaiso_wngl_smoke.athinput 100`
   - Order: right after the ST/STD arms (NOTE fg list).
   - Its line should include the clean_init steady line, dt first/last, and L_top/L_in at cycles 0 / 50 / 100.

-- orion (coordinator, for the USER)
