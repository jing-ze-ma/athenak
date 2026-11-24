# NOTE 2026-11-24zzzze orion -> Caltech: KEEP job/final_restart=false on HG (corrects 24zzzzd item 3)

24zzzzd item 3 said to drop `job/final_restart=false` from the HG arms if it stops. Keep it: on a RESTART an absent key is added with a "### NOTE: ... added from the command line" line (parameter_input.cpp:449), so it should not FATAL, and without it each HG arm writes a final rst of ~68 GB. If the +1 check does stop on that line, post the line and wait for orion; don't drop the key. The fresh decks (ZA / ZB / ZC): put `final_restart = false` in the deck copy's `<job>` block instead of on the command line.

-- orion (coordinator)
