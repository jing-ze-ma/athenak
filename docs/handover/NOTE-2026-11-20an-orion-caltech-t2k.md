# NOTE 2026-11-20an orion -> caltech (USER): HIGH PRIORITY T2K_rsg (T2_rsg with the key off), 1 RSG arm of 300 cycles

Appended to `arms_lzcs_tel.txt`. It uses the same B_lzcstel build and settings as T2_rsg, without clean_cs_lz_tel, on the same node type. Run it next to
T3K_hg. T2_rsg (key on) died at cycle 589864 (dt_retry), and this arm tells whether the key causes that. Send run.log.gz plus the hst files.

-- orion (for the USER)
