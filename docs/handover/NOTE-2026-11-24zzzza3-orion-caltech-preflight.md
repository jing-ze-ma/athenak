# NOTE 2026-11-24zzzza3 orion -> Caltech: PREFLIGHT of 3 fresh stars on the final set (no production start yet)

HG G passes (thanks). Decks: `docs/handover/ia123rc4_decks/candidate_1124/` (MD5SUMS there; tile decks + M35ZAMS wallsym; they name time2_gtime_slope = true, so they need binary **caltech-sdorder3-95813cf85, md5 30d01765**). The x1 grids / tables are the ones the tile set already uses; say if a path differs on your side.
PREFLIGHT (idle nodes, one node each, fresh t = 0, 100 cycles, outputs off except hst + log, `job/final_restart = false` in the deck copy): **bsg, M35ZAMS, m35mms** (NOT t9l62: its steady init warns with thinf = true, T9LATE works on it). Check and post: every candidate key as printed in the log / resolved header, CLEAN WARNING = 0 (including clean_init), FATAL / inf / BE / fallback / retry / clamped reads 0 (first BE step of a fresh run excepted), s/cyc, dt, GPU memory, L_top/L_in at cycle 100.
A failure of any of these: stop that star and post. Do NOT start production; the USER gives the go.
Then S3G / S5G as ordered (24zzzza2) have priority over the preflight for nodes.

-- orion (coordinator)
