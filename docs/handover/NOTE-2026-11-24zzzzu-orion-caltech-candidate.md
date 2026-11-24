# NOTE 2026-11-24zzzzu orion -> Caltech: profile with the CANDIDATE settings (USER); replaces the settings of 24zzzzs / 24zzzzt

CANDIDATE = binary **09d1e4df4 (md5 b8b0d79c)** with every key named: `rad_vet/time2_gtime_slope=true clean_stag_thinf=true clean_stag_lat_sums=true clean_stag_tie_gid=true clean_stag_opac_drag=true vet_gd_allow_clamp=false`, `hydro/vceil=3.0e8`, `rad_vet/clean_pc_sweep=rb`, mg 4 for m35mms / 2 for the others, cfl 0.3, rays as before. (Not adopted yet: for timing only.)
1. **t9l62 (24zzzzs PL / PO / PQ):** same states, CANDIDATE settings, clean_verbose=2; PL needs +40 cycles first. Only ONE setting: no off arms.
2. **M35ZAMS (24zzzzt):** mg 2 / 3 / 4 with CANDIDATE settings, 100 cycles.
3. **Each star, candidate timing (one arm per star, 100 cycles from the latest production rst AND a fresh 100-cycle start from the deck values):** s/cyc over the last 50, dt, its per stage, passes, CW / BE / retry / clamped reads / floor_pay. Fresh starts: the decks in wt_note ia123rc4_decks are OLD: use the keys above on the command line only where the deck names them already (a key absent from a fresh deck is FATAL on the command line: add it to the deck copy), and say what you changed.
4. Projection (24zzzzq item 4), with hours at the CANDIDATE s/cyc scaled per state (label measured vs scaled).
Order: 3 (bsg, M35ZAMS, m35mms, t9l62), then 1, 2. Idle nodes only; HG / S arms keep priority. If anything FATALs, post the first error line and go on.

-- orion (coordinator)
