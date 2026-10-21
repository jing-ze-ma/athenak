# NOTE 2026-10-21 raven -> deltaai: RSG E disk plan (USER: do all three)

1. Move the run's output to /projects/bivj (1 TB free): new links write there (run dir / bin / rst on /projects; keep the old
   /work/hdd files until the move is verified, then delete them).
2. Restart cadence: `<output5>` (file_type = rst) `dt = 432000.0` (5 d) instead of 43200.0 (0.5 d). Output-cadence only, no
   physics key; this is the ONLY permitted deviation from template 4c442452 — report the new deck md5 and the diff (3 @FILES@
   lines + this one). Apply at the next link boundary (restart from the latest rst; no fresh start).
3. Prune: keep the newest 2 rst (not 4).
Report the new disk forecast (GB per wall hour, days to quota) in your next NOTE.
