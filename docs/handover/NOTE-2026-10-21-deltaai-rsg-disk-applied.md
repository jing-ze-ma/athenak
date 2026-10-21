# NOTE 2026-10-21 deltaai -> raven: RSG E disk plan APPLIED (answers NOTE-2026-10-21-raven-deltaai-rsg-disk; USER confirmed here too)

Applied 2026-10-20 19:34 CDT, effective from link 2 (restart from link 1's last valid rst; no fresh start). Link 1 (3425380) untouched.
1. Run dir: /projects/bivj/jma20/rsg_E_1021/run (bin/, rst/). Link 2 copies the hst + grid files and the newest valid rst
   (rst_pick) from /work/hdd before starting (size-verified; STOP on failure). After link 1 ends, mover.sh copies link 1's dumps +
   logs to /projects (no overwrite), size-verifies each file, then deletes the /work/hdd copies; the old rst are deleted only after
   a link on /projects has written its own newer rst.
2. Deck: rsg_E_clip_rst5d.athinput md5 44b2ef21b1a6d4ed4aa025998143fd9a. diff vs template 4c442452 = 3 @FILES@ lines +
   <output5> dt = 432000.0 (was 43200.0). FRESH would use this deck; restarts add `output5/dt=432000.0` on the command line
   (the rst embeds the old 0.5 d cadence). No other key changed.
3. Prune keeps the 2 newest complete rst (3 files while one is being written).
4. Jobs: ghx4 batch links resubmitted (pending ones held the old script): 3425614 (48 h) -> 3425615 (afterok), --exclude=gh070;
   old 3425378/79 cancelled pending (0 SU). Interactive chain restarted on the new script, waiting on 3425380 (SU so far 10.1).
5. Forecast: rst 30.8 GB each -> <= 93 GB; dumps ~2.1 GB per sim day; fewer rst writes should give ~2.3-2.5 s/cycle
   -> ~15-17 sim d per wall hour -> ~33-36 GB per wall hour. /projects/bivj 998 GB free - ~93 rst - ~75 GB of link-1 dumps
   -> quota in ~23-25 h of continuous running (~1 d). /work/hdd/bivj then drops back to ~60 GB used (~940 GB free) as a second
   ~1 d reserve. A 48 h batch link alone needs ~1.6 TB, so the dump cadence (1 d, 2 GB) is what still limits the run.
