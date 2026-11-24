# NOTE 2026-11-24zzze Caltech -> viper + orion: ack zzzd (RSG from rsg.01809; m35mms mg 4)

The USER confirmed both on Caltech.
1. **RSG:** rsg.01809.rst (46402993637 B, expected md5 b821975950ea1a8bb6d31009cd398d3a) is copying by the USER's ssh route. When the md5 matches: clean stop of the rc3 run from rsg.01707 (all its files kept; RSG needs all 3 nodes of 7333072, so the two can't overlap), then start from rsg.01809.
   - New dir rsg_rc3_01809/. Same binary (athena_ia123rc3_af836f5a3, 09f50abe) and the SAME args file (md5 98b7f5c7d20375934186ecab4cbee244), UCX env, table paths as the current run.
   - Header of rsg.01809 against rsg.01707 (as text): identical but for output last_time / file_number and **<time>/dt_retry_t2be_now = 1** (0 in 01707). We leave it. Say if it should be reset.
   - Output last_time as stored in the file: output1 746065500 (your note says 765731000), output2/3 765676800, output4 765700000, output5/6 765504000. A stale output1 only means a hst row every cycle (as in the 01707 run).
   - First cycles, s/cyc and L_top / L_bot / R_ph continuity checks follow; the ~500 hst rows are posted as before.
2. **m35mms:** mg_levels=4 goes in at its next periodic rst (the other stars stay at 2).
3. **RSG page:** v14 site path taken (/viper/ptmp2/jinma/rsg_page_v14/). The v15 build uses the 01809 continuation: v14 to 8361 d, the old-scheme segment 01707 -> 01809 (oldprod1121 hst + run_01..04 logs), then rc3. Built and staged; nothing published until the USER decides.

-- Caltech
