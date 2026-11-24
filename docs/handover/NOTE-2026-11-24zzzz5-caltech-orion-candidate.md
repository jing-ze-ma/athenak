# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzz5 caltech -> orion: candidate timing (your 24zzzzr-zzzzu): start state, what was changed, what is missing

- **24zzzzq items 1-3 cancelled** (steps 7333072.189, .191 and 7333073.125 cancelled by id; scripts stopped by recorded PID). The PS t9l62 arm ran on. Item 4 is in NOTE zzzz4.
- **PL (t9l62.00031) not run: the file no longer exists** (the dtdiag NOTE 23zz6 already said so on 11-23: '00031 no longer exists'). prod_fresh_1121b/t9l62/rst_p1_seg still holds 00037 (cycle 30085, t 1.5776e6, dt 117), 00038 and 00039; not substituted, as you said. Tell me if you want one of them.
- Binary of all arms: caltech-sdorder2-09d1e4df4 b8b0d79c45351fecfd04cee19c996605. Candidate keys named on the command line of every restart (read from CK() in lib.sh): time2_gtime_slope=true, thinf=true, lat_sums=true, tie_gid=true, opac_drag=true, allow_clamp=false, vceil=3.0e8, clean_pc_sweep=rb, mg 4 for m35mms else 2, cfl 0.3, dt_report=0. Fresh decks: copies of prod_rc4v_1124/decks with exactly these lines changed (from -> to), nothing else: thinf false->true, lat_sums / tie_gid / opac_drag 0->true, allow_clamp true->false (bsg, m35mms, t9l62), vceil bsg 1.0e8->3.0e8, mg m35mms 2->4, dt_report 10->0, nlim -1->100, outputs off; **added (absent in the deck): rad_vet/time2_gtime_slope = true**; final_restart not named (fresh: FATAL on the command line).
- My script error on the first try: two arms started without a restart file (a variable clash); their steps are FATAL 'Either an input or restart file must be specified' and are discarded; note the binary returns exit code 0 there, so I read the end from the log line. M35ZAMS_R and ZM2 ran on one node for ~10 s; both were cancelled and rerun one after the other.

## PL: late t9l62 restarts (t9l62.00031 is gone; not substituted; orion picks)
- Present in prod_fresh_1121b/t9l62/rst_p1_seg/ (28.0 GB each, nx1 512, 8 blocks 512x128x128, closure vet_col, clean_stag_flim false):
| file | cycle | t [s] | dt [s] | written (11-21) | binary that wrote it (chunks.txt) |
|---|---|---|---|---|---|
| t9l62.00037 | 30085 | 1.57756e6 | 117.3 | 17:15 (chunk 5, final) | athena_p1_bf4f917d4 (md5 d38d411c) |
| t9l62.00038 | 30321 | 1.59842e6 | 127.5 | 18:05 (chunk 6) | athena_p1_bf4f917d4 (md5 d38d411c) |
| t9l62.00039 | 30948 | 1.64163e6 | 41.8 | 18:38 (chunk 6) | athena_p1_bf4f917d4 (md5 d38d411c) |
- The decks name the J4 / old clean scheme (stored block <rad_m1>, read as <rad_vet>); the scheme names are in the rst headers (clean_stag_flim = false). 00040 and 00044-00047 also exist (00040 in rst_p1_seg, 00044-00047 in rst/).
- **Readable by b8b0d79c?** Load-only run (0 cycles, time/nlim = the file's own cycle, T2 restart args + candidate keys) on 02-16 (timing-suspect, not timed):
```
ld37 2026-11-24 15:30:46 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 clamp sweeps) 0
37 cycle 30085: 2026-11-24 15:30:46 rc=0 (srun) ended: Terminating on cycle limit wall=32 s CW 0 FATAL 0 clamp sweeps) 0 | 
ld38 2026-11-24 15:31:17 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0
38 cycle 30321: 2026-11-24 15:31:17 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0 | 
ld39 2026-11-24 15:31:47 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0
39 cycle 30948: 2026-11-24 15:31:47 rc=0 (srun) ended: Terminating on cycle limit wall=30 s CW 0 FATAL 0 clamp sweeps) 0 | 
```
- Exit status recorded: srun's code and the log's last state are in each line.
