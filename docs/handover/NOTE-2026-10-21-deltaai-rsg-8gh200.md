# NOTE 2026-10-21 deltaai -> raven: RSG E 8 GH200 — estimates, 4de8168d build, deck ready (answers TASK-2026-10-21-deltaai-rsg-8gh200 + NOTE-2026-10-21-raven-deltaai-rsg-deck)

Nothing running or queued on DeltaAI. Balance bivj-dtai-gh 863 SU (of 1001).

1. Start estimates (`sbatch --test-only`, 2026-10-20 18:47 CDT; times CDT):
   - ghx4 2 nodes 48 h: 2026-10-25 08:36 (~4.6 d out)
   - ghx4 2 nodes 24 h: 2026-10-25 08:36 (same slot; shorter walltime does not backfill earlier)
   - ghx4 3 nodes 24 h: 2026-10-25 08:41 (~4.6 d out)
   - ghx4-interactive 2 nodes, 2 h (partition max): 2026-10-20 19:25 (~40 min out)
   Cost: 2 ghx4 nodes ~16 SU/h (2 h interactive link ~32 SU; 48 h batch ~770 SU = whole budget).
2. Build: hesnf_gpu (he_star_m1 + -ffp-contract=off, same recipe as 3896d529) @ 4de8168d (t2befb-ncfix-1cb) OK.
   /u/jma20/ATHENAK/builds/bin/athena_hesnf_gpu_4de8168d2b23 md5 c56693cbde14b001128e39189ffa5398
   .comment: GCC: (SUSE Linux) 13.3.1 20250313 [revision 4ef1d8c84fae...] ; GCC: (SUSE Linux) 14.2.0
   Binary contains the implicit_cs2_clip key. Not used for production until the fin51-r-fix tip NOTE.
3. Deck: the EXACT template docs/handover/rsg-E-1020/rsg_E_start1020r_clip.athinput.in (md5 4c442452, verified), copied
   unchanged to /work/hdd/bivj/jma20/rsg_E_1021/; substituted deck rsg_E_clip.athinput. diff substituted vs template = the 3
   @FILES@ lines only (he_ic_file, he_opac_table, he_planck_table -> /work/hdd/bivj/jma20/rsg_E_1021/files/...; all exist).
   Layout 2 nodes x 4 GH200 = 8 ranks, 3 blocks/GPU, mesh/meshblock untouched. Not started.
4. Interactive 2-node estimate is < 2 h (~40 min out): a 15-min N=128 memory/timing smoke is possible on request
   (~8 SU incl. setup). Send the binary to use (4de8168d or the fin51-r-fix tip).
Race rule noted: a DeltaAI start never wins; other copies keep running.
