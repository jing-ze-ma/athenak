# TASK 2026-10-21 raven -> caltech: queue the He giant GLOBAL (HG global) copy (USER)

Same run as Raven HG R / Viper HG V (fresh t=0, cubed sphere N=128, 704 radial, 24 MeshBlocks 704x64x64).
- Binary: your fin51-r-fix 4de8168d he H200 build (9a56b4a1) — the globals stay on 4de8168d (USER).
- Deck: docs/handover/hgglob-1020/hgglob_fresh_clip.athinput.in (md5 4e6ea902) = Viper HG V deck byte-for-byte except the 3 @FILES@ lines
  (ic / rosseland / planck, in docs/handover/hg145-prod-1018/, md5-identical to Raven/Viper). rst every 5 d, tlim 1e9, longest walltime.
- Layout: 2 nodes x 4 H200, 3 blocks per GPU (~59 GB/GPU by the memcut plan; 1 node x 4 H200 at 6/GPU ~120 GB is too tight).
- Use the UCX setting from TASK 1ba37054 (RSG E nondeterminism / IB fault): do NOT let link 1 start before that test picks the env;
  queue now (priority age) but hold link 1 until then. Same STOP rules as RSG E. NOTE the job ids + run.cfg md5.
