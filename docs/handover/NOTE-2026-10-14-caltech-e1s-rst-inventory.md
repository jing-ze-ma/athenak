# NOTE Caltech -> Raven (10-14 13:15 PDT): E1S rst / bin inventory (reply to 3e8f36eb / 7210fc8b)

Run dir: /resnick/groups/carnegie_poc/jingze/e1s_1014/prod (job 4513330 running, t = 2.21e7 s at 13:10, cycle ~30900).

- **rst**: ALL 45 still present, agcar3d.00000-00044.rst (prune did not remove any; 212 GB total, 3.3-5.4 GB each).
  Cadence output5/dt = 5.0e5 s -> rst N at t = N x 5e5 (00044 = 2.20e7 s, written 12:59). Next (00045, 2.25e7) in ~40 min.
  Newest at remap time = whatever is latest then; nothing extra kept (item 1 withdrawn, noted).
- **bin**: hydro_w and m1, output2/3 dt = 1.0e5 s, indices 00000-00221 (t = 0 to 2.21e7 s), 444 files.
  **float32** (header "size of variable=4"), location doubles (size 8).
- Layout: 536x64x64, 4 MeshBlocks, 1 node x 2 H200; deck agcA_E1_fin4dtr_s3_tlim1e9.athinput (680dadb0), build 6252e58e.
- Health: 0 FATAL. DT_RETRY rate rose from ~1 per 1000 cycles to ~25 per 1000 since t ~2.0e7 (Picard NON-CONVERGED at resid ~5e-8 vs tol 1e-8, recovered by a dt halving each time); dt still ~700 s.
