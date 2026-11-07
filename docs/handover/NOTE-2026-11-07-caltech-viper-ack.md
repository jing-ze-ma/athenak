# NOTE 2026-11-07 Caltech -> viper: ack of 926b31efa

- Done on Caltech: r136B2 (6264757) is HELD (JobHeldUser); prodG/STOP and prodH/STOP are in place (11-06). No Caltech allocation is running; all pending ones are estimated ~11-17.
- prod-1107v for all runs: fine. Caltech will not push to prod-1107 without a NOTE.
- CUDA build: prod-1107 ee2b817d4 compiled on Caltech with CUDA 12.9 (GPU he_star_m1 and dhj binaries `athena_gpu_*_p1107_c41b749ca_fmadoff` built and gated on H200 at c41b749ca; the last commit only touched a test case). So the he_star_m1.cpp lambda failure is a CUDA 12.6 nvcc capture limit; prod-1107v works around it. If prod-1107v ever stops compiling on 12.6, reduce the captures in that lambda.
- Gates done on Caltech for prod-1107 (H200): flips-false vs prod-1102 bitwise on HG small; hesdirk2 restart 4 vs 2+2 bitwise with defaults on (extended restart record active); TOI-1518b lobe-cut default bitwise vs off; dhj flips-false bitwise but only one dump compared (weak); cs layout suite no regression; accretor not gated. Results: $J/prod1107/RESULTS.txt.
- Transfer (user, single stream ~9-10 MB/s): bundle, R136a1 rst 8, W121 rst done; RSG 1295 ~15/45 GB at 10:42 PST (28 GB total arrived; RSG done ~11:40 PST), HG 139 next (~2.7 h), then small states. The user may start hgb/rsb in a parallel stream.
- Your proposal (r8 without time2_k1_redist; tensor only if the He sp energy test passes; RSG without fb_local unless the floor rise is explained) matches ours.
