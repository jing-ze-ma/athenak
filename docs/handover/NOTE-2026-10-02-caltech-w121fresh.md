# NOTE 2026-10-02 (Caltech): WASP-121b fresh 1x / 3x / 10x (face 6 + wall fix) SUBMITTED

For TASK-2026-10-02-caltech-w121fresh.md. Status: QUEUED, not started. A follow-up NOTE is pushed when each arm STARTS;
viper then cancels its copy of that arm (f1x6w 12059646->47, f3x6w 12059648->49, f10x6w 12059650->51->52).

- Dir: /resnick/groups/carnegie_poc/jingze/w121fresh_1002/ (in/, 1x/, 3x/, 10x/, smoke/, prod.sub, rst_args.py, JOBS)
- Layout: each arm 1 node x 2 H200, 2 MPI ranks (1 per GPU), same as w121prod_0928; --exclude=hpc-sm-01-09,hpc-sm-02-16
- Code: rt-integration f24fcfcf (>= 3f3b0414), CUDA build job 3761762 (build_inc.sh, Hopper90, Kokkos pin);
  binary /resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_deep_hot_jupiter_rt_f24fcfcf
- Jobs (link 1 afterok build 3761762; link 2 afterok link 1; 24 h each):
  | arm | link 1  | link 2  |
  |-----|---------|---------|
  | 1x  | 3761953 | 3761954 |
  | 3x  | 3761955 | 3761956 |
  | 10x | 3761957 | 3761958 |
- md5: MD5_IC all 3 OK. MD5_CK: every file the inputs read is OK (1x exo_fms_ck: ck/CE/cia/ray/sw_flux all OK;
  10x ckdata10: all 6 OK; 3x: ck/Premixed_3x_g8_11_hiT2 + CE_tables/FastChem_ck_3x_int_hiT2 OK, copied from the
  09-30 3x package into w121fresh_1002/ckdata3 with cia/ray/sw_flux linked to the verified 1x dirs). Not present
  locally (unused by the inputs): 3x non-hiT2 tables (Premixed_3x_g8_11{,_hiT}, FastChem_ck_3x_int{,_hiT}) and
  CK1X/tools/__pycache__.
- Inputs: bundle files with only @BUNDLE@/@CK1X@/@CK3X@/@CK10X@ filled; nothing else changed.
- prod.sub: fresh start runs an in-job smoke first (smoke/<arm>: leg A fresh -t 5 min, leg B restart from its newest
  rst with last_time reset -t 2 min); gate = rc 0, 0 FATAL/nan/NOT-CONVERGED, startup line
  "ck_sph_face = 6 (VEF), effective: ck_sph_top = 1, ck_vef_every = 150, ..., ck_impl_rowfb = 1, ck_impl_rsec = 4
  (10x: 20, explicit in the input), ck_wall_flux_exact = 1", leg A median dt >= 2 s (logged next to viper's
  15.2 / 12.9 / 12.6 s). Restart links: sanity check on the previous link (refuse on NaN/FATAL, no newer rst,
  median dt < 2 s; writes STOP_<arm>), then -r newest rst with every <outputN>/last_time = floor(t/dt)*dt
  (rst_args.py; t read 232 bytes after the "<par_end>\n" line). A link whose rst is already at tlim exits 0.
- Expected speed: ~14-15 ms/cycle (09-28 w121prod on 2 H200; the 1x 300 rot took ~10.5 h in one link), i.e.
  likely 1-2 links per arm; the 10x multi-rotation check (check_xl.py: per-rotation dt median >= 5 s, KE < 1e34,
  |Etot/Etot0-1| < 1e-2) is done at the ~50-rotation reports.
