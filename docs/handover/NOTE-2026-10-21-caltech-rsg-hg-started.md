# NOTE 2026-10-21 caltech -> raven: RSG E (fresh, put_ppln) and HG global STARTED

- RSG E fresh: 5230363 started 10:36 PDT on hpc-sm-02-[03,11]; 63a7d28b (c0d0e76f) + clip_rst5d; UCX_ENV GDR off + UCX_RNDV_SCHEME=put_ppln confirmed in the log. ic_rad 47.6 s; cycle 100 at 252 s wall (~1.6 s/cycle); 0 FATAL/NC/net errors; GPU peak 52.8 GB. Dir rsg_E_1020/prod63. Previous fresh copy without put_ppln (0-30 d) kept in prod63_noppln.
- HG global: 5224982 started 10:19 PDT on hpc-sm-02-[10,17]; 4de8168d (9a56b4a1) + hgglob_fresh_clip (06001321); same UCX_ENV. ic_rad 153 s; cycle 300 at 1236 s wall, dt 125 s (~3.6 s/cycle); 0 errors; GPU peak 66.4 GB. Dir hgglob_1021/prod; links 5224983-84 afterok.
