# TASK Caltech -> Raven (USER, 11-04): BSG opacity tables and IC files from the ic_star BSG work, for the Caltech BSG wedge setup

Caltech is setting up a new BSG run (L. Ma / Ma+2026 MESA model, M 20 Msun, L 5.581461e38, X 0.7 Z 0.008 GS98) on the
RSG E / He giant global production method: tools/ic_star mode mesa, exact HSE + n2 seed, adaptive MLT closure with the auto
ramp, M1 vet_gd twin off, general EOS. The USER wants OUR OWN opacity tables, i.e. the SAME kind as the RSG E run:
- Rosseland = MESA kap via ic_star opacity.mesa_kap (as rsg-g25-prod/rosseland_mesakap_4c0e8b64d26a.txt, HG opac_mk_afe2dd8da2bd);
- Planck = ic_star ext2 (TOPS + Ferguson 2005 blend + Saha-es / power-law low-rho extension; as RSG
  planck_ext2_x0.653848_z0.019757_d38e5526afef.txt).
Caltech has none of the raw inputs: no ic_star_cache (no TOPS raw, no ferguson05, no MESA kap_data.tar.xz, no MESA install)
and /raven is not reachable from here. Caltech is trying a LANL TOPS fetch (planck_tools/fetch_tops.py) as a fallback for the
ext2 Planck part only; the MESA-kap Rosseland table needs Raven (the kap_data tarball) in any case.

## Please push (to this branch, new dir docs/handover/bsg_ext2_1104/; gzip anything > 5 MB) or give a fetchable path, with md5s

1. Tables for X 0.7, Z 0.008 (GS98, a_metal 16), built with the SAME ic_star recipe/version as the RSG E tables (state the
   commit and the cache hash; include each provenance.json):
   a. the MESA-kap Rosseland table (opacity.mesa_kap with kap_data.tar.xz md5 76857b1d386f411100d07746028ae4be, the RSG join
      log T [6.3, 6.5], X 0.7, Z 0.008); if it was never built for this composition, please build it
      (`<cache>/tables/opac_mk_<hash>/`);
   b. the ext2 Planck AND ext2 Rosseland tables with the TOPS low-density file (RESULTS sect. 8: cache opac_5b8e2fdfefea);
      also the older ext2nolo pair of bsg_newsetup_1015 if it is still there (README check: X 0.7 Z 0.008 rebuilt = ext2nolo,
      every value);
   c. the raw TOPS files used: /raven/ptmp/jinma/ic_star_cache/raw/tops/bsg_x0.7_z0.008_hi.dat and _lo.dat (or the
      `<cache>/raw/tops/<hash>_{hi,lo}.dat` copies), and the ferguson05 ross/planck files of `<cache>/raw/ferguson05/` that
      the X 0.7 Z 0.008 build read (only those files, if the directory is large).
2. The BSG IC and decks of bsgres_1017 / ic_star (RESULTS sect. 9, 9a):
   - ic_bsg.txt (seam-fixed IC, R_ph 3.5877e12, floor 1e-15 rule) and the IC tool outputs of that build (deck_bsg.athinput,
     diag_bsg.txt, prov_bsg.json) from /raven/ptmp/jinma/ic_star_runs/bsg (or wherever the 1e-15 control rebuild lives);
   - /raven/ptmp/jinma/bsgres_1017/decks/bsgres_r496.athinput (f14ecdaa), bsgres_r808.athinput (88ce6a1e), bsgres_g480;
   - the seam deck /raven/ptmp/jinma/bsg_newsetup_1015/decks/bsg_s0p_std_g480.athinput.
3. Result summaries (text only, no dumps): bsgres_1017/NOTE.md, ana/aprime_r496_r808.txt/json and the other ana/aprime_*.txt,
   ana/cmp_*.txt, and the bsg_P10 10-d gate summary (RESULTS sect. 7a).

## For comparison, what Caltech built meanwhile (same pipeline, explicit Ma+2026 TOPS tables, taper off, startup hopf)
- spec /resnick/groups/carnegie_poc/jingze/bsgsp_1103/ic/bsgsp_noseam.json (= specs/bsg.yaml + the RSG/HG startup):
  two-sided fit R_ph 3.57860e12 (yours 3.5786e12), floor 1e-15; seam fixed point 3.58854e12 (2 passes, -0.08 cell).
- r496 grid regenerated from regress/bsg_grid_refine.py: nx1 496, plateau dr 0.0319 Rsun, fit rms 0.012 cells.
Please reply with a NOTE in the same dir (NOTE-2026-11-04-raven-bsg-ext2-tables.md): file list + md5s + the build commit.
