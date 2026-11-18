# NOTE 2026-11-18 Caltech -> viper (USER): newest full-sphere 3-D accretor page (answer to QUESTION 3c5e05596)

## 1. The page

| | |
|---|---|
| Title | **RY Per Accretor 3-D** |
| URL | https://claude.ai/artifact/D8q681tiTgNDnBM6aaJHz3 (Caltech's pre-11-02 account; returns 404 from the current account too, so it cannot be shared) |
| Last update | 2026-11-04 (backup 11-05; the run was not refreshed after that) |
| Main run shown | **spherical-polar full-sphere PRODUCTION**, 500 x 116 x 512 (r x theta x phi), stretched theta with a fine +-20 deg core (58 x 0.703 deg, 29-cell caps per side) |
| Run dir (Caltech) | `/resnick/groups/carnegie_poc/jingze/accsp_1102/prod/fresh/{static,stream}`; launch `accsp_1102/prod/prod_sp_fresh.sh`, resumes `prod_sp_resume.sh`, `prod_sp_resume_accpolar.sh` |
| Physics | 3-D Roche potential, WB tangential tide, Coriolis + centrifugal, implicit M1; static to t 0.01, then 1e-4 Msun/yr L1 stream; polar_ring_avg_gz |
| Binary | launch: `builds_accr/bin/athena_gpu_ry_per_accretor_r2` (ringsmooth-1102 e67c9e967); since the 11-04 02:26 resume: `builds_accr/bin/athena_gpu_ry_per_accretor_prod_6854e53d_fmadoff` (md5 18175a38), extra keys `rad_m1/vet_gd_band_global=true mesh/sp_phi_global=true ...` (full list in `accsp_1102/prod/fresh/prod.log`) |
| Last state on the page | t = 0.1436 code = 27.75 h = 0.3035 orbits, cycle 37450, dt 3.26e-6, ~2.1 s/cycle on 4 H200 |
| Last restart | `accsp_1102/prod/fresh/stream/rst/ryper.00006.rst` (run killed with the allocation, not resumed since) |

Also on the same page, kept at the bottom: the cubed-sphere N = 128 test (`csacc_1031/test128wb`, shown 11-02) and RUN12
(the 10-24..10-30 sp wedge). Other full-sphere accretor pages on the same old account: "Direct Start Stream" P8oPmat5FA2ijxy7hJoK3Q,
"Spin Ramp" 5dYxRviTGLYUnS4a3FxrG3, "Spin Ramp 1.711" HeGZZaXj2mEFC6nkp4Nf7N, "Spin 1 Stream" 6woQzZv4hH28ZJiEuYACbQ (all
fullsph_1027 spin runs, 10-29..10-31, older than the production page).

## 2. Files to republish (pushed, 9.8 MB)

`docs/handover/caltech-to-viper-1118/acc3d_sp_page/`:
- `index.publish.html` (the page, md5 97e4d123) + `img/` (16) + `movie/` (12 mp4/poster) = the exact published file set;
- `republish.json`: one Artifact publish call (title, icon, files map) with paths rewritten to this repo directory; run it from the repo root.
  It creates a new URL on the current account;
- `refresh.sh.caltech` + `src/` (build_page.py, mk_figs.py, job.sh, ...): the Caltech page builder (paths are Caltech's; needs the
  run dir above, so it only refreshes at Caltech).

-- Caltech
