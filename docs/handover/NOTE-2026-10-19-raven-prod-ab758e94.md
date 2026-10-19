# NOTE 2026-10-19 raven -> caltech: production branch fin51-r-fix = ab758e94; USE THE NEWEST PRODUCTION BINARY FROM NOW ON (USER)

## fin51-r-fix fast-forwarded a2713be0 -> ab758e94 (pushed 10-19 ~15:45)
ab758e94 = prodmerge-1019 74811de1 (a2713be0 + accretor v11 acc3d-v11-prod-1019 + m1lpr-1019 + hmstage-1019 (HIP host-staged implicit
halo; NVIDIA unaffected) + cs-open-1018-wip + dt_retry visitor regens) + amdtr-1019 (AmModule state in the dt_retry checkpoint).
Gates (Raven/Viper): strict-md5 gbit bsg/e1/hg BITWISE vs a2713be0; IDORT 20/20 BITWISE; cs RSG/hgglob old-vs-new, restart across
binaries, N==N2, restart BITWISE; accretor 74811de1 == 23e06019 (onset + contact); AM forced-retry / restart BITWISE.

## USER rule (10-19): always run the most updated production binary
- Rebuild your production binaries from fin51-r-fix @ ab758e94 (same build recipe as before; nofma for the accretor) and use them for
  every NEW link / job from now on (switch at link boundaries; record md5s in your NOTE).
- acc3d v11 copy (links 5007430-35, held for the USER): when released, use ab758e94 + deck acc3d_prod_v11am.athinput
  (= v11 + <hydro>/am_transport = true; Raven switches to exactly this at its link 2, ~10-20 06:30 CEST).
- RSG g25 and the global He giant are STOPPED/HELD (cs panel-edge instability in the deep layer, see TASK b2e1bfee). Do not start cs
  production. TASK e3b5c424 (RSG restart DIFF arms) still useful; TASK b2e1bfee (CSEDGE arms) has priority on the H200s.
- Raven He giant wedge (hg145, sp, unaffected) switches to ab758e94 at its link 2 (now).
