# TASK raven -> Caltech (2026-10-18 evening): 3-D ACCRETOR PRODUCTION v3 - HOLD b592b459 LIFTED for this code; START ASAP

Replaces TASK 7e992e66 (186306cb is NOT production any more). USER: start as early as possible. Race vs Raven: first site to
START wins; post link 1's start time at once.

## 1. Build
Branch **accretor-3d-stream-1018 @ a6624af6** (pushed) = 186306cb + the 2-D stream fixes (problem/stream_sector_light) +
problem/env_ic_vz (von Zeipel IC T in the atmosphere, g at the ray's photosphere). Same recipe as your 186306cb build
(PROBLEM=ry_per_accretor, HOPPER90, -ffp-contract=off). Raven A100 twin: athena_a100_ry3d_v10_nofma md5 125b6010.

## 2. Files: docs/handover/acc3d-prod-1018/
- NEW deck `acc3d_prod_v10.athinput.in` md5 5e342048 (sed @FILES@ as before; change no key) = the old deck +
  `stream_sector_light = surface` + `env_ic_vz = true`.
- unchanged: rosseland_mesakap_gs98_x0.70_z0.02.txt d81b7a46, planck_ext2_gs98_x0.70_z0.02.txt 99d668f6,
  column_mesakap_z0.02_p381.txt 5b813886, acc3d_chain.sh 1ad77aed.
- Startup line to report: `ry_per_accretor: env_ic_vz: ... factor a .. b` (expect ~0.97-1.03) and `stream_sector_light surface`.

## 3. Layout: your 2-H200 variant (passed your pre-check): 2 ranks, 4 MeshBlocks per GPU (deck unchanged), gpuwrap.

## 4. Pre-check: NOT required (restart bitwise verified on MI300A and by your 2-H200 run for this deck family); submit
link 1 with no dependency so it starts immediately.

## 5. Chain: SHORT backfill-friendly links
- Links of **12 h** walltime (WT = 12 h - 20 min = 41,999 s), submit **6 links** chained afterok + kill-on-invalid-dep
  (= 3 days; more later). acc3d_chain.sh unchanged: static phase problem/stream=false to t 0.05 (~0.5 s/cycle x ~3,600
  cycles = ~30 min on 2 H200), then the deck as is (stream on, tlim 1e9) from the newest rst; each link restarts from the
  newest rst; rst every 0.01 code time (2.98 GB), pruned IN-JOB to the 3 newest stream rst.
- STOP rule (in the script): only end-of-run summary NON-CONVERGED=[1-9], FATAL, "Picard solve DIVERGED", NaN, rc != 0.
  Expected at stream contact (t ~0.068, i.e. ~0.018 after stream on): dt falls ~1e-5 -> ~3e-6; a few DT_RETRY lines are
  normal (dt_min 1e-8 fallback in the deck) - never stop on those.
## 6. Report: NOTE-2026-10-18-caltech-acc3d-prod.md: binary md5, job ids, link-1 START time, the two startup lines; per link
rc, t, dt, s/cycle, DT_RETRY count, summary NON-CONVERGED, chain.txt tail.
