# xthinfix-1010 on Raven (taken over from viper, which is down since ~07:20 CEST 10-10)

## Context
- Goal (user 10-09): the half-range implicit VET face-flux scheme (implicit_flux = blend, implicit_flux_faces = all,
  implicit_flux_beam = halfrange, idort_f, implicit_blend_xthin) is NOT default because the xthin = 30 override
  destroys convergence. Branch fork/xthinfix-1009 tries fixes. Productions name the hr keys explicitly.
- Caltech results on 413b38af (NOTE-2026-10-10-caltech-arm-matrix / -steep-battery on fork/bsg-files-1009):
  - sq (steep X0 60 + plm + qs 0.1) has the right convergence order (thin pulse, rw_t10, atm in every region),
    but diverges on beams ba0 (cycle 4) and ba20 (cycle 88), costs ~4x on GPU (Picard 45-47 vs 4.3),
    and its plm path segfaulted in 2-D G1/G3 (generic halo path; fixed in e50ec7d6).
  - dc kn arms are cheap but do not converge in the thin top (Kn > 0.3); steep dc breaks rw_t10 / atm.
- viper pushed, untested anywhere:
  - 81cc7f95 `implicit_hr_pos = bound` (per-face Zalesak/TVD bound test, sig = 0 where violated; instead of E<=0 kill)
  - a7ca3f65 `implicit_hr_recon_qs_mode = step` (plm switch + limiter built once per step from E^n / previous step's
    change; operator fixed through Picard) -- aimed at the Picard cost
  - e50ec7d6 plm halo fix + `implicit_hr_recon_pred = true` (dc predictor per step)
- Pending viper->Caltech task: rerun plm-arm G1/G3 gates on e50ec7d6 (viper's own sq check: G1 cfl 1e2 rc 0,
  max|dE/E| 1.74e-2, L1 6.5e-4; G3 cfl 10 rc 0, L1 0.0138).

## Acceptance criteria (user, 10-10)
Goal: an accurate AND fast implicit VET.
1. Second order in space and time overall (X, T and C ladders ~2). First order in SPACE is acceptable only in the
   optically thin part (the thin criterion is one of the variants: xthin modes all/steep/kn, Knudsen Kn0, etc.).
   So: thick / intermediate regions need X order ~2; thin regions (tau < ~0.1 or Kn > Kn0) need X order >= ~1
   (orders ~0 or negative, or errors growing with resolution, fail). T order ~2 everywhere.
2. Passes the simple tests (pulses, rw, atm vs Hopf, G1/G3/G5 gates, beams with no divergence / NON-CONVERGED).
3. Efficient on complex simulations: GPU AG Car A cost close to production `all` (Picard ~4-5, not 4x).

## Workspace
W = /raven/ptmp/jinma/xthinfix_1010
- src/ = git archive of fork/xthinfix-1009 e50ec7d6 + kokkos 4.6.2 (d8e9af03). Do not edit src/.
- xtf/ = docs/handover/xthinfix-1009 from fork/bsg-files-1009 (battery scripts, ana/od.py, inputs).
- Reference Caltech scripts/results: fork clone /raven/u/jinma/ATHENAK/athenak-fork (origin/bsg-files-1009).
- viper's Raven templates: /raven/ptmp/jinma/xthinfix_1009 (build.sh, gpu2.sh), /raven/ptmp/jinma/vgdfuse_1010/raven_time.sh.

## New arms (sq = steep X0 60 + plm + qs 0.1, as in battery_arms.sh)
| arm | keys on top of sq |
| --- | --- |
| sqs | implicit_hr_recon_qs_mode=step |
| sqb | implicit_hr_pos=bound |
| sqp | implicit_hr_recon_pred=true |
| sqsb | step + bound |
| sqsbp | step + bound + pred |
| aqs | aq (mode all + plm + qs 0.1) + step |
References rerun on the same binary: cen, hr, sq, aq (and knp in the gates).

## Rules
- No git push (no credentials on Raven). Never touch the held agcA_rc / agcB_rc jobs or other dirs' runs.
- Short jobs only. Write results into W; do not write into the fork clone.
