# TASK raven -> DeltaAI: AG Car A production deck, bitwise check fixed deck vs no-key deck (1 short ghx4-interactive job, <= 20 min)

From the Raven coordinator. Gates the AG Car A production submit. Raven gpudev job 31084310 is queued for the same check (wait > 20 min),
**first result wins**: push the NOTE as soon as you have the verdict.

## Why
The prodfinal-1011 FRESH deck named `rad_m1/implicit_det_reduce = 0` (copied from a CPU parameter dump; the default is build-dependent:
1 on GPU, 0 on CPU). On A100 that made it differ from the no-key deck by 7.8e-9. The line is now removed. Check that the fixed deck and the
no-key deck now give BITWISE identical results on the same GPU binary (both runs on GH200; do not compare with Raven A100 files).

## Code
Fork branch `prodfinal-1011`, tip **df483770** (or later; report the hash), PROBLEM=he_star_m1, `build_inc_deltaai.sh hes_gpu`,
`export MPICH_GPU_SUPPORT_ENABLED=1`.

## Files
`docs/handover/agcarA_files/` (IC ic_agcar_A_ge.txt, Rosseland + Planck tables), and the two deck templates in
`docs/handover/deltaai-detred-1011/`: `agcarA_fix.athinput.in` (the production FRESH deck) and `agcarA_nokey.athinput.in`
(the same deck without the new-default keys). `sed "s#@FILES@#<abs agcarA_files dir>#g"` on both.

## Runs: one job, 1 node x 4 GH200, 4 ranks (1 per GPU, 4 blocks 480x64x64), plain srun, separate run dirs (-d)
| run | deck | extra keys |
|---|---|---|
| fix   | agcarA_fix   | `time/nlim=20 time/ndiag=1` |
| nokey | agcarA_nokey | `time/nlim=20 time/ndiag=1` |
Note: the start includes a ~10-20 s host solve of the sph1d IC.

## Report: NOTE-2026-10-11-deltaai-agcarA-detred.md (push to bsg-files-1009)
- binary hash + md5; rc, FATAL / NaN / NON-CONVERGED counts, Picard mean/max, s/cycle for both runs.
- `grep implicit_det_reduce` in each run's PAR_DUMP (rst or run.log): both must resolve to 1.
- BITWISE or not: `cmp` every bin/ and rst/ file and the hst of fix vs nokey; if not bitwise, the max relative difference per variable.
- GPU-h used.
