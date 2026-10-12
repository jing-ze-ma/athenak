# NOTE DeltaAI: vet_gd_rotate_slope = correct GPU cost, full AG Car A (TASK-2026-10-12-deltaai-rotslope-cost) -- DONE: PASS (FINC avg / FIN avg = 1.041)

**Picked up** 09:56 CDT. rotslope-1012 **1bfe213a** built: `athena_hes_gpu_1bfe213aa831`, md5 cdacdee5234e405cdc6ca3f075e24afc.
Decks filled with sed only (FIN/FINC differ in the 1 line vet_gd_rotate_slope keep/correct); files md5 4ab9e3f3 / b28e97c5 / 493f9d01.
Runs in ONE ghx4-interactive job **3366761** after TASK defaults (2 x 60 cycles, ~5-10 min), then FINC then FIN, 120 cycles each, commands as in the TASK.
Dir /work/nvme/bivj/jma20/rotslope_1012. STARTED line follows when it runs.

**STARTED**: job 3366761 RUNNING on gh143 since 10:18:55 CDT (defaults part first, ~5-10 min; then FINC, then FIN).

## RESULT (posted ~10:25 CDT): job 3366761, gh143, rotslope part 10:20:30-10:22:18 (job COMPLETED, 3 min 23 s total)
**Verdict: PASS.** Both rc 0, no FATAL or NaN, 0 NON-CONVERGED, 5 correction solves with 0 failed; **FINC avg / FIN avg = 0.2251 / 0.2162 = 1.041** (<= 1.05).

| arm | rc | FATAL/NaN | wall | tana.py run.log 2 20 |
|---|---|---|---|---|
| FINC (correct) | 0 | 0 / 0 | 55 s | ncyc=120 steady=0.2187 (n=112) rot=0.3697 (n=5: 0.380 0.381 0.366 0.377 0.344) avg=0.2251 (n=117) init=20.27 |
| FIN (keep) | 0 | 0 / 0 | 53 s | ncyc=120 steady=0.2133 (n=112) rot=0.2810 (n=5: 0.293 0.282 0.281 0.274 0.274) avg=0.2162 (n=117) init=19.30 |

Picard:
- FINC `implicit transport: solves=248 Picard iterations mean=3.024 max=6 NON-CONVERGED=0`
- FIN  `implicit transport: solves=243 Picard iterations mean=2.893 max=6 NON-CONVERGED=0`

hesdirk2 (FINC; FIN's line is identical up to the slope part, which it lacks):
`stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 time2_vet_col=predict fallbacks=0 source=blend x0=1 time2_pred_cap=0.5 limited cell-builds=0 of 9.358541e8 vet_gd_rotate_slope=correct solves=5 failed=0 (frac 2.928932e-01, reuse on) rotate_align shortened steps=122`
Both end at t = 1.5e5 s, cycle 120.

Notes:
- **N = 5, not 6:** in 120 cycles the rotation steps are at cycles 20, 40, 60, 80 and 100. The rotation at cycle 120 would fall on the step after nlim, so 5 is consistent with the TASK's own count (+5 solves: 248 vs 243).
- Rotation cycles cost +0.089 s (0.370 vs 0.281, x1.32), i.e. one extra Picard solve, as expected.
- Steady cycles are also 2.5 % slower in FINC (0.2187 vs 0.2133). This may be run-to-run noise (one run each), or the post-correction steps may take more Picard iterations: the overall Picard mean is 3.02 vs 2.89, more than the 5 extra solves alone would explain. Not resolved here; per-cycle Picard counts are not in the log.
- `rotate_align shortened steps=122` > 120 cycles in both arms (the same in both); reported as printed.

Build: rotslope-1012 1bfe213a, athena_hes_gpu_1bfe213aa831, md5 cdacdee5234e405cdc6ca3f075e24afc (both arms).
Logs: rotslope-1012/deltaai-run/{FINC,FIN}_run.log. GPU-h: ~0.13 (1.8 min x 4 GH200); the whole job (both TASKs) used 3.4 min x 4 = 0.23 GPU-h.
