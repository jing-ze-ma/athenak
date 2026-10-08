# TASK for DeltaAI: AG Car A and B with the sp BLEND radiation flux (supersedes the central-flux TASK-2026-10-08-deltaai-agcar)

User 10-08: run the AG Car shake-downs where they start fastest. DeltaAI GH200 is ~4x a viper node. Copies are queued
on Raven (A 31011799, B 31011916) and viper (A 12135012); first copy to START wins -- **push a NOTE the moment a DeltaAI
copy starts** (viper cancels its pending twins). Do not start more than one copy of each.

Why blend: the central implicit flux converges to a wrong state in the optically thin atmosphere in multi-D (L above R_ph
collapses to 0-30 % with two-way saturated beams); `implicit_flux = blend` (sp wedge, x1 faces) fixes it, interior
identical (gate: viper /viper/ptmp2/jinma/spblend_1008/gate2/GATE2.md; cost +41 % per cycle).

1. Code: branch **sp-blend-1008 @ 4622f424** (pushed; = rt-integration 37f3d1bd + the sp blend + m1_fs output).
   Build PROBLEM=he_star_m1 for GH200 as in TASK-2026-10-07-deltaai-hegiant sect. 2 (incremental from your he_star_m1
   build dir). Record md5.
2. Files: this branch (`he-ic-eint-from-t`), `docs/handover/agcar-files-1008/` (inputs now: A blend + tfloor 3000 K,
   tlim 8.064e6 s; B blend, tlim 1.3e6 s; MD5SUMS updated) -> `bash SETUP.sh <dir>`.
   **Important:** a RESTART from a central-flux state switched to blend diverges (first solve FATAL) -- run FRESH starts only.
3. Smoke `time/nlim=10` per case, 1 node x 4 GH200, 4 ranks (one 480x64x64 block per GPU). References (A100/MI300A,
   same code; expect ~1e-12, not bitwise):

| | A (viper 12134776) | A (Raven 31011349) | B (Raven 31011830) |
|---|---|---|---|
| t | 1.2525591814773279e+04 | 1.2525591814773266e+04 | 2.4186988199012517e+03 |
| mass | 6.1192879615026949e+32 | 6.1192879615022597e+32 | 1.7891809043881897e+31 |
| tot-E | 1.6552578941765284e+47 | 1.6552578942684888e+47 | 1.0633773267726825e+46 |
| Picard mean/max, NON-CONV | – | 13.6 / 17, 0 | 27.5 / 200, 1 (resid 2.3e-7 at the top cell i=480: accepted) |

   Pass: rc 0, 0 FATAL, NaN 0, numbers within ~1e-10 relative.
4. Production if both pass: A (tlim 8.064e6 s) and B (tlim 1.3e6 s) as chained 1-node jobs with restarts
   (restart output every 5e5 s A / 5e4 s B already in the inputs). Report job ids, start time (relative), s/cycle in
   NOTE-2026-10-08-deltaai-agcar-blend.md on this branch.
