# NOTE 2026-11-11 ORION -> viper: vgdh-sw COMPLETE (sp part) -- the wider smoothstep also helps sp: pulse_u full C
# CAND 1.14 0.45 -> sw3 1.31 1.53 / sw4 1.31 1.41 (F1 0.56 -> 1.73 / 1.92); sp lat twins T no regression

Completes `TASK-2026-11-11-viper-orion-vgdh-sw.md` (interim 32e196f5c = cs + beams; sw5 in 7ffd3cade / 2b93de0e4).
sp on fork/vgdlay-1110 **453f879bd**, pulse_u SP_FULL with ONE angular MeshBlock per level (SP_MBN=1: layout fixed in
physical angle, 1 rank per run), nx1 64 128 256 512 (angular nx2 = nx3 = 4 8 16 32). All runs rc 0, no FATAL.

## pulse_u full refinement C (E / F1 p of 64-128-256, 128-256-512)
| arm | E | F1 |
|---|---|---|
| CAND (453f879bd, same layout) | 1.14 0.45 | 0.87 0.56 |
| CAND-wsx0 | 1.24 0.61 | 1.02 0.80 |
| CAND-wc5 | 0.73 1.44 | 0.39 1.90 |
| **CAND-sw3** | **1.31 1.53** | 1.05 **1.73** |
| **CAND-sw4** | **1.31 1.41** | 1.16 **1.92** |
(sw5 SP_FULL still running; the remaining vgdh / geo arms on this layout too.)

## sp lat twins T (5 dt levels)
| case | CAND | sw2 | sw3 | sw4 |
|---|---|---|---|---|
| lat_th_atm E | 1.93 1.85 1.75 | 1.95 1.92 1.92 | 1.94 1.92 1.93 | 1.96 1.92 1.89 |
| lat_th_atm F1 | 1.53 1.55 1.52 | 1.57 1.60 1.62 | 1.55 1.65 1.58 | 1.56 1.61 1.60 |
| lat_th_atm F2 | 1.39 1.42 1.63 | 1.38 1.40 1.56 | 1.27 1.49 1.47 | 1.39 1.41 1.62 |
| lat_th_mid E / F1 / F2 | 2.02 2.01 2.01 / 2.03 2.02 2.01 / 1.99 1.99 2.00 | = | = | = |
(For reference on the same binary: CAND-w1 lat_th_atm E 0.43 0.66 0.77, CAND-wc5 1.05 0.83 0.91 -> both FAIL on sp too.)

## Reading
On sp the wider smoothstep is also the right direction: the finer pulse_u pair goes 0.45 -> 1.4-1.5 in E and 0.56 -> 1.7-1.9 in
F1, without costing the atmosphere twins (lat_th_atm slightly better than CAND in E). sp is not yet at 2 (the 64-128-256 pair
1.31): the residual is either pre-asymptotic (the coarse angular grid 4x4 at nx1 64) or the radial-face piece wsx0 also
shows (1.24). The cs picture (sw5: 1.90 1.99 on pulse_cs R3) is the cleaner one.
