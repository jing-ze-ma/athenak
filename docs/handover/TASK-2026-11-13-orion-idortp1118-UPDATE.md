# TASK UPDATE 2026-11-13 viper -> Orion: idortp1118 arm q uses EXTRA_HG (USER 11-13 ~10:00 PST)

The user wants ONE uniform key set for all stars and, if it passes everywhere, to make the HG positivity keys default.
So the candidate default is EXTRA_HG (= EXTRA + rad_m1/implicit_vimp_lag_local=true rad_m1/time2_adm_poslim=true
rad_m1/time2_dtr_poslim=true). Change to TASK-2026-11-13-orion-idortp1118 (d2adb6fb6):
- parts (a) and (c): run arm q with EXTRA_HG instead of EXTRA (export EXTRA="$EXTRA_HG" before one_q.sh / the wave runner).
  If q runs with plain EXTRA already finished or are running, keep them (label them q0) and add the EXTRA_HG runs as q.
- part (b): unchanged (both P18 and P18H).
Pass rules unchanged; in addition list every row where q (EXTRA_HG) differs from d6int / q0 and whether the positivity keys
fired (grep the run.log for poslim / vimp_lag_local counters; report counts).
