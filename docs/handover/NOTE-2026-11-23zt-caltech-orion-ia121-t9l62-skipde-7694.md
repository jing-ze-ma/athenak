# NOTE 2026-11-23zt Caltech -> orion: t9l62 IA12.1 skip_de jump 1.11e41 at cycle 7694, with a fallback in the same window (per 23zn)

Per-cycle increment of the floor_audit [update] skip_de: **1.111e41 at cycle 7694**.
floor_audit + CLEAN WARNING + DT_RETRY lines, cycles 7688-7700:
```
<hydro> floor_audit cycle 7688 [update] cell-c2p 1.669668536e+11 e_c<=0 11682 kedom 16581 sub_tfloor 45844118 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 45844118 (kedom 11711) floor_de 1.018516794e+42 <r> 5.487004241e+13 <r>_de 2.484764705e+12 skip(
<hydro> floor_audit cycle 7688 pay: paid 1.367375988e+43 unpaid 0 given_to_rad 2.954361979e+42 cells 25107680 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 7689 [update] cell-c2p 1.671010714e+11 e_c<=0 11740 kedom 16704 sub_tfloor 45871200 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 45871200 (kedom 11769) floor_de 1.018523206e+42 <r> 5.487000994e+13 <r>_de 2.484765029e+12 skip(
<hydro> floor_audit cycle 7689 pay: paid 1.367398995e+43 unpaid 0 given_to_rad 2.954398794e+42 cells 25121389 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 7690 [update] cell-c2p 1.672352891e+11 e_c<=0 11802 kedom 16830 sub_tfloor 45898276 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 45898276 (kedom 11832) floor_de 1.018628865e+42 <r> 5.486997193e+13 <r>_de 2.484772031e+12 skip(
<hydro> floor_audit cycle 7690 pay: paid 1.367472741e+43 unpaid 0 given_to_rad 2.954433201e+42 cells 25135095 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 7691 [update] cell-c2p 1.673695068e+11 e_c<=0 11867 kedom 16955 sub_tfloor 45925343 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 45925343 (kedom 11898) floor_de 1.019551814e+42 <r> 5.486993088e+13 <r>_de 2.484917772e+12 skip(
<hydro> floor_audit cycle 7691 pay: paid 1.368828422e+43 unpaid 0 given_to_rad 2.954467178e+42 cells 25148803 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 7692 [update] cell-c2p 1.675037245e+11 e_c<=0 11934 kedom 17078 sub_tfloor 45952405 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 45952405 (kedom 11965) floor_de 1.021122165e+42 <r> 5.486988873e+13 <r>_de 2.485218629e+12 skip(
<hydro> floor_audit cycle 7692 pay: paid 1.369038337e+43 unpaid 0 given_to_rad 2.95450092e+42 cells 25162504 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 7693 [update] cell-c2p 1.676379423e+11 e_c<=0 12000 kedom 17205 sub_tfloor 45979453 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 45979453 (kedom 12031) floor_de 1.021226346e+42 <r> 5.486984773e+13 <r>_de 2.485233509e+12 skip(
<hydro> floor_audit cycle 7693 pay: paid 1.369081952e+43 unpaid 0 given_to_rad 2.954536876e+42 cells 25176200 | fofc_retest: retested 0 unflagged 0
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12057 kedom 17331 sub_tfloor 46006484 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006484 (kedom 12088) floor_de 1.021319645e+42 <r> 5.48698171e+13 <r>_de 2.485244276e+12 skip(dfl
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
CLEAN WARNING cycle 7694: backward-Euler solve not converged (dt_retry, 19 so far; last change 4.756e-04)
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12021 kedom 17251 sub_tfloor 46006395 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006395 (kedom 12052) floor_de 1.021274924e+42 <r> 5.486985858e+13 <r>_de 2.485238909e+12 skip(df
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
CLEAN WARNING cycle 7694: backward-Euler solve not converged (dt_retry, 19 so far; last change 2.128e-03)
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12021 kedom 17244 sub_tfloor 46006361 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006361 (kedom 12052) floor_de 1.021251293e+42 <r> 5.486985897e+13 <r>_de 2.485236295e+12 skip(df
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
CLEAN WARNING cycle 7694: backward-Euler solve not converged (dt_retry, 19 so far; last change 3.360e-03)
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12022 kedom 17245 sub_tfloor 46006348 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006348 (kedom 12053) floor_de 1.021240884e+42 <r> 5.486985807e+13 <r>_de 2.485235338e+12 skip(df
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
CLEAN WARNING cycle 7694: backward-Euler solve not converged (dt_retry, 19 so far; last change 1.168e-02)
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12017 kedom 17243 sub_tfloor 46006333 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006333 (kedom 12048) floor_de 1.021233477e+42 <r> 5.486986376e+13 <r>_de 2.485234398e+12 skip(df
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
CLEAN WARNING cycle 7694: backward-Euler solve not converged (dt_retry, 19 so far; last change 6.445e-03)
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12016 kedom 17241 sub_tfloor 46006327 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006327 (kedom 12047) floor_de 1.021231398e+42 <r> 5.4869865e+13 <r>_de 2.485234038e+12 skip(dflo
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
CLEAN WARNING cycle 7694: backward-Euler solve not converged (dt_retry, 19 so far; last change 1.678e-04)
<hydro> floor_audit cycle 7694 [update] cell-c2p 1.6777216e+11 e_c<=0 12014 kedom 17238 sub_tfloor 46006322 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46006322 (kedom 12045) floor_de 1.021228962e+42 <r> 5.486986733e+13 <r>_de 2.485233773e+12 skip(df
<hydro> floor_audit cycle 7694 pay: paid 1.38020279e+43 unpaid 0 given_to_rad 2.954573222e+42 cells 25189882 | fofc_retest: retested 0 unflagged 0
CLEAN WARNING cycle 7694 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 20)
<hydro> floor_audit cycle 7695 [update] cell-c2p 1.679734866e+11 e_c<=0 12027 kedom 17264 sub_tfloor 46046455 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 46046455 (kedom 12059) floor_de 1.021234172e+42 <r> 5.486990657e+13 <r>_de 2.485233983e+12 skip(
```
Last cycle line: elapsed=5.050449e+03 cycle=7700 time=4.840446e+05 dt=6.319560e-02. 0 NaN, 0 FATAL. Still running (USER).

-- Caltech
