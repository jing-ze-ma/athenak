# NOTE Caltech: plm-arm G1/G3 gates rerun after the plm halo fix (e50ec7d6)

picked up by Caltech 2026-10-10 05:5x PDT; done 2026-10-10 06:07 PDT.
Answers TASK-2026-10-10-caltech-plm-gates-rerun.md.

## Binary
- fork xthinfix-1009 e50ec7d6, CPU, PROBLEM=built_in_pgens, Release, gcc 13.2 + openmpi 5.0.1 (incremental build_inc.sh, 39 files recompiled; job 4316035)
- /resnick/home/jingze/ATHENAK/builds/bin/athena_cpu_built_in_pgens_e50ec7d6  md5 f684f0c07d411e10a347062a4192a01b

## Runs
- job 4316036 (expansion, 12 cores, 10:44); scripts docs/handover/xthinfix-1009/arms/ (Caltech copy battery_arms_caltech.sh,
  restricted to arms sq aq knp and g1/g3 only; same fixkeys.py vet_tensor workaround as the arm matrix)
- run dir /resnick/groups/carnegie_poc/jingze/arms_1010/run_plm/gates (gates.out = raw lines)

## Results (raw)
| run | rc | max|dE/E| | top5 | L1 | Picard mean/max | solves | NON-CONVERGED |
|---|---|---|---|---|---|---|---|
| g1_sq_1e2  | 0 | 1.738e-02 | 4.567e-03 | 6.514e-04 | 2.06 / 18 | 5119 | 0 |
| g1_sq_1e4  | 0 | 1.632e-02 | 4.491e-03 | 7.236e-04 | 22.0 / 30 | 51 | 33 |
| g1_aq_1e2  | 0 | 1.749e-02 | 4.910e-03 | 1.210e-03 | 2.06 / 18 | 5119 | 0 |
| g1_aq_1e4  | 0 | 1.774e-02 | 5.407e-03 | 1.943e-03 | 21.5 / 30 | 51 | 32 |
| g1_knp_1e2 | 0 | 1.744e-02 | 4.872e-03 | 7.816e-04 | 2.06 / 17 | 5119 | 0 |
| g1_knp_1e4 | 0 | 1.766e-02 | 5.342e-03 | 1.238e-03 | 21.5 / 30 | 51 | 32 |
| g3_sq_c1   | 0 | - | - | 0.0136 PASS | 2.01 / 6 | 51711 | 0 |
| g3_sq_c10  | 0 | - | - | 0.0138 PASS | 2.53 / 10 | 5171 | 0 |
| g3_aq_c1   | 0 | - | - | 0.0137 PASS | 2.01 / 6 | 51711 | 0 |
| g3_aq_c10  | 0 | - | - | 0.0140 PASS | 2.53 / 10 | 5171 | 0 |
| g3_knp_c1  | 0 | - | - | 0.0136 PASS | 2.01 / 6 | 51711 | 0 |
| g3_knp_c10 | 0 | - | - | 0.0138 PASS | 2.53 / 10 | 5171 | 0 |
(G3 L1 = L1(E) = L1(material), gate <= 0.02. No FATAL in any run.log.)

## Verdict
- Segfault fixed: all 12 plm runs rc 0 (were rc 139 at 413b38af).
- viper sq check reproduced: G1 cfl 1e2 max|dE/E| 1.738e-2 (viper 1.74e-2), L1 6.514e-4 (viper 6.5e-4); G3 cfl 10 L1 0.0138 (viper 0.0138).
- G1 cfl 1e4: Picard hits the 30 cap in 32-33 of 51 solves for every plm arm (NON-CONVERGED > 0), but the
  Hopf error stays at the cfl 1e2 level. Viper should decide whether that counts.
