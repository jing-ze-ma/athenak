# NOTE 2026-10-18 caltech -> raven: longer agcE1r nofloor comparison (answer to TASK 8ba77933)

GPU test only; agcE1r not touched (read-only copy of agcar3d.00302.rst in nofloor_1018/R).
- Binary athena_gpu_he_star_m1_4775edf5_nofma (md5 b375ccf1...), 1 node x 2 H200, t_rst 148042909.09 -> t_rst + 5 DT
  (DT 7e5, tlim 1.515429e8), hydro_w + m1 bins every DT (outputs 01481..01485 = 1..5 DT). All three arms reached tlim.
- Jobs: off 4925623 (hpc-sm-01-17), aN 4925624 (hpc-sm-01-11), ref 4925625 (hpc-sm-01-11). rc 0, FATAL 0,
  summary NON-CONVERGED=0 in all arms.
- **Wall: off 2951 s (49 min), aN 4267 s (71 min, +45 %), ref 5000 s (83 min, +69 %).** off ran on a different node
  than aN/ref, so the off-vs-aN wall ratio carries some node noise.

## Verdict
- **The nofloor effect on the solution is not distinguishable from trajectory noise.** At every output time the
  off-vs-aN difference is the same size as off-vs-ref and aN-vs-ref (mean dens 2e-3..1e-2, eint ~1.5e-2, m1_e ~1e-2).
  ref (tol 1e-10) is just a third chaotic member, not a closer answer.
- **Growth with time:** the differences do not grow steadily. Mean dens off-ref 2.3e-3, 6.5e-3, 5.6e-2, 1.2e-2,
  9.1e-3 (1..5 DT); aN-ref 2.6e-3, 5.6e-3, 6.7e-2, 8.3e-3, 1.6e-2. There is a transient peak at 3 DT in both pairs (the
  band 380:480 mean reaches 0.3-0.36; off-aN stays at 3.7e-3 at that time, so off and aN agree with each other and both
  differ from ref there). After that they fall back to ~1e-2. They saturate at the noise level; they do not accumulate.
- **At 5 DT, which is closer to ref:** it is mixed. off is closer in dens (mean 9.1e-3 vs 1.6e-2; band 4.7e-2 vs
  8.3e-2). aN is marginally closer in eint (1.4e-2 vs 1.5e-2 mean, but 4.8e-2 vs 4.1e-2 in band), m1_e (6.7e-3 vs
  7.2e-3) and m1_f1 (9.0e-3 vs 9.6e-3). The closer arm alternates between times: dens mean favours off at 1, 3, 5 DT
  and aN at 2, 4 DT.
- Solver level, nofloor does what it should: tail >1e-2 is 17.6 % in off, 6.4 % in aN and 6.2 % in ref. But this does
  not show up as a solution difference over 5 DT, and it costs +45 % wall.
- Side effect: aN and ref have 87 and 97 DT_RETRY. Each is a per-pass "Picard NON-CONVERGED after 200 passes"
  (resid > implicit_resid_fatal, resid ~1e1, worst cells near i~460, inside the band), recovered by a dt cut (cut 1;
  cut 2 twice in aN, once in ref). off has **0** retries, because the floor accepts before 200 passes. Summary
  NON-CONVERGED=0 everywhere, so none of these is a stop condition.
- Recommendation from these numbers: keep nofloor out of production (no measurable accuracy gain, +45 % wall, and it
  adds the 200-pass retries).

## Tool caveat (cmpbin.py band)
picstop-1017 cmpbin.py ravels and concatenates the meshblocks when n_mbs > 1. The agcE1r bins have 4 MBs of
[64,64,536], so its `--ib` band slice was a **no-op**: the "band" columns of NOTE-2026-10-17 equal the full domain.
Below, the same metrics were computed with mb_data kept as [mb,k,j,i] and the band taken as i 380:480 (script nofloor_1018/ana/
cmp5.py, the same formulas as cmpbin.py). The full-domain columns match cmpbin.py (off-ref at 1 DT reproduces the
10-17 numbers: dens 7.1e+01 / 2.3e-03).

## stat.py
```
arm              rc    cyc  solves  pass   NC retry tail>1e-4 tail>1e-2  max_acc  checks   fails  exh(NC) ce_refus  wall_s t_end
off               0 211577    9847  5.21    0     0   18.62%   17.59%  1.6e+04       -       -        -      -/-    2951 1.515429e+08
aN                0 212157   11008  6.95    0    87   11.61%    6.37%  2.2e+04 0.000000e+00 0.000000e+00 0.000000e+00 6.228000e+03/0.000000e+00    4267 1.515429e+08
ref               0 212225   11146  9.73    0    97   11.12%    6.17%  8.7e+04       -       -        -      -/-    5000 1.515429e+08
```

## Retries / per-pass NC (grep of run.log)
```
off: DT_RETRY 0 recovered 0 per-pass NC(after 200) 0 summary NC>0 0
<rad_m1> implicit transport: solves=9.847000e+03 Picard iterations mean=5.213060e+00 max=2.400000e+01 NON-CONVERGED=0.000000e+00
aN: DT_RETRY 87 recovered 82 per-pass NC(after 200) 87 summary NC>0 0
<rad_m1> implicit transport: solves=1.100800e+04 Picard iterations mean=6.946948e+00 max=4.900000e+01 NON-CONVERGED=0.000000e+00
<rad_m1> picstop: implicit_final_resid_stop=0.000000e+00 (rel 0.000000e+00 x tol_eff) checks=0.000000e+00 failed+continued=0.000000e+00 failed-no-pass-left(NC)=0.000000e+00 solves accepted after a failed check=0.000000e+00 passes after a failed check=0.000000e+00 | implicit_conv_est_guard=0.000000e+00 refused=0.000000e+00 | implicit_conv_est_nofloor=1 refused=6.228000e+03 extra stops=2.000000e+00 | log worst-cell time 3.321323e+00 s
ref: DT_RETRY 97 recovered 95 per-pass NC(after 200) 97 summary NC>0 0
<rad_m1> implicit transport: solves=1.114600e+04 Picard iterations mean=9.726539e+00 max=1.580000e+02 NON-CONVERGED=0.000000e+00
```

## Per-time comparison (max / mean |a-b|/|b|; m1_f1 relative to c*m1_e of the second arm; band i 380:480)
```
== output 1481  t=1.487430e+08 (t-t_rst = 1.00 DT)
off-ref    dens max 7.1e+01 mean 2.3e-03 | band max 7.1e+01 mean 1.2e-02 ; eint max 2.0e+01 mean 1.2e-02 | band max 2.0e+01 mean 2.6e-02 ; m1_e max 4.5e+00 mean 7.6e-03 | band max 4.5e+00 mean 1.5e-02 ; m1_f1 max 9.5e-01 mean 5.4e-03 | band max 9.5e-01 mean 6.7e-03
aN-ref     dens max 2.9e+02 mean 2.6e-03 | band max 2.9e+02 mean 1.3e-02 ; eint max 1.1e+02 mean 2.0e-02 | band max 1.1e+02 mean 3.9e-02 ; m1_e max 5.0e+00 mean 1.3e-02 | band max 5.0e+00 mean 2.3e-02 ; m1_f1 max 1.1e+00 mean 9.1e-03 | band max 1.1e+00 mean 1.2e-02
off-aN     dens max 5.3e+01 mean 1.7e-03 | band max 5.3e+01 mean 8.4e-03 ; eint max 1.3e+01 mean 1.5e-02 | band max 1.3e+01 mean 2.8e-02 ; m1_e max 4.9e+00 mean 1.0e-02 | band max 4.9e+00 mean 1.7e-02 ; m1_f1 max 6.8e-01 mean 7.3e-03 | band max 6.8e-01 mean 8.1e-03
== output 1482  t=1.494430e+08 (t-t_rst = 2.00 DT)
off-ref    dens max 9.5e+02 mean 6.5e-03 | band max 9.5e+02 mean 3.4e-02 ; eint max 2.5e+02 mean 2.0e-02 | band max 2.5e+02 mean 4.5e-02 ; m1_e max 9.0e+00 mean 1.2e-02 | band max 9.0e+00 mean 2.2e-02 ; m1_f1 max 1.1e+00 mean 8.5e-03 | band max 1.1e+00 mean 1.1e-02
aN-ref     dens max 5.8e+02 mean 5.6e-03 | band max 5.8e+02 mean 2.9e-02 ; eint max 1.5e+02 mean 1.3e-02 | band max 1.5e+02 mean 3.1e-02 ; m1_e max 6.1e+00 mean 8.1e-03 | band max 6.1e+00 mean 1.3e-02 ; m1_f1 max 8.6e-01 mean 5.8e-03 | band max 8.6e-01 mean 5.2e-03
off-aN     dens max 1.1e+02 mean 3.0e-03 | band max 1.1e+02 mean 1.5e-02 ; eint max 3.3e+01 mean 2.0e-02 | band max 3.3e+01 mean 3.6e-02 ; m1_e max 1.0e+01 mean 1.3e-02 | band max 1.0e+01 mean 2.1e-02 ; m1_f1 max 1.3e+00 mean 9.9e-03 | band max 1.3e+00 mean 1.0e-02
== output 1483  t=1.501430e+08 (t-t_rst = 3.00 DT)
off-ref    dens max 6.8e+03 mean 5.6e-02 | band max 6.8e+03 mean 3.0e-01 ; eint max 7.0e+02 mean 2.7e-02 | band max 7.0e+02 mean 1.0e-01 ; m1_e max 1.0e+01 mean 9.1e-03 | band max 1.0e+01 mean 2.1e-02 ; m1_f1 max 1.2e+00 mean 5.9e-03 | band max 1.2e+00 mean 9.1e-03
aN-ref     dens max 8.6e+03 mean 6.7e-02 | band max 8.6e+03 mean 3.6e-01 ; eint max 9.2e+02 mean 2.7e-02 | band max 9.2e+02 mean 1.1e-01 ; m1_e max 1.0e+01 mean 8.6e-03 | band max 1.0e+01 mean 2.0e-02 ; m1_f1 max 8.3e-01 mean 5.5e-03 | band max 7.6e-01 mean 7.0e-03
off-aN     dens max 9.9e+01 mean 3.7e-03 | band max 9.9e+01 mean 1.9e-02 ; eint max 3.8e+01 mean 1.6e-02 | band max 3.8e+01 mean 3.7e-02 ; m1_e max 9.2e+00 mean 9.8e-03 | band max 9.2e+00 mean 2.1e-02 ; m1_f1 max 1.3e+00 mean 6.5e-03 | band max 1.3e+00 mean 8.8e-03
== output 1484  t=1.508430e+08 (t-t_rst = 4.00 DT)
off-ref    dens max 6.9e+02 mean 1.2e-02 | band max 6.9e+02 mean 6.4e-02 ; eint max 1.1e+02 mean 1.9e-02 | band max 1.1e+02 mean 5.2e-02 ; m1_e max 7.7e+00 mean 9.9e-03 | band max 7.7e+00 mean 2.0e-02 ; m1_f1 max 5.9e-01 mean 6.7e-03 | band max 5.9e-01 mean 7.9e-03
aN-ref     dens max 1.0e+03 mean 8.3e-03 | band max 1.0e+03 mean 4.3e-02 ; eint max 1.6e+02 mean 1.6e-02 | band max 1.6e+02 mean 4.3e-02 ; m1_e max 8.7e+00 mean 8.4e-03 | band max 8.7e+00 mean 1.8e-02 ; m1_f1 max 6.9e-01 mean 5.6e-03 | band max 6.9e-01 mean 7.1e-03
off-aN     dens max 9.4e+02 mean 1.2e-02 | band max 9.4e+02 mean 6.3e-02 ; eint max 1.6e+02 mean 1.5e-02 | band max 1.6e+02 mean 4.3e-02 ; m1_e max 7.8e+00 mean 7.1e-03 | band max 7.8e+00 mean 1.4e-02 ; m1_f1 max 8.5e-01 mean 4.7e-03 | band max 8.5e-01 mean 5.2e-03
== output 1485  t=1.515430e+08 (t-t_rst = 5.00 DT)
off-ref    dens max 3.9e+02 mean 9.1e-03 | band max 3.9e+02 mean 4.7e-02 ; eint max 6.2e+01 mean 1.5e-02 | band max 6.2e+01 mean 4.1e-02 ; m1_e max 7.0e+00 mean 7.2e-03 | band max 7.0e+00 mean 1.4e-02 ; m1_f1 max 9.8e-01 mean 9.6e-03 | band max 9.8e-01 mean 1.8e-02
aN-ref     dens max 1.8e+03 mean 1.6e-02 | band max 1.8e+03 mean 8.3e-02 ; eint max 2.6e+02 mean 1.4e-02 | band max 2.6e+02 mean 4.8e-02 ; m1_e max 9.8e+00 mean 6.7e-03 | band max 9.8e+00 mean 1.3e-02 ; m1_f1 max 9.8e-01 mean 9.0e-03 | band max 9.8e-01 mean 1.8e-02
off-aN     dens max 7.3e+02 mean 7.1e-03 | band max 7.3e+02 mean 3.6e-02 ; eint max 1.1e+02 mean 1.6e-02 | band max 1.1e+02 mean 3.8e-02 ; m1_e max 1.3e+01 mean 8.1e-03 | band max 1.3e+01 mean 1.5e-02 ; m1_f1 max 6.4e-01 mean 5.6e-03 | band max 6.4e-01 mean 5.5e-03
```
