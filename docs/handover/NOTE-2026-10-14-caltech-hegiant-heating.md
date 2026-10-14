# NOTE Caltech -> Raven: He giant FIN3 base heating is ~400x the star's luminosity; chain stopped after link 1 (USER)

USER decision 10-14: tell Raven; stop the He giant chain after link 1. Links 2-6 (4454172-76) removed from the queue.
Link 1 (4454171) keeps running to its 12 h end (~t 8.5e5 s); all data (bins, rst, hst, logs) kept in
/resnick/groups/carnegie_poc/jingze/hegiant_1013/fin3. Link 1 is numerically clean so far (t 6.67e5 s, cycle ~35.3k, dt 18.18 s,
0.93 s/cycle, 0 FATAL / DIVERGED / NON-CONVERGED / NOT ADMISSIBLE, ~86 GB/H200).

## The problem
- he_base_heat keeps the frozen scaffold's wall-face F_MLT(r_in) as heating (code: F_MLT(r_in) = 4.70239e16 printed at start).
  P_esrc = 2.145e39 erg/s in the wedge (Omega = pi/3), constant = **6.7e6 Lsun full sphere**.
- MESA model (M2pt754_Porb100 profile29): L_ph = **1.62e4 Lsun**. The heating is ~400x L, ~70x a rough Eddington L for 2.75 Msun.
- Source: ic_giant_own.txt (md5 ddc72692) col 7 F_MLT/F = 0.995 at the wall. With col 4 F_r it implies
  L_tot = L_rad/(1-f): 8.0e6 Lsun at 2.0 Rsun, 6.7e6 at 3.1, 5.7e6 at 4.2 (not even constant in r).
  L_rad itself at 3 Rsun is 3.0e4 Lsun > L_MESA, so that layer should be near-radiative with ~0 convective flux;
  the 0.995 fraction looks like an IC-builder artifact (make_ic_giant_mesa.py).
- Wall radiative L_in = 9.83e36 wedge = 3.1e4 Lsun is also ~1.9x L_MESA.

## Run evidence (user.hst, wedge values)
| t [s] | L_top | Etot | KE_int | KEr_int | Mr_int | w_mlt |
|---|---|---|---|---|---|---|
| 0 | 3.62e36 | 1.605e47 | 0 | 0 | 0 | 1.00 |
| 2.5e5 | 1.72e36 | 1.610e47 | 2.15e43 | 2.09e43 | 2.18e37 | 1.00 |
| 5.0e5 | 2.14e36 | 1.615e47 | 1.38e44 | 1.15e44 | 4.51e37 | 0.94 |
| 6.67e5 | 2.17e36 | 1.619e47 | 1.72e44 | 1.31e44 | 4.92e37 | 0.43 |

- dEtot = +1.43e45 = P_esrc * t (energy budget closes: L_in - L_top tiny vs P_esrc). The heat stays in and drives outward
  motion (Mr_int > 0, 76 % of KE radial); e_gas -5 %, E_rad -4 % (expansion). By 30 d the heater would add ~6.7e46 erg
  (full sphere), ~15 % of a rough envelope binding energy (~5e47): the "expansion" would be heater-driven.
- Mass conserved (M_tot 1.5434578814e32 to 13 digits), Mdot_top 0.

## Suggested fix (Raven's call)
Base heating = L_MESA(r_in) - L_rad(r_in) (or 0 if negative) instead of the scaffold's F_MLT(r_in); likewise rebuild col 7
so F_MLT = max(0, F_MESA - F_r). Caltech will report link 1's end and wait for a new TASK before any He giant restart.
