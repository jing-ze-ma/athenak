# NOTE Caltech -> Raven (10-14 21:30 PDT): AG Car E1 FLAG: d105 keeps jumping; vr105 approaching -5 km/s

Shell means from the bin dumps (pagerun = agcE1s to 2.6e7, ROL2U agcE1r after):
| t (s) | vr105 (km/s) | d105 (rho/rho0 at 1.05 R_ph) | rsurf |
|---|---|---|---|
| 2.42e7 | -2.7 | 96 | 1.042 |
| 2.55e7 | -2.7 | 134 | 1.043 |
| 2.59e7 (agcE1s) | -2.3 | 518 | 1.043 |
| 2.64e7 (ROL2U) | -3.4 | 515 | 1.045 |
| 2.65e7 | -3.0 | 1089 | 1.045 |
| 2.67e7 | -4.5 | 851 | 1.045 |
| 2.69e7 | -3.8 | 986 | 1.046 |
- The d105 jump to ~500 happened in agcE1s (FIN4-dtretry) before the ROL2U switch, so it is the evolving state, not ROL2U.
- vr105 trends to -4..-4.5 km/s (flag level -5 not yet reached); rsurf 1.046 (flag 1.15 far). KEt 5.6e41.
- agcE1r: 0.73 s/cycle, 0 DT_RETRY, 0 NON-CONVERGED, Picard ~5. Will flag at once at vr105 < -5.
