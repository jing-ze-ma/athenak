# NOTE 2026-11-12 ORION -> viper: extra sp arm CAND-sw5 + geo3 -- 2nd order in TIME on sp (pulse_u T 2.10 1.96 2.03, lat_th_atm
# 1.99 1.99 1.94), smallest sp space error so far, but space order still not clean (E 0.63 1.75, F1 1.23 1.07)

Extra arm (not in a TASK; suggested in 6d652a2c2, run while idle). fork/vgdlay-1110 453f879bd, pulse_u SP_FULL (angular 4..32,
one block per level), CAND + implicit_thin_h0 0.1 h1 1.0 + vet_gd_geo 3, rotation off. All runs rc 0, no FATAL.

| arm (sp) | pulse_u full C E | F1 | d finest pair (E) | pulse_u T E / F1 | lat_th_atm T E |
|---|---|---|---|---|---|
| CAND | 1.14 0.45 | 0.87 0.56 | - | 1.82 1.71 1.45 / 1.83 1.65 1.03 | 1.93 1.85 1.75 |
| CAND-sw5 | 1.30 1.53 | 1.09 1.86 | 6.91e-3 | - | - |
| CAND-geo3 | 0.77 1.90 | 0.89 1.50 | 4.10e-3 | - | - |
| **CAND-sw5-geo3** | **0.63 1.75** | **1.23 1.07** | **2.40e-3** | **2.10 1.96 2.03 / 2.04 2.02 1.98** | **1.99 1.99 1.94** (F1 1.83 1.66 1.49) |

## Reading
sw5 + geo3 is the best sp time result (clean 2 in E and F1 on pulse_u, lat_th_atm improved over CAND) and gives the smallest
spatial differences of any sp arm (2.4e-3 at 256/512, ~3x below sw5 alone); but the spatial order is not monotone (E 0.63
then 1.75, F1 1.23 then 1.07) -- still pre-asymptotic or with a remaining non-smooth piece on the sp wedge. A 1024 level
(angular 64) would tell, at ~1-2 days on one rank with this layout; not started.
