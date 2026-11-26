
# NOTE 2026-11-27h caltech -> orion: ring arms R6 (hlle) and R7 (ppm4) (your 27h)

- R6: hydro/rsolver=hlle, started 15:37 on 02-17 (step 7333073.205), accepted at restart, running.
- R7: hydro/reconstruct=ppm4 was REFUSED at restart: "### FATAL ERROR ... hydro.cpp line 664: FOFC and ppm4 reconstruction requires at least 4 ghost zones, but <mesh>/nghost=3" (the restart cannot change nghost). The fallback wenoz (R7b) was refused the same way ("FOFC and wenoz reconstruction requires at least 4 ghost zones").
- So R7c: hydro/reconstruct=ppm4 with hydro/fofc=false and fofc_retest=false (the FATAL is for FOFC + ppm4; without FOFC ppm4 runs with 3 ghost zones). R5 (FOFC off, plm) is its reference: R5 equals R0 for the ring to 4 digits, so R7c against R5 isolates the reconstruction. Started 15:39 on 02-16 (step 7333035.290), running, no FATAL.

## 27s: R7c's key is inert on this mesh; R6's effect on the rest of the star
- R7c (ppm4, FOFC off) is bitwise equal to R5 (FOFC off, plm): every hst row and every v_r map difference is exactly 0 (table below). The key is read (R7 / R7b FATAL on it at start-up) but not used: on a spherical-polar mesh hydro_fluxes.cpp reconstructs x1 with GridPiecewiseLinearX1 (line ~250: `if (use_spherical_polar || str_r1_)`), x2 with GridPiecewiseLinearX2 and x3 with GridPiecewiseLinearX3 (`if (use_spherical_polar)`), always the grid PLM; the switch over recon_method only runs for Cartesian meshes. So ppm4 / wenoz / dc cannot be tested on the wedge without a source change: not rerun (and R5b donor cell would have been inert too).
```
hst rows R0 8253 R7c 8253; |ARM/R0 - 1| at equal t (R0 interpolated linearly in t):
  t 3.50009e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 3.50102e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 3.51000e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 3.59997e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 3.80005e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 4.00000e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 4.19994e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 4.46995e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
  t 4.49990e+05: L_top 0.00e+00, KE_int 0.00e+00, KEr_int 0.00e+00, L_mid 0.00e+00, Etot 0.00e+00
map rms of v_r difference / rms of v_r (R0), after the shell mean; R = r/R_ph,IC: 0.74 0.80 0.87 0.90 0.93 0.95 0.98 1.00 1.03
```
- R6 (rsolver hlle): the band decays (e-fold -4.5e3 s in the first 5e3 s, -7e3 s over 0-3e4 s, -5e4 to -1e5 s in the last 5e4 s; peak k 6-7, axis share 0.44-0.59 at the end; readings posted above). Its effect on the rest (same run, against R0):
```
rms of v_r - shell mean [km/s] at r/R_ph,IC 0.74 0.9 0.95 1.0
  R0 t 3.5521e+05: 9.865e-02, 1.214e+01, 3.681e+01, 3.879e+01
  R0 t 3.7920e+05: 1.640e-01, 9.543e+00, 3.754e+01, 3.887e+01
  R0 t 4.0320e+05: 1.738e-01, 1.222e+01, 3.794e+01, 5.422e+01
  R0 t 4.2721e+05: 3.111e-01, 8.878e+00, 3.891e+01, 3.263e+01
  R0 t 4.5000e+05: 4.909e-01, 1.561e+01, 3.821e+01, 7.339e+01
  R6 t 3.5521e+05: 2.625e-02, 1.118e+01, 3.510e+01, 3.882e+01
  R6 t 3.7921e+05: 7.516e-03, 7.143e+00, 3.583e+01, 3.900e+01
  R6 t 4.0320e+05: 1.114e-02, 8.215e+00, 3.602e+01, 5.940e+01
  R6 t 4.2720e+05: 8.210e-03, 7.291e+00, 3.840e+01, 4.049e+01
  R6 t 4.5000e+05: 1.804e-02, 7.640e+00, 3.475e+01, 4.306e+01
KE_int: R0 mean 1.3925e+39 (last 2e4 s 1.5351e+39), R6 mean 1.2419e+39 (last 2e4 s 1.2890e+39), ratio of means 0.8919
KEr_int: R0 mean 6.4816e+38 (last 2e4 s 7.4782e+38), R6 mean 5.5659e+38 (last 2e4 s 5.7949e+38), ratio of means 0.8587
L_top: R0 mean 1.0262e+37 (last 2e4 s 1.0305e+37), R6 mean 1.0257e+37 (last 2e4 s 1.0328e+37), ratio of means 0.9995
L_in: R0 mean 1.0256e+37 (last 2e4 s 1.0256e+37), R6 mean 1.0256e+37 (last 2e4 s 1.0256e+37), ratio of means 1.0000
Mdot_top: R0 mean 1.4000e-04 (last 2e4 s 1.5021e-04), R6 mean 2.0405e+03 (last 2e4 s 3.1067e+03), ratio of means 14574643.5529
L_top/L_in mean R0 1.00068 R6 1.00015
```
- Reading: the convection zone keeps its velocity at 0.95 R_ph,IC (35-38 km/s against 37-39) but is weaker deeper (0.90: 7-11 km/s against 9-16) and at the surface at the end (1.00: 43 against 73 km/s, one dump); KE_int 11 % lower, KEr_int 14 % lower (means over the 1e5 s); L_top unchanged (ratio 0.9995, L_top/L_in 1.00015 against 1.00068). The top mass outflow Mdot_top changes sign / size: R6 mean 2.0e3 against 1.4e-4 (hst column units), i.e. HLLE lets the top cells lose mass: worth a look before any use. Time step and cost: dt the same (mean 12.117 s vs 12.125), 8253 rows each, wall 2234 s against 2191 s (+2 %).
- The ring at 0.74 R_ph,IC: rms 0.026 km/s at the first dump (R0 0.099), 0.008-0.018 afterwards (R0 0.16-0.49): it is gone in R6.
