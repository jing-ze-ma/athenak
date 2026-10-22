# TASK 2026-10-22 caltech -> viper (USER): MLT ramp DEFAULT rule v2 -- start at t = 0, end at the convective saturation time (supersedes fbe76481)

Measured on the fresh RSG E run (N^2 seed l 16-32, 0.97-0.999 R_ph, amplitude 1e-3; MLT adaptive at full strength until 43 d):
  rms v_lat [km/s]   r/Rph 0.98   0.99   0.995  1.00
    30 d             0.003  0.005  0.008  0.004
    40 d             0.005  0.010  0.022  0.011
    44 d             0.006  0.013  0.033  0.017
    48 d             0.009  0.020  0.048  0.025      (v_r there 0.05-0.17 = radial pulsation)
  -> exponential growth with e-fold ~10 d at 0.995 R_ph = the N^2 peak (1/sigma_max 7.35 d; sigma_eff ~0.74 sigma_max), ALREADY at
     full MLT strength: the closure does NOT suppress the instability. Back-extrapolated v_lat(0) ~4e-4 km/s; saturation at v_MLT ~15 km/s
     needs ~10.5 e-folds -> t_sat ~105 d.
New default rule (USER): mlt_ramp_start = 0, mlt_ramp_time = t_sat = ln(v_MLT / v_seed) / sigma_eff, with sigma_eff ~0.74 sigma_max
(or simply ~10.5 / sigma_eff until the seed amplitude -> v_seed mapping is coded); i.e. MLT is gone when the resolved convection saturates.
  RSG E: t_sat ~105 d (Caltech run: ramp 43.0 -> 94.4 d, USER-set mid-run; kept).   HG global: 1/sigma_eff ~0.8 d -> ramp 0 -> 8.4 d
  (Caltech HG fresh deck hgglob_1022_n2seed now mlt_ramp_start 0, mlt_ramp_time 7.3e5 s; not started yet).
Please fold into the n2seed-1022 defaults (auto ramp from the IC N^2; explicit numbers override; -1 = no ramp).
