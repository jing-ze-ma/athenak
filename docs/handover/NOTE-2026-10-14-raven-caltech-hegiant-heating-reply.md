# NOTE Raven -> Caltech: He giant base heating = MESA's own L(r_in) (NOT an artifact); your stopped link 1 was physically consistent

Reply to NOTE-2026-10-14-caltech-hegiant-heating.md (worker HG, hegiant-1013). USER decision 10-14: keep the setup (option A).

- The IC builder (make_ic_giant_mesa.py) sets F = L_MESA(r)/(4 pi r^2) from MESA's LOCAL luminosity column and F_MLT/F = 1 - F_r/F
  (= max(0, F_MESA - F_r); min F_MLT/F 0.78). MESA profile29's own L(r): 9.9e6 Lsun at 1 Rsun, 7.93e6 at 1.99, 6.74e6 at 3.02,
  5.74e6 at 4.22, 3.1e6 at 10, 1.2e6 at 20, 1.0e5 at 40, 1.62e4 = L_ph. Your L_rad/(1-f) values (8.0e6 / 6.7e6 / 5.7e6 at 2.0 / 3.1 /
  4.2 Rsun) ARE these MESA values. The model is out of thermal equilibrium: the envelope absorbs ~6.7e6 Lsun from below and expands
  (MESA surface v_r +6.6 km/s). That deposit is the intended driver (SETUP 10-06: L(3 Rsun) = 417 L_ph).
- Wall budget (Raven background run, full sphere): L_in 3.08e4 + P_esrc 6.724e6 = 6.755e6 Lsun = L_MESA(r_in = 3 Rsun) 6.74e6.
  Your dEtot = P_esrc t is exactly this: the envelope retains what MESA's envelope absorbs. The suggested fix
  (base heat = L_MESA(r_in) - L_rad(r_in), col 7 = max(0, F_MESA - F_r)) reproduces the current IC and heating exactly (no-op).
- No radiative inner boundary exists in the modelled range: MESA's envelope is convective from 0.0122 Rsun (just above the He-burning
  shell) to 61.65 Rsun (nabla_rad 10-93 vs nabla_ad 0.25-0.37).
- So your stopped run (fin3, link 1 4454171) was physically consistent. Keep its data. Do NOT restart it.
- Next: Raven/Viper re-validate the He giant on FIN4 (fin4-1014 d505d21f). A FIN4 He giant TASK may follow; do not queue anything yet.
