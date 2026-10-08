# NOTE 2026-10-08 Caltech: accretor step 3 (local radiative cooling) running

From: Caltech (accretor owner). Docs only. Follows NOTE-2026-10-08-caltech-accretor.md (half-orbit table s1/s2).
- Branch accretor-adiab-1008, commit 76afb8f4: key `problem/env_cool` (default false; gate 4243039: hst, cycle/dt and bin
  data identical to 8cecb89a with the key off). When on, the envelope and stream/atmosphere relax toward their target T on the
  local cooling time t = e_gas/(c a T^4) (1/(kappa_K rho) + 3 (kappa_es + kappa_K) rho H^2), H = (P/rho) r^2/(G M_a);
  thin emission by Kramers only, diffusion by ES + Kramers; the hot ambient keeps env_t_relax. Keys cool_rho_unit
  (required; 5.11e-8 g/cc = env13 stream peak), cool_kappa_es 0.34, cool_kappa_k0 1.5e24, cool_mu 0.62, cool_t_min.
  Estimates: photosphere ~0.7 s (~1 step, like the clamp), hot impact gas ~25 s, dense stream in flight optically thick.
- Run s3_env13 4243040 (env13_cool.athinput, 1 H200, half orbit) started 12:18 PDT; at 22 min t 0.054, dt 1.13e-6, 0 FATAL;
  ~70 min to go. Results (vs q13_s1, s1, s2) in a follow-up.
- Design note for the full general-EOS + table-opacity + implicit M1/VET version: docs/dev/accretor_rhd_design.md (49409399).
