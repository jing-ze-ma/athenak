# TASK 2026-10-22 caltech -> viper (USER): MLT ramp from the IC N^2 growth rates -- make it the DEFAULT (with the N^2 seed, a1fe4fcd)

USER: the MLT closure must switch off at a time set by the convection's growth rate (adaptive alone never switches off where the
resolved flux falls short, e.g. under-resolved surface layers). Rule used on Caltech (N^2 from the IC column with the code EOS +
radiation, docs/handover/n2seed_profile_1022.py; inside the simulated domain):
  mlt_ramp_start = 7 / sigma_max          (ln 10^3: 1e-3 seed -> saturation at the fastest-growing, surface, layer)
  mlt_ramp_time  = mass-weighted median of 1/sigma over the convective envelope (bulk of the envelope has one e-fold to take over)
  RSG E (g25): 1/sigma_max 7.35 d -> start 51.4 d (4.44e6 s); median deep 1/sigma 130.8 d -> ramp 1.130e7 s; MLT off at 182 d
               (was: no ramp, adaptive only).
  HG global:   1/sigma_max 0.59 d -> start 4.2 d (3.63e5 s); median 6.9 d -> ramp 5.96e5 s; MLT off at 11 d (was 14.9 d + 74 d).
viper: make this the he_star_m1 default for fresh starts: e.g. mlt_ramp_start / mlt_ramp_time = auto (computed once from the IC column
with the code EOS, printed with sigma_max, the band and the median; restart state = the two numbers), explicit numbers override,
mlt_ramp_start = -1 keeps "no ramp". Same gates as the seed default (bitwise when named explicitly / off).
Status: RSG E restarted FRESH on Caltech 10-22 09:29 with the N^2 seed (l 16-32, 0.97-0.999 R_ph; log line confirmed) + fix keys
(harmonic, extrap2, accept_fallback, posfix=4) on ff1f52bc; 1 node 4 H200, ~2.5 s/cycle at dt 1000 s, 0 retries at 1.7 d; the ramp
applies from its first restart (~5 d). HG global restarts FRESH with the same rules when its allocation starts.
Side finding: kinetic-energy flux is not in the adaptive closure's resolved flux; for RSG E it is negligible (|L_KE/L| <= 0.002 at 129 d).
