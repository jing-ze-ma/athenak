# TASK 2026-10-22 caltech -> viper + deltaai (USER): seed convection where the IC N^2 says it grows fastest; make it the DEFAULT

Finding (docs/handover/n2seed_profile_1022.py: N^2 = g (dlnP/dr / Gamma1 - dlnrho/dr) from the IC column with the code's general EOS
(eos_gen port) + radiation P = aT^4/3, Gamma1 from the first law; only inside the simulated domain):
  RSG E (g25): fastest growth at 0.993-0.996 R_ph (T ~6.6 kK, H recombination), 1/sigma = 7.4 d, l ~ pi r/H_p ~ 200.
               Three thin spikes at 0.056-0.065 R_ph (1/sigma 10-11 d, l ~ 9) look like IC-column stitching ripples (not seeded).
  HG global:   fastest at 0.992-0.996 R_ph (T ~14 kK), 1/sigma = 0.59 d, l ~ 460; and 0.987-0.989 (24 kK), 1.1 d, l ~ 140.
  The deep envelope grows 5-20x slower (RSG 0.5 R_ph: 160 d; HG 0.5 R_ph: 9 d). The old seed (Y_lm l 1-8, 0.05-0.99 / 0.07-0.94 R_ph)
  put large-scale power deep, where it imprints a cube-locked pattern; the RSG run itself showed onset at 0.99 R_ph by 6 d.
New seed (decks rsg_E_1022_n2seed / hgglob_1022_n2seed): he_seed 1e-3, he_seed_ylm l 16-32 (code cap 32), radial band RSG 0.97-0.999
R_ph (53 cells), HG 0.98-0.999 R_ph (21 cells); everything else = production. Caltech starts RSG E FRESH with it now.
viper (USER "make that seed default"): make the N^2-placed surface seed the he_star_m1 DEFAULT for fresh starts -- e.g. he_seed_band = n2
(band = cells inside the domain with sigma >= 0.3 sigma_max from the IC column + the code EOS; exclude isolated IC ripples), l range
from pi r/H_p at the band clipped to the grid (raise the he_seed_lmax cap 32 -> 64, cost (lmax+1)^2 terms per cell); key off = old seed.
deltaai: for TASK c2c8bf6e (HG global FRESH) use docs/handover/hgglob-deltaai-1022/hgglob_1022_n2seed.athinput instead of
hgglob_fresh_clip_harm.athinput (same @FILES@ substitution, same command-line keys). If already started with the old seed, keep it
running and say so in your NOTE (USER decides).
