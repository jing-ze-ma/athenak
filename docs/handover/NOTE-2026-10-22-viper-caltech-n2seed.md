# NOTE 2026-10-22 viper -> caltech (USER): N^2 seed + auto MLT ramp as he_star_m1 defaults (re TASKs a1fe4fcd, fbe76481)

Branch fork/n2seed-1022 = cf26bc2f7 (on fin51-r-fix 1976bac68; NOT fast-forwarded). CPU gates PASS; GPU gates + 10-day N=32
onset test (seed off/old/n2, scored by corr with dOmega) still queued on viper -- a follow-up NOTE will close it.
Keys (<problem>):
- he_seed_band = n2 (default) | old (= today's seed, bitwise; he_seed default 0). n2: N^2 = g (dlnP/dr/Gamma1 - dlnrho/dr) from the
  code's own IC cells (code EOS, P = p_gas + aT^4/3, Gamma1 from the first law incl. radiation, g = G m(r)/r^2), active cells r <= r_int.
  Band = cells with sigma >= he_seed_n2_frac (0.3) sigma_max, gaps <= he_seed_n2_gap (2) closed, runs < he_seed_n2_nmin (4) cells dropped
  (IC ripples), the run with the largest sigma kept; half-sine envelope across it. he_seed default 1e-3 under n2.
  l window (USER rule): l_hi = min(pi r/H_p, l_8), l_8 = 2 pi/(he_seed_n2_ppw (8) * smallest neighbour angle); l_lo = max(1, l_hi/4);
  cap 64 (clipped with WARNING); flat Gaussian coefficient per (l, m). Explicit he_seed_lmin/lmax in the deck OVERRIDE the window
  (prod/det decks name 1-8, your n2seed decks 16-32: delete those lines to get the automatic window). he_seed_rlo/rhi ignored under n2.
  Needs one MeshBlock across x1. Fresh starts only (restarts/remaps never reseed).
- mlt_ramp_start / mlt_ramp_time = auto (default): start = 7/sigma_max, time = mass-weighted median of 1/sigma over the unstable cells
  below r_int; computed once at the fresh start (with mlt_flux_frozen), stored in the restart (trailer HSRAMP01), reused, never
  recomputed. Explicit numbers as before; -1 = no ramp. A restart with auto and no stored numbers gets NO ramp (running productions
  cannot switch a ramp on).
- Keys not in the deck cannot be named on the command line of this build (FATAL): put them in the deck.
Startup values vs yours:
                  viper automatic                                         Caltech
  RSG E band      0.9756-0.9977 r_int, 42 cells; peak 0.9951, 1/sigma 6.11 d   0.97-0.999 R_ph, 53 cells; 0.993-0.996, 7.35 d
  RSG E ramp      start 3.70e6 s (42.8 d), time 9.29e6 s (107.6 d)            4.44e6 s, 1.130e7 s
  HG band         0.9767-0.9987 r_int, 24 cells; peak 0.9959, 1/sigma 0.775 d  0.98-0.999, 21 cells; 0.59 d
  HG ramp         start 4.69e5 s, time 6.01e5 s                               3.63e5 s, 5.96e5 s
  Bands agree to ~0.005 R_ph; sigma_max differs +20 % (RSG) / -24 % (HG) -- likely centred differences on code cells vs np.gradient on
  MESA points (not tested). RSG ripple run at 0.066 (1 cell) dropped, as you excluded it.
  l window: cs N=16 3-10, N=32 6-22, N=128 (prod) 16-64 (pi r/H_p 341 / 453 > cap, clipped).
CPU gates (cs N=16, 6 ranks, old 1976bac68 vs new, same platform): band=old (+ ramp -1 / explicit) IDENT for RSG E and HG; restart of a
1976bac6 rst IDENT (deck keys, and ramp auto -> no ramp); n2 + auto fresh C2==C IDENT, R==C IDENT (ramp read back from the rst).

## Addendum (GPU gates + onset test) -- all gates pass
GPU (viper apu, cf26bc2f7): old seed bitwise vs 1976bac6 for RSG E and HG; rst restart identical with deck keys and with ramp auto;
n2 + auto C2==C and R==C identical; smokes clean; GPU startup band/ramp = CPU. At N=32 the l window is 6-22 (RSG E and HG).
RSG E N=32, 10 d (seed off / old / n2): surface convection starts by 1.2 d at 0.99-0.999 R_ph in ALL arms, identical for off and n2
(IC/grid-driven at N=32); n2 adds ~1e2 cm/s at 0.98 R_ph (inside its band). The old seed's deep large-scale flow (3.3e3 cm/s at
0.5 R_ph by 10 d, corr(v_r, dOmega) +0.04 constant = the fixed seed pattern, reaching down to 0.05 R_ph) is ABSENT with n2:
0.24 cm/s at 0.5 R_ph, same as seed off. From day 7 both n2 and off show corr +0.04..0.06 at 0.99 R_ph = geometric surface part,
not the seed.
