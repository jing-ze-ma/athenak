# NOTE 2026-10-22 viper -> caltech (USER): the N^2 seed l window (16-32) may need to change

Re TASK a1fe4fcd (N^2 seed default) and fbe76481 (N^2 MLT ramp default): both are being implemented on viper (branch n2seed-1022,
off fin51-r-fix 1976bac68; NOTE with gates follows). Heads-up on the l window, USER decision 10-22:

- Your N^2 analysis puts the fastest growth at l ~ pi r/H_p ~ 200 (RSG E) and ~ 460 (HG). The seed decks use l 16-32 (code cap 32),
  so most seed power sits well below the fast-growing scales; the fast modes then come from grid noise / cascade, and a narrow
  window can imprint its own scale.
- Spectrum today (he_seed_ylm): Gaussian coefficient per (l, m), every m -> power per l ~ (2l+1), flat per coefficient; radial
  envelope one half-sine over the band; pattern fixed by he_seed_rng (default 1234 -- same pattern in every run unless changed).
- New default rule (he_seed_band = n2): l_hi = min(pi r/H_p at the band, l with ~8 cells per wavelength on the grid); l_lo = l_hi/4
  (two octaves); flat per coefficient; cap raised 32 -> 64 (clipped with a warning, not FATAL). Explicit he_seed_lmin/lmax override.
  For cs N=128 the grid limit is roughly l ~ 64, so expect ~16-64 there; numbers per deck come with the gate NOTE.
- So: the RSG E fresh run started 09:29 with l 16-32 is fine to keep (USER decides on any restart); for the HG fresh start, consider
  the new window (or wait for the branch). Also consider a different he_seed_rng per production run.
