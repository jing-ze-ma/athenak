# TASK 2026-10-22 caltech -> viper (USER): test <rad_m1>/implicit_face_opac = harmonic on the ACCRETOR (and BSG if cheap)

## Status
- Seed-free cs HSE test clean (he_seed = 0, HG deck N=32 to 0.3 d): lateral sd of rho, e/rho, E at round-off (<= 4e-16),
  v_lat ~1e-8 cm/s, R^2(dOmega) ~0 -> no well-balance fix needed; the interior pattern was the seed (your 57f45227 agrees).
- USER 10-22: both globals restarted FRESH at Caltech (21:12 PT) with implicit_face_opac = harmonic (binary d3df1401 = f55be371 +
  face_opac; decks hgglob_fresh_clip_harm / rsg_E_1020_clip_harm), inside the held allocations.
- USER 10-22: merge the face-opacity option into production as OPT-IN now: branch faceopac-prod-1022 = fin51-r-fix c92d5b03 +
  de0ee02b (kirchhoff) + d33f5f59 (harmonic); Caltech is building + gating it (key-off bitwise vs c92d5b03, key-on C/C2/R) and will
  fast-forward fin51-r-fix if clean (NOTE to follow). fsub / opac_lag all / opac_relax are NOT merged (diagnostic dead ends).
- Candidate DEFAULT: harmonic for FRESH starts only (restarts keep arith unless named), after it passes on the accretor + BSG.

## Asked of viper
1. Accretor (current prod deck v12lzc / your production setup): short A/B, implicit_face_opac = arith vs harmonic, from a restart
   of the running state (and a fresh smoke): dt / retries / redos, budget (energy/momentum booking you track), L and the thin-region
   fields, s/cycle. Harmonic leans the x1 face rho kappa toward the more transparent cell; check thin / optically-thin-to-thick
   transitions near the accretor surface and the stream for any degradation.
2. If cheap: BSG deck, same A/B.
3. Verdict: GO / NO-GO for harmonic as the fresh-start default.
