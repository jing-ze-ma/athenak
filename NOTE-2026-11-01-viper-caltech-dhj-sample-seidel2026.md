# NOTE 2026-11-01 viper -> Caltech (USER): revised dhj planet sample = 5 planets along the Seidel 2026 wind-Teq trend

Partial reply to TASK-2026-11-01-caltech-viper-more-dhj-planets (bed9c542e), items 1 and the planet choice.
Packages (item 2-4) follow in a later NOTE under docs/handover/dhj_planets_1101/.

**User decision (11-01):** the goal is to REPRODUCE the data points of Seidel et al. 2026. So we pick planets spread
along the paper's trend whose parameters are well constrained. **Sample: WASP-76b, WASP-121b, WASP-178b,
TOI-1518b, WASP-189b.** KELT-20b and WASP-127b are dropped (reasons below).

## The paper
Seidel, Parmentier, Prinoth et al. (32 authors, incl. K. Batygin), "Magnetic field strengths of hot giant exoplanets
consistent with Solar System values", Nature Astronomy, online 2 June 2026, DOI 10.1038/s41550-026-02870-1 (no arXiv
version found).
- **Measurement.** One number per planet: the net blueshift of Fe I at the terminator ("day-to-night wind"). It comes from ESPRESSO /
  MAROON-X transits, one uniform reduction, a 3000 K Fe template on a Kp-vsys map. The signal comes from ~1-100 mbar.
  The paper does not report limb asymmetry, dKp or phase-resolved data; those are in the earlier single-planet papers.
- **Result.** The wind FALLS with Teq. The authors explain it with Ohmic drag and infer **B <~ 2 G** (up to ~6 G for the two hottest planets).
- **Wind values.** The paper prints no table; the values are read from the vector paths of Fig. 1.
- Full summary, per-planet table and archive parameters with errors: docs/handover/dhj_planets_1101/seidel2026.md (this push).

| planet | Teq (paper) | wind [km/s] | M_p rel. err (source) | status |
|---|---|---|---|---|
| WASP-76b | 2229 K | 6.95 +- 0.36 | 1.5 % (Ehrenreich+2020) | **in** (low end, smallest error; puffy) |
| KELT-20b | 2263 K | 5.52 +- 1.30 | upper limit only (< 3.38 M_J) | dropped (mass; albedo disputed) |
| WASP-121b | 2358 K | 4.75 +- 0.30 | 3.7 % (Sing+2024) | **in** |
| WASP-178b | 2398 K | 4.79 +- 1.19 | 7.2 % (Hellier+2019) | **in** (radius err at the 5 % limit; A-star host) |
| TOI-1518b | 2491 K | 3.82 +- 0.50 | 25.7 % (Simonnin+2025) | **in**, marginal: the only well-measured wind at 2400-2640 K; nominal mass, gravity +-26 % as a systematic |
| HAT-P-70b | 2559 K | 2.55 +- 0.74 | upper limit only | dropped |
| WASP-189b | 2638 K | 1.05 +- 1.59 | 7.5 % (Lendl+2020) | **in** (hot end; table coverage of its dayside checked first) |
| WASP-127b | ~1400 K | not in the paper | 11.7 % | dropped (off the relation) |

- **Selection filter:** M_p measured to <= 15 %; R_p, Teq, a/R_* and Teff to <= 5 %; e consistent with 0. TOI-1518b is kept as the one exception.
- **Gap:** 2491-2638 K holds only HAT-P-70b, which fails the filter.
- **KELT-20b:** the 09-25 setup is packaged only as an optional mass-bracket pair (1 and 3.4 M_J), not as a trend point.

## What this implies for the runs
1. **Hydro alone should give winds that RISE with Teq.** The falling trend needs MHD or a drag term, which is the paper's own
   interpretation. Hydro first is the baseline; the comparison with the paper needs the MHD arms.
2. **The paper's B <~ 2 G is below our 10 G WASP-121b seed** and below the 10-100 G deep-dynamo scaling in
   NOTE-2026-11-01-viper-caltech-hj-dynamo-bfield-lit. A B scan of ~1-10 G across the sample would test the paper's
   interpretation directly.
3. **Observable from our output without line RT:** the terminator line-of-sight velocity weighted over 1-100 mbar.
   The Fe I cross-correlation velocity itself needs line RT, which we lack; an approximate Fe-opacity weighting is suggested in seidel2026.md.
4. **A-star hosts** (WASP-178, TOI-1518, WASP-189) have a large short-UV fraction folded into band 10 (KELT-20: 14 %).
   The packages will quantify it per host.

## Status on viper (no production here)
- Setup workers are building the packages:
  - TOI-1518b;
  - WASP-76b + WASP-178b;
  - WASP-189b, with the EOS/ck table-coverage check of its hot dayside first.
- KELT-20b is wrapped up as optional.
- Each package will give: parameters with references, an 11-band sw flux from the host SED, T_int, albedo, R_p at the transit
  pressure, IC, the deck on the hfresh defaults (lhllc, am_vector=true, problem/am_transport=false), grid, tide/g at top and bottom, a/R_*,
  cost per rotation vs WASP-121b, rotations needed, and a dynamo B estimate.
- Separately, the divergence-free MHD remap (C32 -> C64/C128) is being merged and gated on fin51-r-fix. A NOTE will follow when it is pushed.
