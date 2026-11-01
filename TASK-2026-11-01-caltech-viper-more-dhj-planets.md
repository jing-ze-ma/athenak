# TASK 2026-11-01 Caltech -> viper (USER): setups for more dhj planets (Seidel 2026 targets): KELT-20b, TOI-1518b, WASP-127b

**User (11-01):**
- Caltech H200 nodes can likely hold 2-4 dhj runs per 4-GPU node. We are measuring the packing now.
- So we want a few more planets next to WASP-121b, taken from the planets in **Seidel et al. 2026**.
- Candidates named by the user: **TOI-1518b, KELT-20b (MASCARA-2b), WASP-127b**.

Caltech has none of these setups.
- KELT-20b exists only on viper: /viper/ptmp2/jinma/kelt20_0925, ready 09-25, never run.
- Caltech has no Seidel 2026 paper notes and no SED/ck files for these hosts.
- Caltech does have inputs/production/tools/planet_setup.py (fin51-r-fix).

Please:
1. **Seidel 2026:** give the planet list from that paper, with the observables per planet that the paper measures (limb or Doppler wind velocities, species, dKp, ...), so each run has a target to compare against. Say if any other planet there suits our code better than the three above.
2. **KELT-20b:** push the 09-25 package to bsg-files-1009, under docs/handover/dhj_planets_1101/kelt20b/ or similar:
   - deck(s);
   - ic_k20.txt;
   - sw_band_flux_K20_11.txt;
   - the planet_setup.py command line and outputs;
   - run.sub;
   - CPU gate results.
   Rebase the deck on fin51-r-fix 182477553 defaults: lhllc + cb_sensor, am_vector=true, problem/am_transport=false (hydro first).
   Open items from 09-25 needing a recommendation:
   - M_p, unmeasured: 2.0 M_J adopted, range 1.5-3.38;
   - albedo, disputed: 0.32 adopted, range 0.1-0.4;
   - 14 % of the stellar flux below 0.26 um is folded into band 10 (A-star host);
   - IC T is up to 2 % off the table above 1.7e-4 bar.
3. **TOI-1518b and WASP-127b:** the same package built from the literature, using the same pipeline as WASP-121b/KELT-20b: system parameters with references, host SED -> 11-band sw flux, metallicity, T_int, albedo, R_p at the computed transit pressure, IC, keys, grids. Notes:
   - TOI-1518 is a hot fast-rotating A star, so it has the same short-UV band issue as KELT-20.
   - WASP-127b is a much cooler, low-density sub-Saturn. Check that our EOS/ck tables and floors cover its T range, its large scale height and its low gravity. Flag anything that breaks: domain height, the sponge, or the dt the low g and the large H imply.
4. **Per planet:** estimate the cost per rotation relative to WASP-121b on the hydro fresh deck grid (C32 x nx1 112 class), and give the rotations needed (300 like WASP-121b?).
5. **No production on viper**, and only small CPU checks; the runs go on Caltech. Reply in a NOTE, small files only (give viper paths for anything > ~50 MB).

Roche/tide: Caltech is porting the Roche potential to the dhj pgen. Include a/R_* and the tidal ratio (tide/g at the top of the domain) for each planet, so we know where the tide matters.
