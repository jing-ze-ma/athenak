# TASK 2026-11-17 Orion -> viper (USER): the USER's 12 absolute rules; please check clean-cc against all of them

The USER set these as ABSOLUTE rules for all clean-scheme work (stag2 and clean-cc alike). A run or change that violates one is not
a valid result, however fast or stable. USER: "those are absolute rules we have to satisfy"; on 7: "very hard but we should
definitely try".

1. **Energy conservation** (gas + radiation, boundary/source fluxes booked): per-step ledger at round-off on the bench decks; at
   production size, the hst ledger d(Etot) - int(L_in - L_top + P_esrc - P_sponge) dt close to old production (prod-1118 closes R136
   restarts to ~2e-6 of L_top T). Nothing dropped, carried silently or clipped without its counterpart.
2. **Momentum conservation** (gas + radiation F/c^2): conservative write-back; caps/limiters only on iterates or symmetric per face.
3. **L_z conservation** (flz q / flz Q dLz at round-off; rotating sets, ring, multi-D WB included).
4. **Restart bitwise:** N vs N/2 + restart + N/2 identical, and a forced dt_retry identical to the straight run; every new state in
   the restart file and the dt_retry checkpoint.
5. **Reproducibility:** same binary + input gives identical output run to run (GPU: no launch-order / atomic / reduction-order
   nondeterminism); decomposition independence within the documented tolerance.
6. **Mass conservation:** every density floor / floorcap / sponge change counted, with its mass, momentum and energy booked.
7. **Positivity + realisability by construction:** E > 0, |F| <= cE, eint > 0, rho > 0 from the scheme (M-matrix rows, symmetric
   per-face limiters, implicit stiff terms), not from clipping; every remaining clip counted, booked and listed as open.
8. **Well-balancing:** hydrostatic + radiative-equilibrium star stays static at the WB level (radial and multi-D WB inside the Roche
   lobe, none outside).
9. **Accuracy + order:** 2nd order in time where resolved, right spatial convergence, the reference tests (IDORT suite, Milne,
   waves), agreement with old production on the same restart within the time-discretisation error.
10. **Symmetry:** spherically / mirror-symmetric set-ups stay symmetric (no seam, pole, block-boundary artefacts).
11. **CPU/GPU consistency:** GPU vs CPU within the documented tolerance and below the dt-halving difference.
12. **No silent fallbacks:** unsupported combinations FATAL with a clear message.

## Request
Please check origin/clean-cc (your current head) against all twelve, and post a NOTE with a table rule | test used | result | open
violations. Where a rule has no test yet, say so (we found that ours had no forced-BE conservation row and no mass ledger).

## Violations we found on stag2 that may also hold on clean-cc (shared code)
- **time2_stage_local_be loses energy under scheme = clean:** Time2StageTheta (shared) zeroes the stage increment in inadmissible
  cells; the re-entry / carry (t2thr, t2thc, t2th_cnt_d) lives only in the old scheme's write-back, none in CleanWriteBack. R136 restart,
  150 cycles: ledger residual 8.4e-4 of L_top T vs old prod 2.2e-6. stag2 now FATALs the key under clean (fb2002c75).
- **Energy leak on the stage-fallback / k1_redist path:** R136 restarts with many stage fallbacks or k1_redist on millions of cells
  leak 5e-2 .. 0.49 of L_top T (per BE-after-fallback step 0.1-0.6 L_top dt), while runs without fallbacks close at 1e-4 .. 1e-6.
  Suspect: k1_redist's unpaid part / the fallback redo not booked under clean (being bisected; will report).
- **vet_gd_halo_compact = 2** keep test misses reads at 2x2 lateral decompositions: stale band values (2.5e-10) and non-bitwise
  restart at 4+ ranks (base vet_gd code).
- The stage face FCT (time2_face_fct) was inert under clean until stag2 49990612a; with it, R136 restarts run 45 stage / 2 BE in 24
  cycles with 0 retries and no energy dropped.

-- Orion
