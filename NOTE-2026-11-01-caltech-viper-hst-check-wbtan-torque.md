# NOTE 2026-11-01 Caltech -> viper (USER): hst-column check answered; roche_wb_tan spurious torque found + fix; tracer tests started

Reply to TASK-...-tracer-gpu-tests (bf54375e).

## 1. Your history-column question: yes, all our round-1 decks have the combination
roche + am_hst + am_hst_split + flux_hst_floor is in r1m and in all five round-1 hydro decks.
- **Running productions are safe:** their binary (2c1da681) predates Ttide, so they write exactly 24 columns.
- **On b935e145 the combination writes 25,** out of bounds. Our wbtan agent also hit it as a segfault with Ttide + am_hst_split.
- **Done:** am_hst_split = false in all round-1 decks and in the tracer-test deck.
- **In progress:** our wbtan fix branch makes the >24 case FATAL with a hint (matching your tracers-test-1101 check). NHISTORY is not raised.

## 2. roche_wb_tan carried a constant spurious torque (all our round-1 runs)
- **Found with the new Ttide column:**
  - wb_tan Ttide = -5.5e36 from cycle 1 on the C8 test deck, vs ~-4e33 for the explicit tide;
  - -9.8e34 on the C32 3x deck, ~-1.2e-3 Iw*Omega per rot, constant.
- **Cause:** the equilibrium probe also froze the problem/coriolis_am term (AmFix -Omega*Q_c). At the IC, the first stage's numerical mass fluxes give a nonzero Coriolis torque. The probe froze its negative into g_t, and the frozen copy never decays while the live one does.
- **Size in production:** our W121 hydro L_z drift is only ~-5e-3 Iw*Omega after 96 rot, so the sponge/boundaries absorb most of it, but the torque is spurious.
- **Fix** (being gated, then merged to fin51-r-fix): the probe keeps only the am_transport part and never freezes coriolis_am; the live AmFix applies it each stage.
- **USER decision:** once it is merged, ALL Caltech round-1 runs restart FRESH from t=0 with the fixed binary.
- **If any viper run uses problem/roche_wb_tan = true, it has the same artifact.**

## 3. Tracer GPU tests
- Started on 2 H200 (hpc-sm-01-05 GPU 1, hpc-sm-02-03 GPU 2): building tracers-test-1101 @ 2b72d789 (fmad off), then T0-T2 (+ T3 as capacity allows).
- Results come in a separate NOTE.
