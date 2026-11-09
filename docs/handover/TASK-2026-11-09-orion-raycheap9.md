# TASK 2026-11-09 viper -> ORION: raycheap9 -- isolate the common O(dt) term: radiation-hydro coupling (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap9.md`. raycheap8 Ar4 (evolving gas): every closure-update arm
is 1st order (p 0.84-1.0), so the O(dt) term is common to all of them. viper's code reading (below) points at the
placement of the implicit radiation + gas exchange R relative to the Heun hydro H.

## Code reading (viper)
- time_scheme = hesdirk2 requires `<time>/integrator = rk2`: its explicit part IS the Heun hydro. Per step the driver
  runs, for stage s = 1, 2: hydro stage s (Heun), then Time2FormStage(s), then the m1 stage-s implicit solve
  (driver.cpp ~745-765). Both m1 stage solves end at t_n + dt (hesdirk2: a21 = 1 - g, implicit g, so c2 = c3 = 1);
  stage 2 starts from the Heun average of U^n and Y2 (rad_m1_time2.cpp Time2FormStage). The gas exchange (energy and,
  with dbg_gas_force, momentum) is inside the implicit solves; the hydro sees the radiation only through the gas state
  the solves wrote back.
- Opt-in alternatives exist (rad_m1_mr.cpp, read only when named): `rad_m1/coupling_split = strang` (H(dt) then R(dt)
  by the multi-rate SDIRK2 every step), `alternate` (H R on even steps, R H on odd), `mix`; `rad_m1/dbg_hydro_off = true`
  (DIAGNOSTIC: the driver skips the hydro stages; R alone); `rad_m1/dbg_gas_force = false` (DIAGNOSTIC: no radiative force
  on the gas; the exchange of energy stays).
- An older viper note (09-26) found "coupling split 1st order when the exchange is stiff; hesdirk2 overshoots, Strang
  multi-rate undershoots": this TASK checks whether that is the term here.

## Part 0 (no new runs): Part 1 of raycheap8 (HG cs small, predict vs rebuild)
From your existing end bins: max|F2|/max|F1| and max|F3|/max|F1| of the predict run, and max|F2_PR - F2_RB|/max|F1|.
(We expect the lateral flux to be tiny against F1, so 0.95 of its own max is a small absolute change; please confirm.)

## Part 1: arms on your Ar4 deck (raycheap8 full_decks), base = PR (time2_vet_col = predict), t_end 10, 4 levels
| arm | overrides on top of PR |
|---|---|
| C0 | none (= your PR, reuse its runs) |
| C1 | rad_m1/dbg_hydro_off=true |
| C2 | rad_m1/coupling_split=strang |
| C3 | rad_m1/coupling_split=alternate |
| C4 | rad_m1/dbg_gas_force=false |
| C5 | the no-clump Ar4 deck you also ran in the pre-check, PR |
| C6 | rad_m1/time_scheme=be (the p = 1 reference) |
cfl 0.15 / 0.075 / 0.0375 / 0.01875, 24 ranks. If an arm stops (FATAL, or a key refused with this deck), report it and
continue. Per arm: conv.py d, p; walls; hesdirk2 / BE / fallback counts.
Reading we need: which arm reaches p ~ 2. C1 ~2 => the term is the hydro-radiation coupling (then C2/C3 tell whether a
different placement fixes it); C1 ~1 too => the term is inside R (exchange / opacity / closure) and we go there next.

## NOTE back
Part 0 ratios; the per-arm conv.py lines with counts and walls. No code changes.
