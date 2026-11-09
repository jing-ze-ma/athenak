# TASK 2026-11-09 viper -> ORION: raycheap8 -- evolving-gas time-order test of the VET update + rebuild-vs-predict check

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap8.md`. All previous closure-update tests held the gas (so the
predicted, the stage-1 and the lagged builds saw identical inputs; RB == PR == lag by construction). Here the gas evolves.

## Binary
MOD = `fork/raycheap2-1109` head `4f227b777` (code = a1481b185, your raycheap6/7 build: reuse it).

## Part 1 (first, quick): does the stage-2 rebuild reach the operator? (any small evolving-gas cs deck you already have)
Take the HG cs small deck/state you used for the hgpos / hgblend gates (gas evolving, production keys). Run 20 cycles
twice with `rad_m1/time2_vet_col=predict` and with `rad_m1/time2_vet_col=rebuild` (same binary MOD, or the binary
of that TASK if MOD cannot run it), and compare the end states (`tests_raycheap/cmpbin.py` or your rst/bin compare).
Answer: bitwise equal (=> the stage-2 rebuild does NOT reach the operator: a real hole) or different (max rel diff).

## Part 2: stability pre-check, then the 4-level set
Decks `tests_raycheap/evo/evo_<cand>_c<cfl>.athinput`: the cs test star with gas_hold OFF, gas_feedback ON, coupling on,
an EQUILIBRIUM start (`e_out` uniform radiation, `gas_teq = true`: gas at T = (E/arad)^(1/4)), controlled heat
capacities, and a smooth +-20 % inner-flux modulation (period 40). t_end 20 (tlim); also run t_end 10 via
`time/tlim=10` on the command line. Candidates:
- A: kappa 1, ambient rho 1e-6, clump rho 50 (clump chi 50; radiation-dominated gas, follows E), e_out 1;
- B: kappa 1e-3, ambient rho 1e-3, clump rho 5e4 (clump chi 50; gas heat capacity >> aT^4, slow T), e_out 1;
- C: as A with e_out 0.1.
Pre-check (coarse level c0.15, arm PR = `rad_m1/time2_vet_col=predict`, and RB = `...=rebuild`), in the order A, B, C:
the first candidate that reaches tlim with rc 0, 0 FATAL, few or no NOT ADMISSIBLE, AND has RB != PR is used.
If none passes, adjust within the job and say what you tried (e.g. e_out 3, gas_rho 1e-4, a longer modulation period,
or the flux modulation off and only the equilibrium relaxation); report the failure modes.

Arms for the chosen candidate (command-line overrides on the evo deck; the deck itself is X9):
| arm | overrides |
|---|---|
| X9 | none (time2_vet_col = extrap2, x2_fields 9: production) |
| PR | rad_m1/time2_vet_col=predict |
| RB | rad_m1/time2_vet_col=rebuild |
| PRg | rad_m1/time2_vet_col=predict rad_m1/vet_source=gas |
| X9t | rad_m1/vet_gd_tangential=true |
| PRt | rad_m1/time2_vet_col=predict rad_m1/vet_gd_tangential=true |
Levels c0.15 / 0.075 / 0.0375 / 0.01875, 24 ranks, at t_end 10 and 20. Per arm and t_end:
`python3 tests_raycheap/conv.py <c0.15> <c0.075> <c0.0375> <c0.01875>`; plus RB vs PR and PRg vs PR max rel diff at the finest
level, walls, tensor-build counts, hesdirk2 counts. PASS (user's VET criterion): finest p >= 1.8 for PRt (or X9t).

## NOTE back
Part 1 answer; the pre-check log (what ran / failed); conv.py lines; differences; counts; walls. No code changes.
