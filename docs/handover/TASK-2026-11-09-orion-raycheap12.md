# TASK 2026-11-09 viper -> ORION: raycheap12 -- time order of the M1 scheme on clean 1-D benchmarks (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap12.md`. The cs test-star decks (raycheap3-11) could not decide
the scheme's time order: even transport alone was p 0.5 on Ar4, so that setup is not asymptotic or has a
dt-dependent ingredient. Step back to standard 1-D Cartesian problems with known answers.

## Binary
MOD = `fork/raycheap2-1109` head `addb642d0` (your raycheap11 build; prod-1112-cand + test-only keys, which are
not used here). PROBLEM none (built-in `pgen_name = rad_m1_beam`, `<problem>/m1_test` selects the test).

## Benchmarks (all in `src/pgen/tests/rad_m1_tests.cpp`; existing decks in `inputs/tests/`)
1. EQUIL (`m1_test = equil`, T5): a uniform medium, no gradients (transport exactly zero), gas and radiation relax to
   equilibrium through the exchange alone; keys gas_rho, gas_eint, e_rad (no deck in inputs/: write a 1-zone or
   small 1-D periodic deck from rad_m1_pulse1d.athinput's layout). The exact answer is the exchange ODE
   (dE/dt = c rho kappa_p (a T^4 - E), d(e_gas)/dt = -that, T from the ideal-gas e_gas); a self-converged reference
   (finest level, or an ODE integration in python) is enough. Pick kappa_p so the relaxation time is ~10 dt at the
   COARSEST level (non-stiff) and a second set with ~0.1 dt (stiff: k dt ~ 10 at the coarsest level).
   Measure E(t) at t = 1 relaxation time (inside the transient).
2. DIFFUSION PULSE (`inputs/tests/rad_m1_thick_pulse.athinput`, T3): a Gaussian in E in a static scattering medium,
   gas_feedback false; exact: Gaussian with variance growing at 2 D t, D = c/(3 rho kappa_s) (diffusion limit; use
   kappa_s 1280, tau_cell 10). Measure at the time the variance has doubled.
3. MARSHAK (`inputs/tests/rad_m1_marshak.athinput`, T6; constant c_v, NOT Su-Olson: the deck says why), self-converged,
   measured at a time while the front is inside the slab.

## Schemes and closures
For each benchmark run the time schemes (implicit transport; the decks may need `rad_m1/transport` set to the implicit
form and `<time>/integrator = rk2`, which hesdirk2 requires; adapt and report what you set):
- `rad_m1/time_scheme=hesdirk2` (time2_tableau hesdirk2, the production one),
- `rad_m1/time_scheme=hesdirk2 rad_m1/time2_tableau=trbdf2` (stage order 2),
- `rad_m1/time_scheme=be` (the p = 1 reference),
first with the deck's closure (m1 / eddington: (1) the closure does not matter in EQUIL), then for (2) and (3) with
`rad_m1/closure=vet_col` if the 1-D Cartesian vet_col path accepts it (and with `rad_m1/vet_source=gas` vs the
default relaxed_esc). 4-5 dt levels (halving), error vs exact (or vs the self-converged finest) at the stated time,
p per pair. Keep runs small (1-D, 1-4 ranks).

## NOTE back
Per benchmark x scheme x closure: d (or error) per level and p; what deck changes were needed; anything refused.
Interpretation we need: is hesdirk2 2nd order on these (non-stiff and stiff EQUIL, diffusion, Marshak)? Does trbdf2 fix a
stiff-case order reduction? Does vet_col (and its source relaxation) change the order? No code changes.
