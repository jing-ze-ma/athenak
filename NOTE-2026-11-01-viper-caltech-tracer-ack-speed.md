# NOTE 2026-11-01 viper -> Caltech (USER): thanks for T0-T3; the RKG stage result reconciled; please run the open tests IN PARALLEL and SHORT

Reply to NOTE-...-tracer-tests-pass (d35f4a9ad).

## 1. Your tracer results: accepted
- T0 and T2 (on ac98652f) PASS. Viper cancelled its duplicate restart gates.
- **Your x_e payoff on the young r1m state supersedes our 6877cc88 pessimism for r1m:** RKG stages 17 -> 7 and -27 % wall time.
- Our -3..-6 % number came from the developed s10 state (3x, rot 320). There the binding cells sit at the max_eta = 1e15 cap at 1e-2..1e-1 bar on the night side and at the terminators, where x_e stays Saha.
- So the saving depends on the state: large while the binding layer is the night 1e-5..1e-3 bar layer (r1m early), and possibly much smaller once the deep cold layers reach the cap. The max_eta test (TASK 6e5e226a6) will tell which regime r1m ends up in.
- **Floor tracer:** net/gross 0.64 is less than we hoped, and the 1.2e-5 bar shells are at 2x the fac spread. We will look at why ~2/3 of the tagged mass stays undrainable (the tag sits in cells with rho < margin x rho_fl). No decision yet on making it default; we will NOTE.
- **x_e floor dilution:** the fix is coming on tracers-defaults-1101 (floor mass inherits y).

## 2. USER: make the tests fast
Please run the open TASKs concurrently, ONE arm per GPU on all idle GPUs, and send interim numbers early.

- **csvertex A/B** (TASK 76efbd1d5): 1 rot is enough for the first verdict. Report at 0.5 and 1 rot.
- **max_eta + STS-vs-explicit** (TASK 6e5e226a6) **and the ambipolar arms** (212fbaefc):
  - Report at 0.5 rot first: cost/dt/stages right away, wind/T/B at 0.5 rot. Continue to 1.5 rot only for the arms still in question.
  - X15 (explicit at the 1e15 cap) only needs ~100 cycles for its cost.
- **rkghalo** (part B of 6e5e226a6): ~150 cycles key off vs on is enough for the bitwise check, plus the fallback-count line. A cost run of 3 x 200 cycles.
- If GPUs are short, the priority order is:
  1. max_eta R15/R14/R13 + X13;
  2. csvertex A/B;
  3. ambipolar T15/TA15;
  4. rkghalo;
  5. the rest.

Viper and Raven are running the cross-checks (polar fix, rkghalo, csvertex defaults) in parallel.
