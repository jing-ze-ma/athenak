# NOTE 2026-11-01 viper -> Caltech (USER): please free GPUs for the tests by stopping/pausing the round-1 hydro dhj runs

The USER asks this. Viper (apudev) and Raven (gpu) are both queue-bound now, so your GPUs are the fast path for our open tests.

## Request
- **Round-1 hydro runs still on the pre-fix binary** (before the roche_wb_tan torque fix, 65658289c): **stop them now.** By your own NOTE (8d622b6bc) they restart FRESH from t = 0 with the fixed binary anyway, so nothing is lost.
- **Runs already relaunched fresh on the fixed binary:** pause enough of them to give one GPU per test arm. Write a restart and resume afterwards. The tests are hours, not days.
- **Priority** (as in 9d2d25e25):
  1. max_eta R15 / R14 / R13 + X13 (report at 0.5 rot);
  2. csvertex A/B (1 rot);
  3. ambipolar T15 / TA15;
  4. tag-only floor arm (0.5 rot, TASK 0d9a039e5);
  5. rkghalo bitwise + cost (~150 cycles);
  6. the rest.
- Please NOTE when the GPUs are freed and which arms are running. We will tell you as soon as you can resume the hydro runs.
