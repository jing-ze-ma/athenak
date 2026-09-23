# Making the implicit-M1 speed-ups the defaults: gates and verdict

Date 2026-09-23, viper. Base: HEAD a6202226 (rt-integration). NOT COMMITTED.
- **Patch:** `bench/m1_defaults_0923/defaults.patch`. It passes `git apply --check` on HEAD.
  It touches `src/rad_m1/rad_m1.hpp` (comments only) and `src/rad_m1/rad_m1_implicit.cpp`
  (the GetOrAdd defaults and their comments). cpplint gives 800 errors on the two files
  both before and after the patch, so it adds none; no line is over 90 columns.
- **Work directory:** `bench/m1_defaults_0923/`, with `base` = HEAD and `new` = HEAD + patch.
  - Binaries are in `bin/`:
    - `athena_base_*`: HEAD, CPU and GPU (GPU md5 0e93b207).
    - `athena_new_*`: the first, unconditional patch (superseded).
    - `athena_new2_*`: the final patch (CPU only).
  - Build directories and every run live on `/viper/ptmp2/jinma/m1_defaults_0923/`, because
    the viper u2 inode quota was full (a test file loop failed at 496 files). `bench/.../cpu`
    and `bench/.../runs` are symlinks to ptmp2.
- **Tools:** `cmp.py` (statistics against a reference run), `arms_*.txt` (the CPU arm
  lists), and `gpuA-D.sh` + `gpu_common.sh` (the GPU jobs).
- **The switches.** ON = `implicit_line_solver=pcr implicit_bcg_sync=1
  implicit_lres_test=false implicit_conv_est=true implicit_predictor=step
  implicit_lin_ew_max=1e-2`. OFF = the HEAD defaults.
- **Method.** Every ON/OFF comparison uses the same binary, with the keys set on the command
  line. The inputs `inp/slab2d_k`, `inp/box3d_k` and `inp/vetm1_k` name the keys, because a
  command-line override of a key that is absent from the input file FATALs.
- **Round-off controls.** "+-ulp" means `hydro/rad_flux_inner` changed by 1 ulp
  (2475202000000000.5 and ...1999999999.5).

## Verdict

**Not every gate passes for the full set with every closure. Gate 4 fails for closure = m1
(Levermore).** Three findings:
- `implicit_lin_ew_max` alone raises the NON-CONVERGED steps.
- The full set raises them from 1-6 (the ulp controls) to 100-126.
- On the lag = pass slab, `implicit_lres_test=false` together with `conv_est` / EW stops
  the Picard loop before it has converged: KE1 is off by 2.3e-4, against 6e-9 for the ulp
  control.

The README_PICARD argument that "the lresid test is redundant under bicgstab" holds only
when the closure is frozen for the step.

**Shipped patch (the safe subset):**
- `implicit_line_solver = pcr` and `implicit_bcg_sync = 1` for every configuration.
- `lres_test = false`, `conv_est = true`, `predictor = step` and `lin_ew_max = 1e-2` only
  for closure = eddington | vet_sc | tau.
- m1 / minerbo / kershaw keep the old Picard defaults.
- EW defaults to 0 when `implicit_bcg_sync = 0` is set explicitly, so such an input does not
  FATAL.

For eddington, vet_sc and tau, every gate passes.

## Gates (all CPU runs: seeded 2-D slab, 84x32, tlim = 100, 621 steps, 1 rank unless stated)

Columns: Picard mean; total inner iterations; NON-CONVERGED (NC); relative differences from
the OFF run at t = 100 (F1top/Fin is the mean over t > 50). Source: `python3 cmp.py cpu ...`.

**(4) Closures.**

| closure | OFF: Picard / inner / NC | ON: Picard / inner / NC | ON vs OFF: F1top, KE1, KE2, totE | +-ulp control, largest |
|---|---|---|---|---|
| eddington | 4.07 / 33823 / 0 | 2.31 / 20869 / 0 | 1.2e-11, 2.1e-9, 1.0e-8, 1.9e-11 | 5.6e-12, 1.1e-8, 8.4e-9, 1.2e-12 |
| eddington + gas_newton + eos_cache | 4.64 / 34974 / 0 | 2.31 / 20890 / 0 | 2.3e-11, 7.3e-9, 3.3e-9, 1.1e-11 | 2.0e-12, 2.6e-9, 9.2e-9, 7.3e-13 |
| tau | 4.05 / 36847 / 0 | 2.82 / 20635 / 0 | 7.8e-11, 3.1e-8, 2.4e-8, 1.9e-11 | 1.1e-11, 7.6e-9, 2.3e-8, 1.8e-12 |
| vet_sc uniaxial | 4.07 / 35686 / 0 | 2.82 / 20148 / 0 | 4.6e-11, 2.8e-8, 1.1e-8, 1.2e-11 | 9.5e-12, 1.4e-8, 1.1e-8, 6.1e-12 |
| vet_sc full | 4.06 / 35940 / 0 | 2.82 / 19964 / 0 | 7.4e-11, 1.7e-8, 2.3e-8, 4.3e-12 | 1.1e-11, 5.9e-9, 2.5e-8, 4.7e-13 |
| **m1, lag = step** | 12.83 / 333092 / **6** | 24.68 / 237580 / **126** (1348 steps, not 689) | chaotic run, see below | ulp: NC 1-5, KE 3e-2 to 0.8 |

- These are the same tolerance-level shifts as in README_PICARD (F1top at most 8e-11, totE at
  most 2e-11). KE is within 1 to 4 times the ulp spread.
- The large "hstmax" numbers come from the 2-mom / 3-mom columns. Those are net momenta of
  about 1e-10 to 1e-12 of |1-mom|, i.e. noise, and they are excluded.

**m1 closure, each switch alone** (`m1s_*`: OFF plus one switch; `m1x_*`: ON minus one):

| arm | Picard | NC | note |
|---|---|---|---|
| pcr / bcg_sync / predictor / conv_est alone | 11.6 / 12.83 / 12.83 / 12.83 | 2 / 6 / 6 / 6 | sync, pred and conv are bitwise OFF (the predictor is inactive for m1) |
| lres_test=false alone | 7.92 | 0 | |
| lin_ew_max=1e-2 alone (+sync 1) | 18.38 | **25** | |
| ON without lres (lres_test=true) | 21.8 | **38** | |
| ON without conv_est | 13.1 | **11** | |
| ON without EW | 5.25 | 0, but 1516 steps | 505k inner iterations: dt collapsed |

The seeded m1 slab is a violent regime: F1top/Fin is 0.42 and the ulp controls differ by
5-80 % in KE. Only the counts are meaningful there.

**m1, lag = pass** (`tests_m1/runs_3j_vet/he_slab_vet.athinput`, nlim 150; `pm_*` arms):

| arm | Picard | KE1 end vs OFF |
|---|---|---|
| +-ulp | 79-83 | 6e-10 to 5.6e-9 |
| pcr alone | - | 1.1e-8 |
| lres=false alone | 11.6 | 3.6e-8 |
| **ON** | 12.8 | **2.3e-4** (F1top 4.1e-7) |
| ON without conv_est | - | 4.4e-6 |
| ON with lres_test=true | 91.8 | 6.6e-9 |

**(1) implicit_partition = gather** (eddington, x1 split into 2 blocks (nx1 = 42) on 1, 2
and 4 ranks (4 = also nx2 = 16), and into 4 blocks (nx1 = 21) on 1, 2 and 4 ranks):
- **pcr is not used under gather.** `ImplicitTridiagSolve` sends `part_nblk > 1` to
  `ImplicitGatherSolve` (the serial Thomas on the gathered column) before it looks at
  `impl_line_solver`. So ON with pcr and ON with `line_solver=thomas` are **bitwise
  identical**, on 2 blocks / 2 ranks and on 4 blocks / 4 ranks (`g_*_onT`).
- Every gather arm has NC = 0.
  - OFF: Picard 4.069-4.074, inner 33850-33900.
  - ON: Picard 2.306, inner 20793-20876.
  - This matches the 1-block runs (4.066 / 33823 and 2.306 / 20869) within the ulp spread
    (33823-33881).
- Differences from the 1-block OFF run: F1top at most 1.5e-11, KE at most 1.4e-8, totE at most
  1.8e-11. The ulp control on the gather layout gives 1.3e-12 and 5.6e-9.
- With the final patch, the patched binary on 4 blocks / 4 ranks, run on the plain input, is
  bitwise equal to ON (`nd2_x4b_r4`).

**(2) Multi-rank splits.**
- CPU, x2 split (nx2 = 16 on 2 ranks, nx2 = 8 on 4 ranks), eddington and vet_sc full: NC 0.
  Counts per rank count are within the ulp spread. The largest difference is KE1 2.7e-8
  (vet_sc, ON, 4 ranks).
- CPU, 3-D slab 84x16x8, x3 split (1 rank; nx3 = 4 on 2 ranks; nx3 = 4 and nx2 = 8 on
  4 ranks), tlim 20: NC 0. OFF Picard 4.20-4.24, ON 2.653. Differences at most 5e-11 on
  F1top and 1.8e-8 on KE1.
- GPU: see gate 2 below.

**(3) Boundaries and untouched paths.**
- No M1 input uses a non-periodic x2 or x3 boundary with an implicit transport. The only
  non-periodic x2 inputs are `rad_m1_beam` (vacuum) and `rad_m1_shadow` (reflect), and both
  are explicit.
- A reflecting-x2 test slab was run anyway (`rf_*`). ON vs OFF: F1top 1.3e-11, KE1 3.7e-9,
  totE 6e-13, NC 0. The ulp control gives 1.4e-12 and 2.9e-9.
- Every explicit and implicit_x1 input without bicgstab is **bitwise** unchanged by the final
  patch (gate 5). The exception is the 1-D implicit_x1 He columns, where pcr is active.

**(5) Every `<rad_m1>` input in the repo** (49 files: `inputs/hydro/he_box_m1_1d`, the 15
`inputs/tests/rad_m1_*`, and 33 under `tests_m1/`), nlim = 20, HEAD vs the final patch binary
(`cpu/g5/{base,new2}`):
- **No new FATAL.**
- Seven inputs FATAL identically on HEAD and on the patch, for reasons unrelated to the
  switches:
  - The two `rad_m1_radshock` inputs: `ref_m2.txt` is missing.
  - The five `pulse_md*` inputs: the x2 slice is out of the mesh.
- 32 of the 46 CPU-run inputs are bitwise identical (the 7 FATAL ones included). The rest differ at the solver level.
  - The m1-closure slabs keep Picard 77.9 -> 77.0 (pcr only).
  - The V3edd slabs go from 4.0 to 2.9.
- The three 3-D GPU inputs (3k, 3l, 3n sc_3d) did not reach cycle 1 in 900 s on one CPU core.
  On the GPU (job 11945189, committed input vs a copy with the ON keys written in) all three
  run: NC 0, Picard 4.55-4.60 -> 3.10.
- **Regression tests in `tst/` that touch rad_m1: none** (grep for rad_m1 / m1 over
  `tst/test_suite`). So no tst reference changes.

**(6) Long 3-D GPU run** (job 11945067): `he_slab_m1_3d` (84x104x104, eddington, gas Newton,
EOS cache). ON (GPU 0) and OFF (GPU 1) ran side by side on node vipa1327, same binary,
12.5 min wall.
- **ON:** reached tlim = 1000 in 6206 steps. **OFF:** reached t = 318 in 1973 steps.
- NON-CONVERGED 0 in both. Picard 2.85 (ON) vs 4.91 (OFF).
- ON, t = 500-1000: F1top/Fin mean 1.000002, range 0.999981-1.000014.
- Over the common time t = 0-318, with ON interpolated onto the OFF output times:
  - max relative differences: totE 1.2e-11, KE1 1.8e-8, KE2 1.6e-8;
  - F1top/Fin mean over t = 159-318: 1.00000213 in both.

**Restart.**
- OFF restart at t = 50: bitwise equal to the continuous run.
- ON restart: NOT bitwise. KE1 differs by at most 1.6e-9, totE by 7e-12, F1top by 1.8e-10.
- The cause is the predictor. Its increment is not in the restart file, so the first step
  after a restart is cold. ON with `predictor=none` restarts bitwise (`rs_onnp` / `rsr_onnp`).
- **Consequence of the flip:** M1 runs with eddington / vet_sc / tau are no longer
  restart-bitwise.
  - A restart-reproducibility gate should set `implicit_predictor = none`.
  - Alternatively, `ipred`, `pred_dt` and `pred_ok` would have to go into the restart
    (not done).

## GPU timing (apudev, MI300A, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, 3-D box, ms/cycle over cycles 20-120)

Jobs 11945065 (A) and 11945068 (B). The "b" arms are repeats; the GPU runs are
deterministic, so the repeats are bitwise identical.

| config | OFF (HEAD default) | ON |
|---|---|---|
| 1 GPU, eddington | 350.4, 367.1 | 115.8, 114.4 (**3.1x**) |
| 1 GPU, tau / vet_sc full | 357.8 / 354.1 | 110.6 / 117.0 |
| 2 GPUs, eddington | 324.5 | 93.5, 93.5 (**3.5x**) |
| 4 ranks on 2 GPUs (apudev MaxNodes = 1, so 2 nodes was not possible) | 328.2, 330.0 | 92.2, 94.2 |
| 4 ranks on 2 GPUs, vet_sc full | 340.6 | 107.8 |
| gather, x1 split into 2, 2 GPUs / 4 ranks | 294.9 / 314.1 | 158.8 / 165.2 (pcr unused) |
| long run (gate 6), 1 GPU each | 375.6 | 113.3 |

GPU accuracy at cycle 120 (ON vs OFF, 1, 2 and 4 ranks, gather, tau, vet_sc): F1top at most
1.2e-11, KE1 at most 2e-8, totE at most 7e-12, NC 0. The ulp controls give 2.2e-12, 2.8e-10
and 5.5e-12. The iteration counts are the same on 1, 2 and 4 ranks: ON Picard 2.858 with
7530-7571 inner iterations; OFF 4.57-4.58 with 13396-13426.

## What changes bitwise with the final patch

- Every implicit (multi-D bicgstab) input changes.
  - All 2-D and 3-D He slabs in `tests_m1/runs_3b4 ... 3n`, including the m1-closure ones:
    pcr and sync 1 are round-off.
  - The eddington / vet_sc / tau slabs change at the tolerance level.
- The 1-D implicit_x1 He columns change (`runs_3a/3b2/3g he_box_m1_1d_impl`): pcr.
- Every reference number quoted in earlier READMEs that used the HEAD defaults changes.
  Inputs that set the keys explicitly are unchanged. These are the 3k and 3l GPU inputs (pcr
  written in), 3n `he_slab_sc` (sets bcg_sync 0), and every picard / sync gate input.
- Explicit-transport inputs are bitwise unchanged, including all `inputs/tests/rad_m1_*`.
- No tst test changes.
- Bitwise identities checked on the CPU:
  - the new binary with the OFF keys equals HEAD;
  - the new binary with its defaults equals HEAD+ON for eddington, tau, vet_sc and gather
    4x4;
  - for m1 (lag step and lag pass) it equals HEAD with pcr + sync 1.
