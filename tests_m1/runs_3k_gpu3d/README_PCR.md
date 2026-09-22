# Speed-up #1: parallel cyclic reduction for the implicit-M1 x1 line solve

Date 2026-09-22, viper. Base: HEAD 92d36db3 (rt-integration). NOT COMMITTED.
Patch: `bench/m1_pcr_0922/pcr.patch` (applies with `git apply` to HEAD; files
`src/rad_m1/rad_m1.{hpp,cpp}`, `src/rad_m1/rad_m1_implicit.cpp`).
Snapshot and builds: `bench/m1_pcr_0922/{snap,build_cpu,build_gpu}`, GPU binary
`bench/m1_pcr_0922/athena_gpu` (md5 cb6f40ea...).

## What changed

- `ImplicitTridiagSolve` now dispatches. The old Thomas / cyclic-Thomas body moved unchanged
  into `ImplicitThomasSolve`. The new `ImplicitPCRSolve` does parallel cyclic reduction with
  one Kokkos team per (m,k,j) column:
  - The a/b/c/r rows load into team scratch, with the team's threads running along i. These
    loads are coalesced.
  - It runs ceil(log2 nx1) rounds (7 for nx1 = 84). The rounds are double-buffered, with one
    team barrier per round. Then x = r/b.
  - Periodic x1 uses the same Sherman-Morrison split as the Thomas code, with the second
    right-hand side carried through the same elimination.
  - It writes only M1_IW_S2.
  - The gathered stack sweep (`part_nblk > 1`) is unchanged and is always Thomas.
- New inputs in `<rad_m1>`:
  - `implicit_line_solver = thomas | pcr`. The default is thomas, so existing configurations
    stay bitwise unchanged.
  - `implicit_pcr_team`: the GPU team size. 0 means the next power of 2 >= nx1, capped at
    256. The host build uses AUTO.
  - `implicit_pcr_check`: a diagnostic. It runs both solvers on every call and prints
    max|x_pcr - x_thomas|/max|x|.
- Style: cpplint reports 69 errors on the three files both before and after the patch, so the
  patch adds none. No line exceeds 90 columns.

## Gate 1: CPU, 2-D seeded slab, 200 steps

Input: `runs_3g_newton_T/he_slab_m1_2d_V3edd.athinput` plus the three new keys, run as
`bench/m1_pcr_0922/cpu/he_slab.athinput`, `time/nlim=200`. Build: gcc/14 openmpi/5.0,
Release, 1 rank. Run directories: `bench/m1_pcr_0922/cpu/*`.

- **Line solutions:** max over 12 694 calls of max|x_pcr - x_thomas|/max|x| is **1.66e-15**.
  Run `ck_thomas` with `implicit_pcr_check = true` gives this number, and `ck_pcr` gives
  1.66e-15 too.
- **Iteration counts and state, compared with a round-off control.** The control is Thomas
  with `hydro/rad_flux_inner` changed by +/-1 ulp (runs `ctl_thomas_*`):

| run | Picard passes | inner its | breakdowns | hst max diff vs thomas (all rows / last row) |
|---|---|---|---|---|
| thomas | 805 | 6347 | 0 | 0 |
| pcr | 808 | 6361 | 1 | 5.0e-7 / 5.7e-6 |
| thomas, flux +1 ulp | 805 | 6345 | - | 3.4e-7 / 8.2e-6 |
| thomas, flux -1 ulp | 803 | 6335 | - | 5.1e-7 / 8.9e-6 |

  (hst columns are normalised by the column max; the last row is relative per column.) The
  seeded slab amplifies any round-off difference. The Thomas-vs-PCR spread in state and in
  iteration counts (+0.4 % Picard, +0.2 % inner) is the same size as a 1-ulp input
  perturbation of the Thomas run, so the counts are identical up to round-off, not bitwise.
  With `nx2 = 16`, 2 MeshBlocks, both solvers reproduce their 1-block counts exactly.
- **CPU cost:** 25.36 s (thomas) vs 25.62 s (pcr) wall for 200 steps. PCR is O(n log n) with
  a team size of 1 on the host, so keep thomas on CPU.

## Gate 2: GPU, apudev, 1 MI300A, job 11943021

`bench/m1_pcr_0922/submit1.sh`, using the same recipe as `bench/m1_gpu3d_0922/submit1.sh` arms
S and P1. Both solvers ran from the same binary, interleaved in the order S_T, S_P, P_T, P_P.
Logs are in `bench/m1_pcr_0922/runs/<arm>/run.log`, and `summarize.py` does the split.

| arm | grid | ms/cycle, cycles 0-250 | zone-cyc/s (whole run) | Picard mean | inner its / solve |
|---|---|---|---|---|---|
| S_T thomas | 84x52x52 | 315.1 | 7.20e5 | 4.42 | 22.5 |
| S_P pcr | 84x52x52 | **180.0 (1.75x)** | **1.22e6** | 4.65 | 22.0 |
| P_T thomas | 84x104x104 | 421.0 | 2.15e6 | 4.53 | 24.3 |
| P_P pcr | 84x104x104 | **269.4 (1.56x)** | **3.24e6** | 4.69 | 24.6 |

The ms/cycle column is (elapsed at cycle 250 - elapsed at cycle 0)/250, so both solvers are
timed over the same cycles. The whole-run counts cover different cycle ranges (1:30 or 2:00
walltime limits).

Same-cycle iteration counts come from the profiled arms (nlim = 20, rocprofv3 kernel stats in
`runs/PP_*/prof/rank_0_kernel_stats.csv`):

| arm | Picard passes | inner its | line-solve kernel | calls | avg per call | share of kernel time |
|---|---|---|---|---|---|---|
| PP_T thomas | 92 | 2771 | ImplicitThomasSolve | 5542 | **719.6 us** | 39.5 % |
| PP_P pcr (team 128) | 91 | 2763 | ImplicitPCRSolve | 5526 | **47.9 us** | 4.2 % |
| PP_P64 pcr (team 64) | 91 | 2763 | ImplicitPCRSolve | 5526 | 55.3 us | 4.8 % |

- **Kernel:** 15.0x faster, from 1.26 to 18.9 Gcell/s on 908k cells. The automatic team size
  (128 threads for nx1 = 84) beats 64.
- **Wall:** saves 152 ms/cycle on the production box, matching the RESULTS.txt estimate of
  ~155. The half-width box gains more (1.75x) because the Thomas cost did not scale with the
  number of columns.
- **Stability:** NON-CONVERGED stays 0 in every arm, and iteration counts over the same
  20 cycles agree to 0.3 %.
- **Default:** passes both gates, but left at thomas so that existing inputs and bitwise
  regression gates do not move. Flipping the default is one line (the GetOrAddString default
  in `ImplicitInit`). The recommendation is `implicit_line_solver = pcr` for every GPU M1
  input.

## Can the line factorisation be reused? (coordinator question)

**Where the coefficients are built.** In the HEAD line numbers of
`src/rad_m1/rad_m1_implicit.cpp`, TA/TB/TC are written only by the assembly kernel
`m1_impl_asm` (l. 3202, stores at l. 3333-3335), once per Picard pass. TB there carries:
- SRCB, the emission/absorption diagonal built from the lagged T' at l. 3195;
- TDIA, the transverse diagonal from `ImplicitTransverseTerms`, l. 1401.

**Across Picard passes:** they change, through SRCB(T') always and through the opacities
under `implicit_opac_update`.

**Within one BiCGStab solve:** they are fixed. `ImplicitPrecond` (l. 2246, 2265) overwrites
only TR. `ImplicitApplyOp` reads TA/TB/TC. Nothing writes them between assembly and the end of
the solve.

**Frequency.** From P1 in RESULTS.txt: 76 891 inner its x 2 preconditioner calls over 3117
Picard passes gives 49 applications per assembled matrix.

**Estimated saving.**
- One option is to factor once per pass and store the per-round PCR elimination factors
  (2 x 7 per row), so each application only updates the right-hand side. That removes the
  divisions and the a/b/c loads per round but adds global reads of the stored factors.
- The upper bound is the whole PCR kernel: 2 x ~115 inner its x 47.9 us ≈ 11 ms/cycle
  = 4 % of the 269 ms cycle. A realistic gain is about half of that, ~5 ms, or 2 %.
- Doing the same for Thomas would not help: that kernel is latency- and coalescing-bound,
  not division-bound.
- **Not worth doing now.** A cheaper gain of similar size: fuse the `m1_impl_prein` and
  `m1_impl_preout` copies of `ImplicitPrecond` into the PCR kernel (read rc and write zc
  directly). That removes 2 of the 3 launches per preconditioner call. This is inferred, not
  measured.
- After this patch the next targets are RESULTS.txt items 3 and 4: the host syncs and
  ImplicitHaloCopy. This patch did not change their absolute cost (not re-profiled), so they
  now take a larger share of the shorter cycle.
