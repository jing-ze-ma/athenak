# T1 stall: why `ck_implicit` on tm does not converge in the prod4 configuration (2026-09-22)

Follow-up to `README_T1.md` §3. Nothing is committed. The work is in
`bench/impl_stall_0922/`: a `git archive HEAD` (a6cc3298) snapshot plus `t1.patch` plus a
diagnostics patch, the CPU arms in `cpu/`, the GPU arms in `gpu/` and the bitwise gate in
`gate/`. Restart: `cs_mhd_prod3/rst/dhj.00567.rst`, with the prod4 `<problem>` block
(`impl.athinput`: spherical, beam_sph, tm).

## Verdict

* **Cause: the two-level split `ck_impl_arat` (default 2).** The top cell of the day-side
  columns changes class every pass: thin, then thick, then thin again. The thin solve
  (`CkThinSolve`, `E ~ (e/e_n)^4`) and the Newton row (`c_v = e/T`) both badly
  underestimate how steeply that cell's emission depends on e. Each overshoots, and the
  class flips back. The result is a period-2 limit cycle, not slow convergence.
* **Fix: an input change, with no code change.** Set `problem/ck_impl_arat = 1.0e30`,
  which turns the split off. On GPU (same binary, 400 calls, both ranks counted):
  * non-converged calls fall from 400/400 to 122/400;
  * the gap falls from max 1.7e-4 / median 6.4e-6 to max 5.4e-6 / median 6.6e-7;
  * the cost falls from 5.41x to 4.66x the semi-implicit step.
* **What still fails, and cannot be fixed locally:** one night-side column (near the
  terminator) that has an odd-even temperature structure in the restart state. Its cold
  cells have a net source below minus their own emission (`A = S + E < 0`). Their
  backward-Euler root lies beyond the `ck_impl_demax` bound, and removing the bound does
  not make them converge. This column accounts for all remaining non-converged calls
  (measured on CPU) and for the gap max. **Gate (c), gap <= 1e-6 max, is therefore still
  not met.**
* Job 11943097 (maxit 20, started by the previous agent): tm `i20` and four-pass `f20`
  both finish with 200/200 calls not converged. `i20`: res max 5.62e-3 / median
  1.92e-5, gap max 1.24e-3 / median 3.6e-6. `f20`: res max 5.69e-3, gap max 1.14e-3.
  `ana.py` over `bench/impl_t1_0922/gpu/{i20,f20}`. More passes do not help.

## 1. Evidence

Diagnostics come from `ck_impl_debug = 2`: each rank writes `ckdbg.<rank>.txt`, with the
residual split by class, alternative norms, the class-flip count, the worst-residual and
worst-step cells, and one tracked column row by row. They were taken on CPU with 24 ranks
(one MeshBlock each), maxit 20. Summaries use `cpu/dsum.py <arm>`.

**The Newton loop is rank-local.** `CkImplStep` does not reduce over MPI, so the
`### ck_implicit` line reports rank 0's verdict only. With `ck_impl_debug = -1` every rank
prints its own line, and `gpu/ana.py` counts a call as converged only if both ranks
converged.

**(a) Where the stall is.** Arm `cpu/d1`, job 11943210, call 1, passes 18/19. The
largest steps are in the top cell (`i = ie = 129`) of day-side columns on 8 ranks. In
those cells:

* the beam deposit is `Qs` = 1.1-1.3e-3;
* p ~ 1e-8 bar, well above pfloor (1e-9 bar);
* the cells are far above the ck cut (`icut = 27`).

The same cell (rank 13, m 0, k 17, j 8) alternates between two states:

| pass | e | T [K] | p [bar] | E | A/E | class used |
| --- | --- | --- | --- | --- | --- | --- |
| 18 | 0.0930 | 2843 | 9.4e-9 | 4.86e-5 | 26.5 | thin |
| 19 | 0.0993 | 4055 | 1.34e-8 | 1.57e-3 | 0.94 | thick |

A 6.8 % change in e gives a 43 % change in T and a 32x change in E. From these two
iterates:

* `d ln T/d ln e = 5.4`. The gas is atomic hydrogen, and its e carries the H2
  dissociation energy (~2.2 eV/H against 1.5 kT = 0.45 eV at 3500 K), so the true c_v is
  about e/(5.4 T).
* `d ln E/d ln T = 9.8`. These are Wien-side bands.
* `d ln E/d ln e = 53`. The thin solve assumes 4, and the Newton row, with `c_v = e/T`,
  assumes about 10.

The class flips every pass. In `d1` call 0, 214-232 cells flip on every pass from pass 5
on, and `dstep` holds at 0.19 / 0.24 without shrinking. In `P2` call 2 (job 11943254),
19000 cells flip and 21000 cells take steps larger than 1 %.

**(d) The discontinuity is the split, and removing it fixes the day side.** The runs are
5 cycles = 10 calls, maxit 20 (job 11943294, `cpu/{P5,A5}`, per-rank lines):

| arm | ranks with a non-converged call | other ranks: passes mean / max, res max |
| --- | --- | --- |
| P5, `arat = 2` (production) | 10 ranks; 9 of them fail 10/10 calls | res max 5.7e-3 |
| A5, `arat = 1e30` | **rank 13 only** (9/10 calls) | 7.74 / 17, **<= 1.0e-8** |

A hysteresis variant was also tried: a cell stays thick for the rest of the call once it
has been thick (`cpu/L2`, job 11943254). It removed the flips (0 per pass), but the cells
that stay thin still oscillate at the 0.5 step cap (about 18700 cells take steps larger
than 1 % in call 2). The thin solve itself is unstable for these cells, so this variant
was dropped.

**(b) Incomplete Jacobian: not the cause.** With the split off, the T1 tridiagonal
(nearest-neighbour entries, negative per-chain entries dropped) converges every column
except the odd-even one. Tracked in `A5`, and in `a30` of job 11943210: 5-11 passes in
calls 1-3. The runs that include the dropped entries were not needed, and were not done.

**(c) The norm: not a floor.** In stalled passes the diagonally scaled norm
`|R|/((1+4 bdt E/e)(e+eps e_max))` is within about 15 % of the production norm (for
example 4.58e-4 against 5.18e-4 in `d1`). The per-cell steps stay at 11-24 %. Columns
that converge reach 1e-8 in 5-11 passes, so the tolerance is reachable.

**(e) What prod4 adds.** Two arms isolate it (job 11943210, `cpu/b0` and `cpu/f0`):

* With `ck_beam_sph = false`, only 58-64 cells are thin instead of about 83000. The split
  is then almost inert, and every call converges except for the rank-13 column. The same
  holds for the four-pass form without the beam (`f0`).
* The spherical beam heats the tenuous day-side top to A/E >> 2, which puts those cells
  on the thin branch. That is why phase 4 (`ck_beam_sph = false`) converged and T1 does
  not.

**The column that remains** is rank 13 (gid 13 in a 24-rank run; rank 1 of the GPU runs),
m 0, the MeshBlock corner column k = ks, j = je. It is night side with a grazing beam at
the top. For i = 65-70 (p 2-17e-6 bar), T already alternates cell to cell in the
restart's e^n: 1838 / 4986 / 5417 K at i = 69 / 68 / 70. The cold cells have `A/E` from
-24 to -107. The sweep takes more out of them than their own emission, which is what a
linear-in-tau source function gives between hot neighbours.

Cell i = 69 in stage 1 (`bdt` = 12.9 s):

* its backward-Euler root is below `0.5 e*`, so it is pinned at the `ck_impl_demax`
  clip: `e = 7.408 = e*/2` on every pass, residual 5.2e-4;
* with `ck_impl_demax = 0` (arm `A5d`, and `A2d` with the column tracked) calls 0, 1 and
  3 converge only linearly;
* in call 2, cells are driven toward zero with `|R|/e` up to 67.

This is a defect of the state and the discretisation, not of the solver.
**Recommendation: look at this seam/corner column in the production run separately.** It
is not measured here whether the semi-implicit path clips it in the same way.

## 2. Gate numbers, GPU (job 11943285, `bench/impl_stall_0922/gpu/`)

`apudev`, 2 ranks. One binary for every arm: `athena.gpu.dbg2`, which is T1 plus the
diagnostics patch and the inert `ck_impl_latch` knob (default off, removed from the
patch afterwards); the diagnostics are inert unless `ck_impl_debug` is not 0. 200 cycles = 400 RT calls from
the rot-283 restart. The arms were interleaved as s, i1, a1, s, a1, i1. All implicit arms
use maxit 8, `reuse_jac = 1`, `seed = 2`. `ana.py` counts both ranks.

| arm | setting | cpu time used, r1 / r2 [s] | x semi | non-conv calls | passes mean (max over ranks) | res max / median | gap max / median |
| --- | --- | --- | --- | --- | --- | --- | --- |
| s | implicit off | 15.04 / 15.13 | 1.00 | - | - | - | - |
| i1 | T1 as briefed (`arat = 2`) | 81.60 / 81.52 | **5.41x** | **400 / 400** | 8.00 (cap) | 1.15e-3 / 1.74e-5 | 1.72e-4 / 6.35e-6 |
| a1 | + `ck_impl_arat = 1e30` | 70.52 / 70.10 | **4.66x** | **122 / 400** | 6.67 | 1.09e-3 / 1.0e-8 | 5.37e-6 / 6.58e-7 |

* In a1, rank 0 has 19 of 400 calls not converged, all near misses (res <= 3.1e-7), with
  a mean of 5.89 passes and a gap of max 1.76e-6 / median 1.83e-7. Rank 1 has 122, and it
  holds gid 13 (the odd-even column).
* In a1 the gap median meets the 1e-6 target and the max does not: 5.4e-6, set by the
  odd-even column. The mean of 6.67 passes compares with 5.67 in phase 4, which was run on
  a different configuration (no beam_sph, four-pass, frozen_op).

## 3. Patch and bitwise gate

The fix needs no code. Put `ck_impl_arat = 1.0e30` next to `ck_implicit` in the prod4
`<problem>` block. The split exists for the 10-100x dt cold-start tests of phase 2, and
prod4's dt is hydro-bound (DESIGN_tm §0).

`bench/impl_stall_0922/stall_debug.patch` goes on top of T1 and is diagnostics only:

* `ck_impl_debug = 2` writes the per-rank dumps; `ck_dbg_rank/m/k/j` picks the tracked
  column;
* `ck_impl_debug = -1` prints the report line from every rank, tagged `rank=`;
* the rank-0 line gains ` rank=0`.

**Off path** (CPU, `gate/gate.sh`, `gate.log`, 20 cycles; `.hst` bytewise, `.bin` and
`.rst` payloads compared past the parameter dump, since the new `GetOrAdd` keys make the
dump longer) against `impl_t1_0922/athena.cpu.new`:

| gate | settings | result |
| --- | --- | --- |
| `a_prod_tm` | production form, implicit off | BITWISE (25 files) |
| `t_imp_tm_rj` | implicit tm, `reuse_jac = 1`, `seed = 2` | BITWISE (25 files) |

The latch knob was tested (`L2`) and then removed; it is not in the patch.

## 4. Files

`bench/impl_stall_0922/`:

* `build.sh`, and the binaries `athena.cpu.dbg`, `.dbg2` and `.final`, `athena.gpu.dbg2`;
* `stall_debug.patch`;
* `cpu/`: `submit.sh`, `arms.sh`, `dsum.py`, and the arms `d1 a30 b0 f0` (job 11943210),
  `L2 A2 P2` (job 11943254) and `A2d A5 A5d P5` (job 11943294);
* `gpu/`: `submit.sh`, `ana.py`, `log.out.11943285`, and the arms `s_r1 s_r2 i1_r1 i1_r2
  a1_r1 a1_r2`;
* `gate/`.
