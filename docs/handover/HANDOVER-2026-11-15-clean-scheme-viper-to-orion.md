# HANDOVER 2026-11-15: clean radiation scheme, viper -> Orion

From now on the Orion session owns the development of the clean implicit grey radiation scheme (code and CPU
tests). viper commits nothing more to `fork/clean-cc`. Times are PST.

## 1. Where the code is
| branch (fork) | head | what |
|---|---|---|
| `clean-cc` | dc67dba7b | YOURS. Clean scheme with both layouts: `clean_layout = staggered` (E cells, F faces) and `cell` (E, F/c cells). Cell is the layout under development. |
| `clean-1114` | c9d1a898e | base: clean staggered scheme + restyle modules 1-2 (cubed sphere, spherical polar). viper's restyle worker merges further modules here; merge `clean-1114` into `clean-cc` when you need them. |
| `clean-cc-geo` (local on viper, not pushed) | 51b860d62 | spherical-polar support for the cell layout (results in progress/wt_geo_PROGRESS.md); ask viper to push it when you get to geometry. |
| `clean-cc-pc` | merged | multigrid V-cycle preconditioner, already merged into clean-cc (7e4c9bacb). |

Files: `src/rad_m1/rad_m1_clean.{hpp,cpp}`, `rad_m1_clean_closure.cpp`, `rad_m1_clean_solver.cpp`,
`rad_m1_clean_cell.{hpp,cpp}`; the formal solution is `rad_m1_vet.cpp` (VetShortCharSweep). Build: built-in pgens,
MPI on, `-ffp-contract=off` as in your earlier TASKs. Key list with defaults: progress/wt_cc_PROGRESS.md, section
"STATE AT HANDOVER" (top of the file) — read it first.

## 2. The scheme (cell layout, current defaults)
- Unknowns E and y = F/c per cell; SDIRK2 (hesdirk2) implicit stages; gas momentum/energy updated conservatively
  from the same face fluxes (total energy and momentum conserved to round-off, verified: 1e-16).
- Face fluxes: E flux = central F - alpha (a/2)[E] + (1-alpha) Rhie-Chow; F flux = c^2 central(D E) - alpha (a/2)[F];
  a = c sqrt(D_max); alpha = A(tau) = (1-e^-tau)/tau of the face optical depth (tau enters ONLY the dissipation);
  jumps from van Leer PLM with lagged slopes, Picard passes folded into the Newton loop.
- Closure from the rays only: short-characteristics formal solution -> moments J, K, H+-, <|I-S|>; D = K/J,
  h+- = H+-/J, C1. Cells where the rays carry no radiation take a filled closure from neighbours.
  `clean_closure = extrap` (default): ray moments extrapolated in time from converged steps, plus ONE closure rebuild
  from the rays of the stage-1 iterate (stage 1 needs the closure of its own state Y1): ~2 formal solutions per step,
  time order 2.00. `iterate` (rebuild every Picard pass, ~9-13 solves/step) is the reference; `eddington` = fixed
  D = I/3, no rays (for the order gate).
- Positivity limiter: cells with E <= 0 or |F| > cE mark their faces low order (first-order states, alpha = 1, speed c,
  D = I on all stress rows); up to 3 re-solves, then backward Euler with marking, then BE with every face low order.
  Failed linear solves fail the stage (dt_retry).
- Boundaries: open faces use the rays' outgoing half-range stress (c k_out E) + exact prescribed inflow in the
  high-order solves; the all-low-order solve uses the upwind stress (b4dc0e75).
- Linear solver: right-preconditioned BiCGStab, multigrid V-cycle with block line relaxation (4-12x fewer iterations).

## 3. The user's rules (binding)
- Code style: original AthenaK style (upstream src/hydro on main): clean, human-readable, concise comments on every
  function and main step, NO history/patch/tag comments, no dead code, 90 columns, cpplint clean.
- Closure and regime ONLY from the rays; the optical depth chi*dx may appear ONLY in the numerical dissipation.
- Gates: (a) ORDER in space and time = the linear-wave (rw) suite; formal gate with FIXED D = 1/3 (as Jiang 2012,
  Jiang 2021, Ma 2025 do), both layouts, space ladders (N 64-512, nt = 2N) and time-only ladders (N256, nt 512-4096);
  the ray-closure wave (ref vetqs_d) is the stricter diagnostic. Amplitude 1e-6 for linear waves; time-only ladders
  may use 1e-4 with a 1e-5 linearity check (successive differences / amplitude must agree). (b) IDORT decks
  (nltek, hcyl, shd/shdg, cb/cbg/ws66) = ACCURACY gate: close to the exact solution, beams realizable; no order
  required. (c) Total energy AND momentum conserved to round-off on every test. Old production scheme is not a
  reference. Check a reference before trusting a floor (nltek's old reference was the Eddington approximation).
- Physical fixes preferred over numerical patches. No commits/PRs to main. Push branches to the fork only.

## 4. Status (numbers in progress/wt_cc_PROGRESS.md)
- rw ray closure: 2nd order in space (tau 0.1, 10) and time (P 3.3 and 1/3: 2.01/2.00 with extrap + 1 stage-1 rebuild).
- rw fixed D (eddington), cell layout, space: 1.87-2.36 in every case (P 1/3, 3.3; tau 0.1, 10, 1000). Staggered:
  invalid so far (its Newton loop stalls at clean_tol_nl 1e-13: use >= 1e-12). Time-only fixed-D ladders: not run.
- IDORT, cell layout (exact references in ana/): nltek1 9.4e-5 vs exact transport (= the reference's own accuracy,
  matches staggered); hcyl x0.5 0.047 (staggered 0.025); shd x0.25 0.264, shdg 0.256 (staggered ~0.02-0.03);
  cb/cbg/ws66 x0.5 1.68/1.62/1.35 and NOT realizable.
- 1-D prototype (proto2/README.md): cell scheme + limiter realizable to round-off in the full sweep; staggered with F
  evolved everywhere is more accurate on smooth waves but not realizable (88 violations) -> cell layout chosen.

## 5. Open problems, in priority order
1. 2-D beam realizability: every high-order stage solve on cb/cbg/ws66 is non-realizable (~3100-3300 bad cells;
   shd/shdg ~6500), so every step falls back to BE; after the all-low-order BE still 238/538/320 bad cells. Cause not
   found (filled-closure outflow faces were tested and refuted; the boundary stress explains only the low-order part).
   The 1-D prototype is realizable, so look for what is multi-D: transverse stress rows, the filled closure in cells
   without ray radiation, the low-order flux with D = I in 2-D.
2. shd 0.26 and hcyl 0.047 vs staggered 0.02 / 0.025 (references checked: hcyl exact formal solution converged to
   0.4 %; shd exact beam attenuation, thermal part < 0.5 %). Likely tied to 1 (limiter keeps fronts low order).
3. Fixed-D order matrix (both layouts, space + time-only) on the current head: the formal order gate.
4. Why production's one-solve predictor is 2nd order in production (prod-1118 on the same time ladder: 2.09/2.03)
   but 1st order in the clean scheme (cheap closure currently needs ~2 solves per step).
5. rw/conserve.py reads one row (slice_x2) of a 4-row mesh and misses x2 exchange: fix the script (code is fine).
6. Later: multi-MeshBlock support (prerequisite for production decks), spherical polar closure from vet_gd
   (clean-cc-geo), cubed sphere (plan in progress/wt_geo_PROGRESS.md), GPU validation (viper on request).

## 6. Package `docs/handover/handover_1115_clean/`
- `progress/`: wt_cc_PROGRESS.md (cell layout, all results; START HERE), wt_pc_PROGRESS.md (multigrid),
  wt_geo_PROGRESS.md (spherical polar + cubed-sphere plan), rad_clean_scheme_design.md.
- `proto2/`: the 1-D Python prototype (README.md with all results, scheme.py, stag.py, tests).
- `ana/`: scorers with the exact references: score.py (IDORT metrics), nlte_exact.py (exact two-level-atom transport,
  replaces the Eddington formula; the eps = 1e-4 reference is not converged to better than ~1 %), hcyl/shd exact npz.
- `rw/`: run_radwave.py (wave driver; ROOT, BINCONV and EXE are Raven paths: set RW_EXE/RW_RUNS/RW_PRODKEYS or edit),
  radwave_disp.py (dispersion relations: edd, vetqs, vetqs_d, transp), common.py, conserve.py (see problem 5),
  ladcmp.py, rwfarm.py + README.md (run-list front end).
- `farm_req/`: run lists for the current head (rw time/space ladders, IDORT incl. beam1d/beam45/beam45n).
- The IDORT decks + scorer package are in `task_1115_cc_idort_head/` (TASK of the same name).

## 7. Reporting
NOTE-<date>-orion-viper-<topic>.md on bsg-files-1009 for results and questions the user should see; viper relays and
can run GPU checks on request. The user may also talk to you directly.
