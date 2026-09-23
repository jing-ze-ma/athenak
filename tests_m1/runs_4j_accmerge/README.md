# runs_4j_accmerge: m1-accel merged into rt-integration, its levers default for be, and the inert conduction sweep skipped

Date 2026-09-24, viper. Branch `m1-accmerge`, worktree `/viper/ptmp2/jinma/wt_accmerge`.
Build and run tree: `/viper/ptmp2/jinma/accmerge_0923`. Every gate compares binaries built
from `git archive` snapshots there (`build.sh`, `gbuild.sh`). CPU gates ran on the login
node (nice 10, at most 12 ranks, `sched.py`). The GPU jobs ran on apudev with
`HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.

## Commits

| commit | content |
|---|---|
| 5d628e82 | merge m1-accel 0841ea5d (restart-header conflict resolved) |
| 11d3b7a2 | the runs_4a_accel levers default on for `time_scheme = be` |
| ac1e2892 | conduction: the inert tau-blend sweep is skipped when every weight is provably 0 |
| 2781b09d | `implicit_halo_overlap` defaults on where it is valid |
| 545706f2 | merge rt-integration 7066f13f (`implicit_bc_advect` default on) |

### Merge conflict

Conflicts were in `src/outputs/restart.cpp` and `src/pgen/pgen.cpp` only. rt-integration's
EINTRST1 header and m1-accel's M1ONEP01 header were both placed "behind the M1TIME2A one".
The resolution uses one fixed order, the same for writing and reading:

`M1PRED01 -> M1TIME2A -> M1ONEP01 (9 Reals, no slabs) -> EINTRST1 (last)`

- EINTRST1's w0(IEN) slabs stay the last slabs of the tail.
- Each header is optional and is detected by its magic. Files from rt-integration (no
  M1ONEP01) and from m1-accel (no EINTRST1) therefore both read.
- The predictor header now accepts 1-5 channels.

### Defaults (11d3b7a2, 2781b09d)

The levers default on only with `transport = implicit`, `time_scheme = be` and closure
eddington, vet_sc or tau:
- `implicit_fast_kernels = true`;
- `implicit_vimp_fold = true`, where implicit_vimp and implicit_op_stencil are on;
- `implicit_one_pass = 8` and `implicit_predictor_order = 2`, where the predictor is on.
  runs_4a_accel measured one_pass = 8 slightly faster than 4.

`implicit_halo_overlap = true` is the default where implicit_halo_mpi is on and the run has
more than one rank.

These defaults follow the implicit_enthalpy pattern:
- a restart (`global_variable::restart_run`) whose file lacks a key keeps the old value;
- the resolved values are always echoed into the parameter dump and restarts;
- explicit input values override.

hesdirk2, implicit_x1 and the other closures keep the old values.

### Conduction skip (ac1e2892): a code-level skip, bitwise

- **What a pgen sets.** A problem generator that reads nothing derived from `rad_tauf` sets
  `Conduction::rad_tauf_unread`. box_convection does so under `<rad_m1>` when the top
  sponge `vdamp_top_tau` is off.
- **The bound.** BuildRadWeights then bounds every face tau from above: it takes a parallel
  team sum of |kfac kappa rho dr| per column and adds a rounding margin of 1e-10.
- **When it skips.** If the bound is <= rad_tau_lo, every weight the sweep would write is
  exactly 0. rad_w is then zeroed once and the serial sweep is skipped.
- **What it does not cover.** A NaN or inf in the sum falls back to the sweep. rad_tauf is
  not filled while the sweep is skipped, and the dt-collapse diagnostic prints -1 for it.
- **Why not drop the conduction block.** Dropping the block would change box_convection's
  wall_noflux energy form (`diff_flux_`), so that route is not bitwise.

## Gate table (CPU unless stated)

The slab is the 2-D He slab: `slab2d_plm_vimp_*`, plm + vimp, vpert 1e-2, t = 100, with
restarts at t = 50. The box is the 3-D box: 84x32x32 in 4 blocks, 40 cycles. The radwave
runs are the runs_4a_accel G0 inputs. The compare is hst/tab byte for byte; bin and rst
are compared on the payload behind `<par_end>`, because the echoed keys change the
parameter header. Summaries: `GATES.txt`, `GATES_E.txt`.

| gate | result |
|---|---|
| merge (5d628e82, levers unnamed) vs 5c9e4432: slab be Eddington 1/2 ranks, vet_sc 2 ranks, hesdirk2 1 rank, hesdirk2 + vet_sc 2 ranks, box 1/4 ranks | PASS, 7/7 runs SAME, 0 DIFF |
| final (ac1e2892), levers named off, vs 5c9e4432: same 7 runs + radwave edd be, vet be, vet hesdirk2 | PASS, 10/10 SAME (34 tab files per radwave run) |
| final (545706f2), levers + halo_overlap named off, vs 7066f13f: slab e1/e2/e4/v2/hv2, box 2/4 ranks | PASS, 7/7 SAME |
| conduction skip in every M1 fin run | all builds skipped (e.g. 1243 skipped / 0 swept on the slab, 81/0 on the box); bitwise per the rows above |
| new defaults (no keys): slab e1/e2/e4/v2 and box 1/2/4 ranks | NON-CONVERGED 0; positivity, line_jacobi, vimp and newton fallbacks 0; one_pass period 8 echoed |
| new defaults, radwave (be) | NON-CONVERGED 0; last m1 tab vs levers off within 6.0e-12 (edd) and 3.5e-13 (vet), relative. edd_fd shows line_jacobi fallbacks = 329, a count of the radwave line solver |
| halo_overlap default vs named false (levers on in both), 2/4 ranks | round-off: last-row KE <= 6e-8 and tot-E <= 6e-10 relative (transverse momenta ~1e-2 relative, noise around zero); NON-CONVERGED 0 |
| restart, levers on: be Eddington 1 rank, be vet_sc 2 ranks, hesdirk2 + L 1 rank, hesdirk2 + vet_sc + L 2 ranks; plus levers named off | PASS: t = 100 bin and rst identical; 50/50 restarted hst rows identical |
| restart with halo_overlap on (2 ranks, be Eddington and vet_sc) | PASS, same criteria |
| old restart files (5c9e4432 slab be 1 rank, hesdirk2 1 rank, vet_sc 2 ranks; 7066f13f slab 2 ranks) read by the final binary | PASS: continues bitwise with the old behaviour (one_pass 0 and predictor_order 1 echoed), no warnings |
| deep hot Jupiter, read in place: hydro `bench/cs_hyd4_prod/rst/dhj.00147.rst` (t = 2.2418e7 s, 73.5 rot), MHD `bench/cs_mhd_prod4/rst/dhj.00115.rst` (t = 1.7538e7 s, 57.5 rot; P = 3.05e5 s), 6 ranks, 3-4 cycles | PASS: ac1e2892 vs 5c9e4432 7/7 SAME (hydro) and 7/7 SAME (MHD) |
| MHD dhj restart (EINTRST1 path), restarted from its own cycle-2 file | PASS: rst and bin identical, 3/3 hst rows identical |
| GPU job 11955024 (1 GPU, full box 84x104x104): ac1e2892, levers named off, vs 5c9e4432 | PASS: be 15/15 SAME, hesdirk2 15/15 SAME |
| GPU job 11955208: 545706f2 vs 7066f13f, same gate + timing | see `gpu2/log.out.11955208`, `python3 gsum.py gpu2` |

`rt_profile.bin` (box_convection's horizontally averaged profile dump) does not match across
a restart, with levers on or off and for files written by 5c9e4432. Its dump schedule and
averages are not carried by the restart file. It is a diagnostic; the state (bin/rst/hst)
is bitwise.

## GPU cost (job 11955024, 1 GPU, 3-D He box, 120 cycles, ms/cycle over cycles 20-120)

The arms are interleaved and repeat 2 runs in reverse order. "before" = 5c9e4432 defaults;
"after" = ac1e2892 defaults.

| arm | ms/cycle (a, b, c) | mean |
|---|---|---|
| be Eddington, before | 51.8, 52.7, 52.1 | 52.2 |
| be Eddington, after, levers named off (conduction skip only) | 46.6, 48.4, 49.0 | 48.0 |
| be Eddington, after (defaults) | 41.3, 39.7, 40.6 | **40.5** |
| be vet_sc, before | 57.7, 58.2, 58.3 | 58.1 |
| be vet_sc, after (defaults) | 55.6, 48.8, 46.7 | **50.4** |

- Picard passes per solve fall from 2.125 to 1.233 (Eddington) and from 2.733 to 1.242
  (vet_sc).
- The conduction skip is worth about 4.2 ms/cycle here (52.2 -> 48.0).
- NON-CONVERGED 0 in every arm.
- **Flag.** With the new defaults the per-cell gas-coupling Newton falls back to the
  bracketed root find 1 (Eddington) and 2 (vet_sc) times in about 1.35e8 cell-passes; the
  old defaults have 0. m1-accel's own be + L arms show the same 1-2 (accel_0923/runs/*beL*).
  All implicit-solve fallbacks are 0.
