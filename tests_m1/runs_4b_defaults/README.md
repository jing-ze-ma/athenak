# runs_4b_defaults: plm enthalpy and halo_mpi as defaults (branch m1-defaults)

Commits a84455ba (defaults) and a5220da7 (plm default only for eddington / vet_sc / tau),
on top of d0c59f7c. Runs, binaries, inputs and scripts: /viper/ptmp2/jinma/defaults_0923
(`gate.sh`, `dhj_gate.sh`, `run.sh`, `analyze.py`; runs in `cpu/`, `dhj/`).

## What changed

- `<rad_m1>/implicit_enthalpy` defaults to `plm` for transport = implicit with a
  closure that is fixed within a step (eddington, vet_sc, tau). m1 / minerbo / kershaw
  keep `upwind` unless the input names `plm` (not gated there; see below).
  - The key is now always read and echoed.
  - d0c59f7c read it only when named, so its restart files usually lack it. On a restart
    (`global_variable::restart_run`, set in main() for `-r`) a missing key defaults to
    `upwind`, so an old restart keeps its old scheme.
- `<rad_m1>/implicit_halo_mpi` defaults to true when implicit_halo_direct is on and the
  mesh is single-level without cubed-sphere seams or polar boundaries; otherwise false
  (the ordinary exchange), silently. It is bitwise the ordinary exchange.
- Explicit input values override both. implicit_vimp, time_scheme, implicit_krylov_pipe
  unchanged.

## Gates (CPU, gcc + openmpi, box_convection)

Inputs: the runs_3w_krylov `slab2d_def` / `box3d_def` as `_old` (names `implicit_enthalpy
= upwind`, `implicit_halo_mpi = false`) and `_nd` (names neither). 2-D slab 200 s; 3-D box
84x32x32, 4 blocks, 60 cycles. bin/rst compared from `<par_end>` on (the header only
echoes the input). Gates 1-5 were run on a84455ba.

| # | test | result |
|---|---|---|
| 1 | explicit old values, ref d0c59f7c vs new: slab Edd 1/2/4 ranks, slab vet_sc 2 ranks, box 1/4 ranks | bitwise (hst, bin) |
| 2a | new, no keys vs ref with explicit plm (halo_mpi off): slab 1/2/4 ranks, vet_sc 2 ranks, box 4 ranks | bitwise (hst, bin, rst payload); halo_mpi ON in the new 2- and 4-rank runs |
| 2b | solver health, new defaults | NON-CONVERGED 0, positivity fallbacks 0, line_jacobi fallbacks 0 |
| 2c | plm vs upwind, relative, last hst row | slab Edd totE 2.2e-6, KE1 3.6e-3, KE2 5.0e-3; vet_sc totE 1.3e-6, KE1 2.1e-3; box (t = 9.7) KE1 2.0e-5, KE2 1.6e-3 |
| 2d | Picard mean, upwind -> plm | slab Edd 2.268 -> 2.048; vet_sc (2 ranks) 2.686 -> 2.623; box 2.233 -> 2.233 (4 ranks) / 2.217 (1, 2 ranks) |
| 2e | box new defaults, 1 vs 2 vs 4 ranks | differs as upwind already does: KE, totE, dt within 1e-9; only the near-zero 2-/3-momentum sums move at large relative size |
| 3 | restart, new defaults, 2 ranks, t = 40 full vs restart at t = 20 | bitwise (final bin, final rst, 20 hst rows) |
| 4 | d0c59f7c restart file (no enthalpy key) continued by the new binary vs d0c59f7c full run | bitwise; the continued run is upwind, halo_mpi off, and echoes `implicit_enthalpy = upwind` |
| 5 | deep Hot Jupiter hydro, cs_hyd4_prod dhj.00117.rst (t = 1.78426e7 s, cycle 905665, rotation 58.5), 3 cycles, 8 ranks | bitwise, all 6 output files |

Re-checks of a5220da7 (closure-limited default):

| test | result |
|---|---|
| slab Edd, 1 rank: no keys (a5220da7) vs ref d0c59f7c with explicit plm | bitwise (hst, bin, rst payload); NON-CONVERGED 0, positivity fallbacks 0, Picard 2.048 |
| slab `closure = m1`, 1 rank: no keys (a5220da7) vs ref d0c59f7c (upwind) | bitwise; both NON-CONVERGED 6, positivity fallbacks 1448, Picard 11.98 |
| restart echo, no keys | Eddington: `implicit_enthalpy = plm`, `implicit_halo_mpi = 1`; m1: `implicit_enthalpy = upwind`, `implicit_halo_mpi = 0` (halo_direct is off for m1) |
| cpplint rad_m1_implicit.cpp | only the pre-existing build/include_subdir errors |

## m1-closure finding (pre-existing)

The seeded 2-D slab with `closure = m1` and the OLD default (upwind) already gives
NON-CONVERGED 6 and 1448 positivity fallbacks (Picard mean 11.98) over 200 s; plm gives
NON-CONVERGED 7 and 1421 (Picard mean 12.20). This is a pre-existing m1-closure problem,
not caused by plm; it is why the plm default is limited to eddington / vet_sc / tau.
