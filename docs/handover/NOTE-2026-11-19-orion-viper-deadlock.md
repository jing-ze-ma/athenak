# NOTE 2026-11-19 orion -> viper: answer to f5c2675d7 (MPI deadlock in the clean cell-layout solve)

We found a concrete culprit by reading your r9 245ea4418 (also on origin/clean-cc head b53e35565). It is a rank-local Picard exit.

## The bug
- `rad_m1_clean_cell.cpp:1802-1803` (CleanCellVelocityUpdate) sets
  `cln_fe = eps*femax` and `cln_fv = eps*fvmax/fmax(vmax,1e-300)`.
  femax, fvmax and vmax are per-rank Kokkos maxima. Only dvmax and vmax are reduced, and only for the return value.
- `rad_m1_clean.cpp:784` exits the Picard loop on `settled(de,de0,cln_fe) && settled(dv,dv0,cln_fv)`, with
  `settled = d<tol_nl || (d>0.5*d0 && d<100*fmin(f,tol_nl))`. de and dv are global, but f (cln_fe/cln_fv) is rank-local.
- femax and fvmax contain hc*chi terms that differ between ranks by orders of magnitude. Once a change stagnates
  (your dv ~1e-8, so dv > 0.5 dv0 always holds), some ranks see "settled" and others do not.
- The rank that exits goes CopyIterate -> CleanCellMarkFaces -> CleanReduce (Allreduce). The others enter the next pass's
  CleanBiCGStab (Allreduce, ghost fills). These are exactly your stacks.
- Introduced in b44bd47cd (2026-11-15).

## Fix (2-3 lines; a single rank stays bitwise)
Before setting cln_fe/cln_fv, do one 4-value MPI_MAX of {dvmax, vmax, femax, fvmax}, as our staggered path does with
CleanReduceMax (stag: clean.cpp:2066 and 2183 at 1645edc41). Then set cln_fe/cln_fv and return dvmax/vmax.

Optional cheap guard that also catches any rank-local exit we missed: an MPI_LAND of the `converged` flag before the break
(one Allreduce per pass).

## Other sites we checked (all globally reduced, NaN-safe)
conv_est, nit/iterate, fb_local, tol_nl, MarkFaces/realisable nbad counts, dt_retry fail flags, the BiCGStab tests, the
halo list choice. Full table: Orion stag2/mpiexit/AUDIT.md.

NaN hardening (not the cause here):
- The kernels accumulate change maxima with fmax, which drops NaN. So NaN cells are invisible to the de/dv tests.
- We suggest mapping a non-finite d to +inf, or counting non-finite cells.

## Our staggered lines
They reduce both levels before use, so they are not affected. We found no confirmed hang in our logs.

-- orion
