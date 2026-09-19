# Handover 2026-09-19 B (viper, written ~23:45)

Supersedes the He4 part of HANDOVER-2026-09-19.md. Everything is on branch he4-presn-global,
worktree bench/wt_he4, HEAD = "tests_r13: fitgrid.py coarse-atmosphere target ..." on top of
341e6b4b, NOT pushed. Full chronology with 7 retractions: memory note
he4-r11-floor-arms-verdict; the agreed next design: memory note
he4-rad-eos-force-split-design (both copied to docs/handover/claude-memory-2026-09-19-B/).

## Fixed today (all default off unless noted, default-off builds bitwise)
- problem/mlt_split_deposit: MLT flux divergence deposited in the Strang half steps BEFORE the
  two-stream sweep (was inside the RK stages: O(dt) force dipole under the Fe bump, 4.3 -> 0.4 %
  of rho g).  27e0c387.
- problem/cs_max: sound-speed ceiling as a density floor (hot evacuated cells kept LTE radiation
  pressure: c_s 4e9 cm/s, dt 17 -> 0.3 s).  problem/e_ledger, problem/rt_no_heat diagnostics.
- REAL BUG, always on: the transverse P_rad grad w force read an UNFILLED ghost temperature at
  every MeshBlock edge (T_g filled over active angular cells only).  de1d41e4.  It was a 3-D dt
  killer on the cubed sphere too; earlier 3-D "convection onset" and "hot voids" are suspect.
- Implicit transverse radiative diffusion (ADI) admitted on a spherical-polar WEDGE off the axis,
  periodic in theta and phi.  47f80807.  New input inputs/hydro/he4_presn_sp.athinput.
- problem/vpert_sp_rand (periodic random entropy seed on the wedge, a61b29a8),
  problem/vdamp_all_mean_only (damp only the shell-mean v_r, 341e6b4b).

## Established physics / numerics
- The expansion is THERMALLY driven (adiabatic arm quiet); user: treat the expansion/pulsation
  as physical.  Energy ledger: the star gains 0.2-0.45 % L, all through the two-stream sweep.
- Wedge 90x90 deg, 96x192x192 (wedge2/3): clean 3-D keeps dt 3.2 s through 2 turnovers.
  Entropy seed: convection's horizontal KE grows x~3 per 0.25 turnover, independent of the MLT
  closure; withdrawing the closure (wedge4) only feeds the radial expansion.
- The radial motion IS a pulsation but violent: outward to ~2.75 turnovers, photospheric gas
  reaches the top at 1.014 R, 17 % of the domain mass leaves, then fallback at Mach ~0.3.
  No boundary trick at 1.014 R conserves the envelope (tests_r13 b0 / b_in / b_lid).
- Tall domain (x1max 1.25 R, stretched grid, floor atmosphere rho 3e-13; tests_r13 b_tall,
  b_tall2): conserves mass, shows a radiatively driven OUTFLOW above the photosphere, but the
  EOS temperature gate radiation-loads shock-heated transparent gas (dt 0.6 s at cs_max 1e8).

## Running at handover
- tests_r11/wedge3 (job 11850864, 4 nodes): fallback phase, dt ~0.6 s; stress test only.
- tests_r11/wedge5 (job 11852628, 4 nodes, athena_v19): wedge3 + mean-only radial damping until
  4 turnovers (18820 s), then release, to 6.  READ THIS FIRST: does KEh saturate, does the
  closure let go, is the released pulse milder?
- tests_r13/b_tall2 (apudev chain 11852759-61): tall domain, nx1 112, cs_max 3e7, to 4 turnovers.

## Next (agreed with the user)
Radiation redesign, see he4-rad-eos-force-split-design: EOS weight density-only with a very low
window (1e-11..1e-10, T gate off) + optical-depth-gated FORCE in two_stream_rt.hpp (default-off
switch, s=0 must reproduce today's formula bitwise) + IC eint regenerated with the new w above
~0.97 R; validate at t=0 with tests_r12/fbud.py on the 1-D sp column (tests_r13/r13_arm.sh), then
the tall domain, then 3-D.  Opus agent quota is back: delegate the implementation, verify diffs.
Parked: port 4bcdc855 + 94c7165d + de1d41e4 to rt-integration, bvals_fc seam defects, pushes.
