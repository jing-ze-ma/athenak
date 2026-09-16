---
name: red-giant-floor-energy-creation-fix
description: 2026-09-10 02:30 — the energy-creating floor FOUND and FIXED (src_floor, bin_floor, floorfix.diff): (1) on the cubed sphere the floor fires in Coordinates::GnomonicEquiangleRaiseVel (defer_cons_floors), which rebuilt E = e_floor + KE keeping KE when eint<0 -> donates |eint| (91% of the +1.7e39 erg was ONE floored cell's KE); (2) the tabulated EOS pins the T solve at 10^(logt_min-3) = 0.1 K, so EnergyFromPressure(pfloor) returns e(0.1 K) = 1.42264e7 rho, a pressure 4.7e9x pfloor. Switches hydro/efloor_from_ekin and hydro/eos_floor_consistent (default off, control bit-identical); FL1 with both on passes the 1.43114e6 death; efloor_de event-log column
metadata:
  type: project
---

Agent a87c9fc892f825a2c (session 293c860c). Probe: bin_floor/eosprobe2.cpp; diff (637
lines, +262/-30: eos.cpp/.hpp, eos_table.cpp/.hpp, general_c2p_hyd.hpp, general_hyd.cpp,
coordinates.cpp, prolong_prims.cpp, mesh.hpp, eventlog.cpp) at bin_floor/floorfix.diff;
reviewed by me: the RaiseVel hunk is a parallel_reduce that, when eint < e_floor and
efloor_from_ekin, rescales m1..m3 by sqrt((E-e_floor)/KE) (writes them back), keeps E,
creates energy only if E < e_floor, and sums the created energy into ecounter.efloor_de.
EOS facts: the table (logt 2..7.5) is NOT clamped at 100 K (log-log extrapolation, p/e
consistent); the CLAMP is the root-find bracket 10^(logt_min-3)=0.1 K: below e/rho =
1.42264e7 T=0.1 K, p=4.71e-3 (rho/1.35e-9), cs=2100.9, Gamma1 1.265 independent of e (that
constant e/rho=1.423e7 is the fingerprint of every "floored" cell in the censuses).
eos_floor_consistent: for e < e_min(rho)=e(rho,T_min): T,p scale linearly with e, Gamma1 at
T_min, exact inverse in EnergyFromPressure (pfloor means pfloor); hot-path cost one divide.
Not done: EOS_Data::ThermoAt (WB background walk) has no sub-floor branch.
Tests: FL0_ctrl (195283, switches off) dies 1.43114e6 at (0,9,9,295) T 1.766e7 = D1 exactly,
dt history identical at all 14 diagnostic cycles; FL1_floorfix (195284, both on) NO death,
0 collapses, past 1.4539e6 at 02:30, efloor_de per log interval 3e18 -> 6e8 (4.9e9x less),
vertex max|v_r| 4.2e5 (was 9.95e7), vertex-column E smooth. Which switch carries it is NOT
separated (mechanism says efloor_from_ekin); FL2 (ekin only) requested.
OPEN: the MOMENTUM anomaly at the floored vertex cell (v_r -3e5 -> -2.7e7 in one cycle,
2.6e4x more force than grad p or rho g) is capped by the fix, not explained: suspect the
cubed-sphere geometric/WB source term at a floored cell.
Related: [[red-giant-vertex-floor-hydro-trigger]], [[red-giant-vpert-seed-chart-imprint]].

## RESOLVED 2026-09-10 ~02:10: efloor_from_ekin carries it ENTIRELY
FL2_ekin_only (195292) survives (tracks FL1 digit for digit, efloor_de 6.31e8); FL3_consistent_only
(195293) dies at 1.43089e6 (tracks FL0, efloor_de 3.14e18). eos_floor_consistent is a
correctness fix (pfloor really 1e-12) but not the cure. Ship efloor_from_ekin=true.
## MOMENTUM KICK EXPLAINED: the WELL-BALANCED background (vtxcheck_floor.py)
wb_rmax is ABSENT from rg.athinput -> WB active at 1.082 R (wellbalance_dynamic, polytropic,
cache every 10). Budget of the vertex cell 301 at the last pre-event dump (dyn/cm^3): observed
d(rho v_r)/dt -2.1e-4; -rho g -2.1e-9; grad p +4.0e-9 (outward, the cell is a p MINIMUM);
curvature +3.7e-11; WB source -2.06e-6 (973x gravity, inward) = the only anomalous term. The
hydrostatic half-cell walk p_l = p exp(rho g dr/(2 p chi_rho)) has exponent 9.86 for this
74.8 K cell (0.15-0.35 for normal cells) -> p_l = 6406 vs the cell's p = 0.333 (19219x); the
log gain dln p_l/dln p = -10.9, so a 1.5x drop of p between cache build and trigger gives the
missing 100x; after the floor clamps p (0.333 -> 4e-4) the exponent is 8063, exp overflows,
WBGuard flattens the background and the cell gets NO gravity. The shared 300/301 face gets
p_L ~ 23 vs p_R ~ 6406 in the background alone. => a cold floored cell under WB receives a
spurious inward force ~1e3 g. Same family as [[red-giant-wb-kills-ambient-medium]].
NEXT TEST (one input line): <hydro>/wb_rmax = 3.3e12 (just inside the star) on the FL0 config
without the floor fix -> does the 1.43114e6 death vanish? and with the fix.
CONFIRMED 01:37: FL4_wbrmax (wb_rmax=3.3e12, floor fix OFF) passes 1.43114e6 with no collapse
(at 1.4413e6, dt 14.3); FL5 (wb_rmax + efloor_from_ekin) identical dt. => the WB background
at a cold cell is THE trigger of the open-top deaths; wb_rmax inside the star cures it
with no code change; efloor_from_ekin is the second, independent safety. FL1 at 1.547e6,
FL2 at 1.487e6, both clean. Production inputs MUST set hydro/wb_rmax (e.g. 3.3e12 for x1max
3.52e12); check the dhj inputs too ([[red-giant-wb-kills-ambient-medium]]).
FINAL 02:20: FL1/FL2/FL4/FL5 all reached tlim 1.56e6, ZERO collapses. wb_rmax=3.3e12 (FL4)
is THE cure: efloor_de unchanged (3.0e18 = control) yet no death, and the cold-cell census
HEALS 192 -> 8 cells by 1.56e6 (coldest never < 2e-2 of median), whereas the floor-fix-only
runs keep 150-300 cold cells (FL2's worst 8e-6). FL5 (both) = FL4 bit-for-bit in dynamics
(dt 5th digit, census 4 digits), efloor_de 5.5e8. wb_rmax runs also end with higher dt (15.6
vs 15.0). => the WB background above the star does not just kick floored cells, it SUSTAINS
the cold-cell population (it fights every cold cell with a spurious force). prod12: set BOTH
(orthogonal, FL5 costs nothing). The knee-cooling story remains the source of the cold
cells' first drift; wb_rmax removes their amplification.

## wb_rmax A/B DONE 2026-09-10 ~02:40 (FL4/FL5, jobs 195300/195301, verified by me)
FL4 = wb_rmax 3.3e12 ALONE (both floor switches off): NO death, 0 collapses, reaches tlim
1.56e6 (cycle 69528), dt tracks FL2 within 1%, while efloor_de stays at the FL0/FL3 level
(2.9-3.1e18 per interval, 2.4e21 total = the energy-creating floor fully active but harmless).
FL5 = wb_rmax + efloor_from_ekin: identical trajectory to FL4 to 5 digits in dt, efloor_de 6e8.
=> the WB background source term at a cold floored cell is what turns the floor event into
the conduction-dt collapse; either remedy alone cures 1.43114e6. Ship BOTH in prod12:
hydro/wb_rmax=3.3e12 (WB has no business above the star anyway) + efloor_from_ekin=true.
Ported into the main tree 02:00 (build_port; verified byte-identical to src_floor + a ThermoAt
sub-floor hunk); random symmetry-breaking seed added (problem/vpert_rand_amp etc.); smoke
tests V5_portcheck (bit-identity to V4 at amp=0) and V6_randseed running.
