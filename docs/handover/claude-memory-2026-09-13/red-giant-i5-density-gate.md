---
name: red-giant-i5-density-gate
description: "I5 DENSITY-GATED OPACITY (09-11 01:00): I3's corona config with the radius cutoff replaced by hydro/rad_gate_rho=1e-16 (dex 0.5), from V9f rst 5e5: SURVIVES the 5.7e5 ejection (5 transient hydro dips, dt 30.63 flat from 6.23e5), tot-E flat, and the post-burst emergent luminosity SETTLES AT 1.15 L (I2/I3 with the cutoff: 2.7-4 L; lidded: 1.3-1.5 L) with L_cut back to 0.06 -> the cutoff artefact is CONFIRMED and the open-top surface radiates ~1.15 L. Cautions: L_out 0.35 L right after the restart (unexplained), 15k two-stream RESCUE events in the gated corona base (i 415-422) during the burst, ~1e3 per 1e5 s after"
metadata:
  type: project
---
Job 195989, dir I5_gate, binary build_impl/athena_pin3 (gate code in conduction.cpp ~122+,
shared RadGate at every kappa site incl. two_stream_rt.hpp; first launch 195980 died on the
new gate+rmax fatal because the RESTART EMBEDS rad_kappa_rmax -- set it to 0.0 explicitly).
Threshold 1e-16 is a decade BELOW the corona base (1e-15): the corona base is radiatively
active, which is why RESCUE fires there; the agent was asked to justify/adjust.
L_out/L timeline: 0.35 (5.02e5) rising to 0.61 (5.51e5), 1.32 (5.64e5), burst values 5e2..4.5e10
(5.99-6.14e5, single hot cells, hst untouched), 519 (6.23e5), then 1.14-1.17 (7.17-7.29e5).
Cap burst 5.8-5.9e5 only. Gate comparison runs: J3_cut/J3_gate/J3_g1e30/J3_nocut (job 195979,
two steps failed on the embedded-rmax fatal). See [[red-giant-open-top-kappa-cutoff-artefact]],
[[red-giant-i4-active-medium]], [[red-giant-implicit-radial-diffusion]].

## FINAL 09-11 01:40: I5 ran to tlim 9e5 (job 195989 COMPLETED). Corrections to the note above:
- The 0.35 L at restart is NOT the gate: I5 == I3 (same opac_tmin=-1) at 5.208e5 to 7 digits
  (0.3975326 vs 0.3975325). CLAMP OFF lowers L_out right after the restart (0.40 vs I2's
  0.62) -- the grain/molecular opacity makes the cool top layer thicker until it readjusts.
- The RESCUE cells (i 415-455, rho 2e-14..1e-13, G=1) are EJECTED STELLAR GAS next to runaway
  cells, not the corona; the corona (r>3.85e12, rho 1.2-1.6e-17, T 2e4) has median G 3e-4 to
  8e-4 and is inert as designed. The corona's real density is 8e-17 -> 1.5e-17 (bg_rho=1e-15
  only LOCATES the join; rho_c = 8.1e-17 pressure-matched). Threshold 1e-16 justified: 99.4%
  of corona cells G<0.1, 95.9% of atmosphere cells G>0.9; 1e-14 would inert 49.8% of the
  atmosphere (worse than the radius cutoff's 39.3%).
- L_out/L I5 vs I3/I2: 1.19 vs 2.55/2.70 (6.95e5), 1.17-1.20 vs 3.12/3.17 (7.3e5), 3.32 vs
  5.35/5.08 (8.2e5), 2.36 vs 3.87/3.70 (8.8e5), 1.98 at 9e5; L_cut 0.05-0.10 vs 1.36. So the
  cutoff artefact is confirmed, but the open-top surface still fluctuates 1.2-3.3 L (lidded
  1.3-1.5). Mass above 3.6e12 x65 (5.2->6.0e5), same as I3. Photosphere above the old cutoff:
  0.0000 at every dump (so the "photosphere crossed the cutoff" wording was wrong; it is the
  tau=100 CUT and the ejected gas that crossed it).
- NEW, from the hydro-cell collapse report: photospheric cells reach T = 1.4e16-1.3e17 K
  (EOS table max 3e7!) with dt only dipping to 3-6 s. tot-E 6-digit flat cannot see 1e41 erg.
  Conduction can no longer create a hotter-than-neighbour cell (implicit + cap), so the
  source is HYDRO (vacuum-Riemann / shock into an evacuated cell) or the RT source
  (semi-implicit relaxes toward a neighbour-set equilibrium; rescue clamps). "Stability was
  bought; the runaway itself is untouched." NEXT: per-task max-T scan (RGNanScan style) on a
  restart from 5.0e5 to catch the first operator that lifts a cell above its neighbours.
- Gate sites: 11 (conduction tau, x1 explicit, cap c2/c3, x2/x3, implicit, NewTimeStep x2,
  two-stream grey + ck); conduction gates the OPERATOR (flux/stiffness x G), the opacity
  sites use G*kappa + (1-G)*kappa_above; ck path gained a gate it never had.
- Pins: athena_pin (implicit only) 2a6712ed, pin2 (+cap) 5d9be6a1, pin3 (+gate, hydro cell
  report) 16e8203d; build_impl/src/athena == pin3 == working tree.
