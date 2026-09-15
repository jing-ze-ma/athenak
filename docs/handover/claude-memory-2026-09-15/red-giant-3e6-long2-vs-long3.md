---
name: red-giant-3e6-long2-vs-long3
description: 2026-09-13 comparison at 3e6 of the two complete red-giant runs (bench/RG_fofc_long3/analysis/VERDICT.md): long3 (wb_rmin 6e11) base convection SATURATED (rms/v* 2.0 at i=5, e-fold at i=10 3.3e6 = flat over 2.5-3e6, zone extending outward) while long2 (WB everywhere) still exponential; long3 deep KE_r 3x, i 40-99 113x; mean rho/e at i 0-40 unchanged <0.4% since 9.2e5 (the KE grows on an unchanged background); floors/fofc/tclamp per cycle HALVED in long3; FOFC fires EXCLUSIVELY in the corona (i>=428, r>3.83e12, 100% at i 445-479 in both dumps, seam/vertex enrichment 0.99/0.82 = none); corona drained 3x further in long3 (top shell 2x dfloor); tot-E drift 8.8e-5 in long3 (front-loaded +5313 L -> +418 L) vs 1e-6 in long2. VERDICT: long3 config = production basis, but first attribute the E drift and settle the top boundary (sponge or corona floor)
metadata:
  type: project
---
Tables in cmp23.txt, fofcloc.txt. eos_fail 0, efloor_de 4e-38 of E, mass exact in both. dt long3 30.6 -> 15.8 s
(saturated convection at the base binds dt), long2 -> 26.5. The 4 fofc/tclamp episodes in long3 (1.04-1.13e6,
1.57-1.73e6, 1.94-2.17e6, 2.48-3.0e6) are corona events. Next: (a) E-drift attribution (corona/top boundary first,
then the i 0-2 wall mean flow, then the T clamp), (b) the top boundary: sponge vs a corona density floor.
**E DRIFT ATTRIBUTED + FIXED (2026-09-13 ~12:30, EDRIFT.md):** the wall_noflux correction returned dm*h with the
CELL-CENTRE potential while the face flux removed it with Phi(face)=0 -> |dm|*Phi(i=is) (3.9% of h, always heating)
left in the wall cell every stage; residual = 0.40-0.47 x |dM|*Phi in every bin of long3 AND long2 (the counter's
stage-1 bias); wb_rmin was the TRIGGER (wall mass flux 90x larger once WB left the wall cell), not the source. FOFC,
tclamp, dfloor, the WB source, the top face and radiation all ruled out (code + budget). Fix committed (h from the
gas total energy + face potential), unverified by a run (RG runs stopped by the user). Not the cause of the deep
mode (long2 drifted 1e-6 with the same growth rate) but plausibly of long3's 3x more vigorous base convection.
**Wall fix, take 2 (d2788159):** the gas-enthalpy + face-potential form (05130338) drifted tot-E by -6.8e-5 in 2e5 s on
the no-WB run RG_v4 (wall mass-flux counter ~1e6 L: without WB the wall cell is far from balance, so any mismatch in the
cancellation is amplified 100x). Any h rebuilt from the cell state fails; the correction now subtracts the ACTUAL
hydro energy flux bdt*A*F_IEN/V at the face (exact). RG_v4 relaunched from scratch on it (11670700/1); gate: E drift
< 1e-5 by 3e5 (first attempt kept in RG_v4/old_11670310).
