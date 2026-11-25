# NOTE 2026-11-24zzzzd2 orion -> Caltech: t9l62, poloidal momentum removed by floor_from_ekin (USER approved the counter)

**Changes to the draft below:** `time/nlim=5350` (71 cycles, so the 5346 payment is inside), branch name `caltech-sdomom-373e2cc1c`, scorer in docs/handover/caltech_tools/mscore.py. **Priority: after zzzzb5, zzzzc6 and zzzzd1**; one arm, one node, ~5 min.

**Question.** With `hydro/floor_from_ekin` a cell below the temperature floor takes the deficit from its poloidal kinetic energy (E, L_z, mass
fixed). How much poloidal momentum does that remove per cycle, where in the step, and how does it compare with the losses already booked
(floorcap, dfloor_flux_limit)? In your SDO_CELL logs (printed cells only) 55 of 57 cells with e(U* + d) <= 0 were filled this way.

## Binary
Branch `vet-sdorder-mom` **373e2cc1c** = a000d803c + one diagnostic commit (4 files, not for merge, not for the release). Env `SDO_MOM` unset:
bitwise a000d803c on Orion CPU (sp with floor_from_ekin with and without vceil, cs; slope key off and on; 6 files per pair, 6 pairs; with
SDO_MOM=1 also bitwise, sp, 4 pairs). Counts equal the floor_audit skip counters exactly (4 sp arms); 2 ranks = 1 rank (145 and 73 lines, 0 differ).
**Not compiled with CUDA on Orion**: if nvcc rejects SdoCause / MomPost / MomReport, send the first error, don't edit. PROBLEM = he_star_vet, your
CUDA 12.9 + g++ 13 recipe, separate dir, tiny gate 3.686e-11, .prov checked. Spherical polar only (SDO_MOM is FATAL on the cubed sphere).
`SDO_MOM` needs `hydro/floor_audit` or floor_pay on (both are in the t9l62 ARGS).

## Run (one node, 4 H200, 4 ranks, nothing on production)
Start: /resnick/groups/carnegie_poc/jingze/compton_1124/rst/t9l62.00015.rst. Args: SDO_CELL/Con/run.log ARGS word for word, with only
`time/nlim=5350` (the restart starts at 5279: 71 cycles) and env `SDO_MOM=1` (no SDO_CELL; `floor_audit_worst` as in Con).
**Not a precondition, but report it:** hst rows vs your Con run; a000d803c adds two clean_stag_thinf commits, so a difference is possible.
The 9.9e32 payment (5346) and nothing after 5329 is covered; say if you want 71 cycles instead (nlim 5350, same cost).
Outputs: no bins, no rst. ~5 min.

## Return
log.0.gz, user.hst, run.log with ARGS; the stdout of `sdorder/mscore.py log.0.gz` (below); CW / fallback / BE / retry counts.

## What the lines mean
`SDOMOM total cycle`: |p_pol|, signed p_r, p_th, KE_pol of the whole mesh at the step's end. `SDOMOM sum cycle <tag> <cause>`: cells, signed and
absolute dp_r, dp_th (volume integrated), KE removed, eint gained, dL_z, |p_pol| of those cells before, dE. Cause: dfloor (rho changed), vceil,
ekin (E fixed), other. Tag: update1 / update2 (after each Heun hydro stage), which2 (register U* + d), post2 (after the stage-2 solve), flux3r (flux
state, register restored afterwards: not persistent), flux3 (key off), other. `SDOMOM worst` = per rank, per conversion, the ekin cell of largest
|dp_pol| V (tag 0 update1, 1 update2, 2 which2, 3 post2, 4 flux3, 5 flux3r, 6 other). A loss of update1, which2, flux3 enters the final state with about
weight 0.5 (inferred from the Heun average).

## Reading, fixed before the run (no tuned number)
- dE = 0, dL_z = 0 in every ekin line, and dKE = deint to round-off: identity check; a failure means the counter or the claim is wrong.
- **Comparable to the booked losses:** weighted ekin |dp_r| + |dp_th| per cycle <= floorcap |dp| per cycle + dfloor_flux_limit x1 momentum per cycle of the same
  run. **Larger:** above it. Both are reported with the share of |p_pol| (no threshold).
- Where it happens: the share of which2 (made by the scheme's register increment) against update / post2 (hydro and solve).
