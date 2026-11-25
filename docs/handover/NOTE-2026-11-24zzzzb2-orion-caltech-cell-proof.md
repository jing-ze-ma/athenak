# NOTE 2026-11-24zzzzb2 orion -> Caltech: t9l62 one-cell proof (USER approved)

Branch to build: `caltech-sdocell-883d79d97`. Priority: idle node, after the HG repeat (24zzzzb1) has its nodes; production first. Text by SDORDER:


**Question.** In S5G one cell paid 9.917e32 erg (on arm, cycle 5346) and 1.242e32 (off arm, cycle 5312): a paid raise in a floor_audit [refresh]
with e_c <= 0 (+1), KE-dominated (+1), r ~ 2.47e12 cm in both. Which state of the step has the negative eint first?

## Binary
Branch `vet-sdorder-cell` **883d79d97** = 95813cf85 (= your 30d01765's source) + one diagnostic commit (3 files, +115 lines, not for merge).
It only reads: env `SDO_CELL` unset = bitwise 95813cf85 on Orion CPU (sp + cs shells, key off and on, 24 dumps / histories); env set =
bitwise too. **Not compiled with CUDA on Orion.** If nvcc rejects SdoCell / SdoEint, send the first error, don't edit. PROBLEM = he_star_vet,
your CUDA 12.9 + g++ 13 recipe, separate dir, tiny gate without vet_gd_ls = 3.686e-11, .prov checked. One binary for both arms.

## Arms (one node, 4 H200, the same node layout and ranks as S5G: 4 ranks; nothing on production)
Start: /resnick/groups/carnegie_poc/jingze/compton_1124/rst/t9l62.00015.rst. Args: the S5Gon / S5Goff ARGS line of SDO_1124/S5G{on,off}/run.log,
word for word, with these changes only:
| arm | key | nlim | added |
|---|---|---|---|
| Con | `rad_vet/time2_gtime_slope=true` | `time/nlim=5350` | `hydro/floor_audit_worst=true`, env `SDO_CELL=0` |
| Coff | key not named | `time/nlim=5316` | the same |
Run Con first; run the same on a second node or after. `job/final_restart=false` kept (restart: absent key is added, NOTE line).
Outputs: no bins, no rst. hst as in S5G. ~10 min each.

**Precondition (say it first):** CLEAN FALLBACK floor_pay of Con at cycle 5346 = 9.917e32 and of Coff at 5312 = 1.242e32, to the printed digits
(each run is deterministic across your S5 / S5F / S5G binaries: off 1.242e32 in all three). If not, the events moved: report the list of cycles with
floor_pay > 1e30 and stop.

## Return
- log.0.gz, user.hst, run.log with ARGS; `grep -aE '^SDO(CELL|BOX)|floor_audit|CLEAN FALLBACK'` lines for cycles 5344-5348 (Con) and 5310-5314
  (Coff) as text; CW / fallback / BE / retry counts.
- Where the cell is: `SDOCELL <point> ... gid .. kji .. r .. th .. ph` (the rank's cell of smallest eint at that point) plus
  `floor_audit_worst ... [refresh] rank gid kji r th ph dE rho E eint_pre eint KE123` for every c2p. The two must name the same gid / kji / r (2.47-2.49e12).
  The 5^3 box (rho, M, E, eint, E_rad, offsets dk dj di) is printed as SDOBOX lines under each SDOCELL; units code.
- Points, in step order, for each rank (a point prints only where its smallest eint <= 0):
  1. `stage1-Y` u0 entering the fluid task after the stage-1 solve (Y = S + d);
  2. `which2` the same cell of the smallest eint(U* + d), with **U^n, U*, S, Y, U* + d** (5 states: rho, M, E, eint, E_rad) (key code: rows 0-5);
  3. `which3-entry` the register after FloorGasCap and the first refresh's floors;
  4. `stage2-postsolve` u0 after the stage-2 solve, before its refresh.
  A cycle-0-type BE step prints only point 4.

## Criteria, fixed before the run (per event, both arms)
Let e_X = eint of state X in the event cell, dM = M(U* + d) - M(U*), dE = E(U* + d) - E(U*), v_X = M_X / rho.
- **PROVED (U* + d cross-term):** at cycle 5346 (Con) / 5312 (Coff) the first point with e <= 0 is `which2` with e(U* + d) <= 0 in the paid cell, while
  e(U^n), e(U*), e(S), e(Y) are all > 0; and |e(U* + d) - [e(Y) + e(U*) - e(S) - (v_U* - v_S).dM]| <= 1e-6 x max of the terms (identity for constant rho;
  if rho differs by more than 1e-6 between the states, say so and skip this check). The paid raise must be the paid cell's (same gid kji).
- **REFUTED:** any of: e(U^n) <= 0 or e(U*) <= 0 (the hydro / register was already bad), e(Y) <= 0 at `stage1-Y` (the solve's write-back), or `which2`
  prints nothing and `stage2-postsolve` is the first with e <= 0 (stage-2 increment or write-back), or the event cell is not the SDOCELL cell.
- **INCONCLUSIVE:** precondition fails; or no SDOCELL prints at the event cycle though floor_pay does (then the negative state is made inside a c2p's own
  refresh: say so with the floor_audit_worst line).
Both arms are scored by the same rules; a different verdict on Con and Coff is a result, not a failure.
Do not continue past nlim, do not change any key. Nothing else runs on the nodes.

-- orion (coordinator), drafted by SDORDER
