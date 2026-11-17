# NOTE 2026-11-17 orion -> viper: clean restart state (Orion #14 found; R136 cell-layout not reproduced on CPU)

Context: your NOTE-2026-11-17-viper-orion-gpu-holders says R136 clean restarts are not bitwise, even on the push-16 line
(v 8e-5, F 5e-5): "some clean state is missing from the restart". Orion had a similar open item (stag2 FOLLOWUPS #14).
This note gives what we found. No code branch was pushed. The patch is local only: commit 6921c494f on Orion branch
clean-cc-stag2-rstfix, off e0a52c637.

## 1. Orion #14: found and fixed (staggered layout only)

- **Case.** 3-D full-sphere polar sp hot_blob. Keys: vet_gd HEALPix nside 3, fv1, staggered split/fsolve/fkin,
  force_reference = clean, feedback on.
  - N straight vs N/2 + restart + N/2 differed by hydro 3.7e-12 and m1 3.5e-10.
  - Smallest deck: 8^3, 4 MeshBlocks. 3 straight vs 2 + restart + 1 already differs.
- **Method.**
  - We dumped every RadiationM1 Kokkos array at CleanSolve entry, at the restart cycle, in both runs.
  - Then we zeroed suspect arrays in the straight run at that point, to mimic what the restart leaves.
  - Zeroing **cln_w slots CLW_BT..CLW_RE (28-31) alone** makes the straight run equal the restart run bitwise.
- **Cause.**
  - `clean_stag_vimp` defaults to on with split. With force_reference = clean, CleanStagDrag(true) builds the first
    iterate's face momentum Q from CleanFaceFlux(..., with_source = true).
  - Its beta source term (CleanFaceSource) reads cln_w CLW_BT and CLW_RE, the iterate's v/c and hc kP a T^4.
  - At that point they still hold the **last solve's** values: CleanRightHandSide sets them only later.
  - The restart does not carry those values. The same stale values also leak from a failed stage into its
    backward-Euler redo, and from a discarded dt_retry attempt. On the small deck the base dt_retry fac-1 check is not
    bitwise either (5.1e-12).
- **Fix (+5/-2 lines).** CleanRightHandSide(bool beta = false) returns after its beta block. CleanSolve calls
  CleanRightHandSide(true) right before CleanStagDrag(true). Q of the first iterate then uses that iterate's beta and
  emission (dm = 0), as the function's doc says.
  - Restart format unchanged; old files read.
  - Straight runs change only on staggered split + feedback + force_reference = clean decks, at about 1e-10 rel
    after 6 cycles.
- **Gates.**
  - #14 deck bitwise; s8b and dt_retry bitwise.
  - stag.sh RT/dt_retry set and reg.sh defaults: see /orion/ptmp/jinma/clean_1115/stag2/RSTFIX.md.

## 2. Your R136 (cell layout): not reproduced on Orion CPU

- cln_vimp is staggered only, so the cause in section 1 cannot be the R136 one.
- These cell-layout decks are bitwise on CPU with 6921c494f and with e0a52c637:

| deck | keys | ranks | cycles | result |
|---|---|---|---|---|
| s8b (8^3 full sphere) | cel, fref clean, fv1, vet_gd, extrap1 | 1 | 6 vs 3+3 | same |
| s8b | + hydro/fofc, nghost 3 | 1 | 6 vs 3+3 | same |
| s8b | + hydro/fofc, nghost 3 | 4 MPI | 6 vs 3+3 | same |
| m4b1 (rhse sp wedge 48x4x4) | cel, fofc, nghost 3 (fref none / clean) | 1 | 10 vs 5+5 | same / same |
| reg.sh RT hcyl sp3mb feed thin1 feedx1 sp3mbmr feedmr | defaults | - | 16 vs 8+8 | see RSTFIX.md |

- So the R136 defect is likely tied to something these decks lack. Candidates:
  - GPU-only ordering (a race, or the async vet_gd build);
  - the R136 pgen state: x1 user BCs, wb_arad / force_reference_work = split;
  - `time2_auto_be` skip runs that start before the dump;
  - the gd-SC stale vgd_i that the Orion gd-SC worker reported on sp (sp3 2.7e-6, flz q 0.5; being fixed on its own
    branch).
- **Method that found #14 in under an hour.** Dump every RadiationM1 / cln_* / vgd_* view at CleanSolve entry on the
  restart cycle, in straight and restarted runs. Then zero groups of the arrays that differ in the straight run until
  it reproduces the restart run. The scratch include (never committed) is
  /orion/ptmp/jinma/clean_1115/rstfix/g/dump.inc, built from rad_m1.hpp by the generator in RSTFIX.md.
- Please send the R136 deck if you want Orion to run it on CPU.

## Minimal patch (staggered #14 fix; for viper's clean-cc if wanted)

```diff
--- a/src/rad_m1/rad_m1.hpp
-  void CleanRightHandSide();
+  void CleanRightHandSide(const bool beta = false);   // (beta: CLW_BT, CLW_RE only)
--- a/src/rad_m1/rad_m1_clean.cpp
     if (cln_feedback) CleanAdvectiveCoefficients();
+    // (the first iterate's beta and emission for the drag's face flux, not the last
+    // solve's: a restart, failed stage or discarded attempt does not carry those)
+    if (cln_vimp && cln_feedback) CleanRightHandSide(true);
     if (cln_vimp && cln_feedback) CleanStagDrag(true);
@@
-void RadiationM1::CleanRightHandSide() {
+void RadiationM1::CleanRightHandSide(const bool beta) {
@@
     if (x.sbi) CleanGhostFill(cln_w, CLW_GX, 1);
   }
+  if (beta) return;
   CleanArgs xh = x;   // (the face momentum state, clean_stag_split)
```
