---
name: cs-seam-multiblock-halo-fix
description: cs cell-centred seam halo was WRONG with >1 MeshBlock per panel (2x2 11x, 4x4 100x): clamp to the sending BLOCK's range on an extrapolating stencil + x2x3 edge ghosts never resampled; fixed 94c7165d on cs-seam-4x4; bvals_fc.cpp (MHD) still has BOTH defects
metadata:
  type: project
---

09-17 (wt_seam4/tests_seam4/README.md). Symptom: diffusion-test seam residual 3x worse at 4x4
blocks/panel, for the explicit operator too. Static halo scan (new `<problem>/seam_halo_scan`
in cs_test iprob 15, nlim 0) localised it to block ends facing the panel centre (the resample
pulls the donor position toward the centre by (layer+1/2) dx sin 2a; zero at a = 0, pi/4).
Defects in src/bvals/bvals_cc.cpp PackAndSendCC do_cs: (1) `b` clamped to the sending block's
[kl, ku-2] and the value clamped unconditionally -> at an interior block end the quadratic
extrapolant is replaced by the nearest node (O(dx)); fix: skip the value clamp when |u| > 1
(bvals_fc already had that guard). (2) x2x3 edge buffers (slots 40-47) across a seam were a
plain copy; fix: decode the flanking seam face and resample with cs_seam = 2/3; clamp to the
source's active range. After: 2x2/3x3/4x4 = 1x1 to printed digits; 1x1 BITWISE; 4 cs tests
pass; blast dt history unchanged. Residual 1.24x at n=32 4x4 = the extrapolation (O(dx^3)).

**OPEN:** bvals_fc.cpp (face-centred B) has BOTH defects -> every cs MHD run with >1 block per
panel (dhj cs MHD productions!) carries them; needs MHD gates. Corner slots 48-55 still plain.
See [[cs-seam-colocation-fixed]], [[cs-narrow-block-resample-degeneracy]], [[he4-presn-global-plan]].

**RETRACTED 09-21 night (tests_seam_fc/README.md, commit after 0eb7e2c4):** the "OPEN: bvals_fc has BOTH
defects" statement above was a code reading, not a measurement, and is WRONG. Static scan of an analytic
div-free B (cs_test iprob 11, `seam_halo_scan_fc`): bvals_fc guards the clamp (fpos in [0,2]) and resamples
edge AND corner slots 40-55; 1x1/2x2/4x4 identical at n=64,128; controls with the guards removed are
164-858x worse; dynamic MHD gates identical, div B 2.5e-14. The cs MHD productions (2x2) never carried a
B seam-halo defect. STILL OPEN: bvals_cc.cpp CORNER slots 48-55 are a plain copy (fc control suggests
it matters: 400-900x in the corner bin) -> fix after the cs-restart agent finishes (same tree).
09-22: bvals_cc CORNER slots 48-55 FIXED (commit after 2db7095c, tests_seam_cc_corner/): corner-seam bin
1.4e-4 -> 8e-8 (2x2), 7.4e-4 -> 5e-6 (4x4) at n=32 = the edge-seam bin; vertex fill untouched; blasts
bitwise in 8 layouts. No current dynamic gate reads these ghosts (matters for multi-D stencils at a radial
block boundary on a seam + prolongation/AMR). Remaining cs halo item: the 3rd-order interior-block-end
extrapolation residual (needs the receiver-side restructure).
