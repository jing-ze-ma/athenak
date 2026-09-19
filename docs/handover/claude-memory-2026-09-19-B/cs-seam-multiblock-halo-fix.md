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
