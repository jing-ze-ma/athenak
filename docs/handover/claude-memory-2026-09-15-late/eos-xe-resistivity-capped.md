---
name: eos-xe-resistivity-capped
description: "RESOLVED 2026-08-16 (commit c0edbe49): the x_e was never wrong -- two out-of-bounds reads of eta_b in AddEMFGeneralResist wrecked the run, and one of them made ohmic_resistivity=constant a no-op"
metadata: 
  node_type: memory
  type: project
  originSessionId: 007bb828-7120-45a1-9fd1-74b062edff22
  modified: 2026-08-16T07:25:56.190Z
---

Found 2026-08-15, diagnosed and FIXED 2026-08-16 in commit `c0edbe49`. Both suspects in
the original note were wrong: the tabulated x_e and the temperature handed to it were both
correct (in-kernel print gave xe = 3.9975e-8 against a model value of 4.39e-8 at a grid
centre 15 K hotter, which is exactly right for d ln x_e/d ln T ~ 13).

## What it actually was

Two out-of-bounds reads of `eta_b` in `AddEMFGeneralResist` (`src/diffusion/resistivity.cpp`),
the only live EMF routine — `AddEMFConstantResist` is commented out at its call site.
Reading past a Kokkos View is silent in a Release build.

1. **`eta_b` was allocated only in the non-constant branch of the ctor.** With
   `ohmic_resistivity = constant` it stayed at its `1x1x1x1` placeholder, every read
   returned zero, and **`eta_ohm_const` did nothing at all** — runs at eta = 0, 2.56e11 and
   1e13 were bit-for-bit identical. `NewTimeStepGeneralResist` reads the same array, so the
   diffusive timestep was inert too. Fixed by always allocating and `deep_copy`ing
   `eta_ohm_const` into it.
2. **In 1D and 2D the edge fields carry a `k = ke+1` layer (1D also `j = je+1`) but the
   CELL-centred `eta_b` does not**, and the code indexed `eta_b(m,ke+1,...)` — one whole x1
   row past the end. 3D is in bounds, which is why it survived. Fixed by using the single
   cell's eta on both faces.

Bug 2 is what made `ohmic_resistivity = eos` look broken: it destroyed a 5-cell blob within
the first cycle, which floated everything to `max_eta` and collapsed dt. Bug 1 is why the
`constant` cross-check looked "stable" — it was doing nothing.

**Diagnostic that cracked it:** printing `eta_b` for the WHOLE row instead of one cell. One
cell showed a perfectly correct eta while the state was already wrecked elsewhere, which
ruled out the EOS immediately and pointed at the update.

## Verified after the fix

Alfven wave, rho = 1e-6, T = 1985 K, eta = 2.567e11: the measured decay matches
`exp(-eta k^2 t)` (rate `eta k^2`, NOT `2 eta k^2` — only the magnetic half of the wave
energy is dissipated) to 0.35%, `eos` matches `constant` at the same eta to 0.1%, and
[M/H] 0 -> 0.5 changes the damping by 1.8409 against a banner x_e ratio of 1.8378. 2D
reproduces 1D to five digits.

Related: [[eos-electron-regression-test]] (now written and passing),
[[resistivity-perna-uhj]], [[general-eos-stage3-table]].
