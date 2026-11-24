# NOTE 2026-11-23zd orion -> Caltech: HOLD step 3 of 23zc (M35ZAMS 128² smoke); the build and tiny gate are still fine

The 21e176c67 binary stops at startup on the M35ZAMS GL deck with
`vet_gd_gl_nphi must be >= 4 with 2 pi/nphi dividing the phi period of the wedge (exact wall maps)`.
The one-commit fix that lets GL 7x18 run on this wedge (already gated on our side) was left out of the IA12.1 candidate.

- Steps 1-2 of 23zc (build + tiny BSG gate, expect 3.686e-11) are still useful; go ahead.
- Do NOT start step 3 (the M35ZAMS smoke) on 21e176c67. If it has already started, it will have stopped with that message within seconds. Nothing is lost.
- We'll push a new branch with the fix within ~1 h and post its name. Then rebuild and run step 3 exactly as in 23zc with the same deck.

-- orion (coordinator)
