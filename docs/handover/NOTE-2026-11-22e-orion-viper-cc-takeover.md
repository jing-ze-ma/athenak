# NOTE 2026-11-22e orion -> viper: USER — Orion takes over clean-cc development

Thanks for NOTE 22d (+ correction e434b06db). USER (11-22): Orion now owns clean-cc (cell-centred clean) development, starting
from fork/clean-cc-int-viper 06c9ee3cc. Goal: one cc scheme for all production (RSG/HG/R136/T9/T19/BSG) and the test suite,
conserving total energy and momentum (L_z too on sp); on the cubed sphere only radial momentum, L_z and total energy are required.
Plan: port the layout-independent stag2 fixes (positivity by construction, WBFCT write-back limiter, floor_pay, xpat=pos,
conv_max 1e-5, gas Newton switch, harmonic face opacity, ray ring) into cc, make F eliminate locally if the joint solve is the
bottleneck, fix cs radial momentum, then full gates (incl. speed vs prod).
No action needed on viper; please don't push to the clean-cc branches meanwhile. We will ask if we need run dirs or history.

-- orion (coordinator)
