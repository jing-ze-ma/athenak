# NOTE 2026-11-20cy orion -> caltech (USER): trim the beam arms; no unnecessary tests

- **Drop BFXC20x2_hg/_rsg and BFXC19x2_hg/_rsg.** The USER rejects the xcorner=2 partial, because beams must cross without artifacts.
- **BFX19:** when a star's BFX20 arm ends without FATAL, stop or skip that star's BFX19 arm. Run BFX19 for a star only if its BFX20 fails.
  - BFX19_hg is already running, so stop it if BFX20_hg passes.
  - Start BFX19_rsg only if BFX20_rsg fails.
- **From now on, beam arms are 100 cycles, 1 node each (NR=4), all concurrent.**
- Keep G_on_hg / G_off_hg (600 cycles, p300 partials), then the v20.2 gate (NOTE cx), then the paused arms.

-- orion (for the USER)
