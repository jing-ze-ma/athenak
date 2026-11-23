# NOTE 2026-11-23m viper -> orion: CORRECTION to 23l item 2 — the constraint to remove is the block-size / reach one

USER (11-23 ~10:45 PST): "I mean cannot use 128 on 4 gpus".
- The constraint the USER wants gone is the vet_gd near-tangent LATERAL REACH limit, not (first) the GL phi-period rule:
  the paper's M35ZAMS wedge, 384 x 128 x 128, cannot run exactly on 4 GPUs. Splitting it into 4 blocks of 384 x 64 x 64
  gives blocks narrower than the reach (127 lateral cells here: wl = ceil(2 sqrt(r^2 - r_{l-1}^2)/(r_{l-1} dtheta)) + 2 with
  dr/r 4.9e-4, dtheta 5.0e-4), so lateral reads are clamped (rad_m1_vetgd.cpp ~7729-7742 on bd9553537: FATAL unless
  vet_gd_allow_clamp, and then the sweep is not exact). The only exact layout is one 128 x 128 block on one GPU.
- Please make the sweep exact for MeshBlocks smaller than the reach (a band deeper than one neighbour block, e.g. the raw band
  over several blocks / ranks, or the global band of vet_gd_band_global made to work on a wedge), so 128 x 128 runs exactly
  on 4 GPUs (2 x 2 blocks of 64 x 64), with a gate: 1-block vs 4-block run of the same deck identical to round-off (or
  layout_invariant bitwise), plus the speed on 4 GPUs.
- The GL phi-period rule of 23l item 2 is secondary: only if cheap.
- Item 1 of 23l (smoke IA11_m35zams 256 x 256 on Raven) stands unchanged.

-- viper
