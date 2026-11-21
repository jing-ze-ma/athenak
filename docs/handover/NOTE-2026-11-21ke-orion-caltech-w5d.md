# NOTE 2026-11-21ke orion -> caltech: W5d = rotation ON with the cost/dt fixes (one arm, next free test node)

W5 showed dt pinned at 55 (block alignment on a restart without vet_gd_rotate_stretch) and 1.42x cost (clean_rot_fuse is a
fresh-run default only). One arm, same setup as W5b (B_v207, copy of the same t9l62 rst, + conv_max 1e-4), ON keys:
  rad_m1/vet_gd_rotate_dt=5500 rad_m1/vet_gd_rotate_t0=<t_rst - 2750> rad_m1/vet_gd_rotate_stretch=0.02 rad_m1/clean_rot_fuse=true
  rad_m1/vet_rotate_avg=4 rad_m1/vet_rotate_avg_ema=128 rad_m1/vet_rotate_avg_ema_wb=true rad_m1/vet_rotate_avg_ema_change_form=step
  rad_m1/vet_gd_rotate_xfade=true rad_m1/vet_rotate_advance=true rad_m1/clean_rot_ring=live
150 cycles, rst off. Report as W5 (CW, s/cycle vs W5a 1.95, dt range vs W5a 59-67, GPU mem, KE/L/Lz series) -> V211/W5d.
Priority: after W4e/W4f.

-- orion (for the USER)
