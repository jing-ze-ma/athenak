# NOTE 2026-11-21kb orion -> caltech: W5 = rotating rays on t9l62 (ROTWB), on the first node free after W4

USER wants rotating rays for wedge-star production. One test node, v20.7 binary (B_v207 86d31021), from a COPY of the current
t9l62 production rst (as V8a), keys from the rst + rad_m1/clean_conv_max=1.0e-4 + rad_m1/clean_stag_thinf_held_gate=true. NO clean_pass_damp.
  W5a OFF : as is, 150 cycles.
  W5b ON  : + rad_m1/vet_gd_rotate_dt=1100 rad_m1/vet_rotate_avg=4 rad_m1/vet_rotate_avg_ema=128 rad_m1/vet_rotate_avg_ema_wb=true
            rad_m1/vet_rotate_avg_ema_change_form=step rad_m1/vet_gd_rotate_xfade=true rad_m1/vet_rotate_advance=true
            rad_m1/clean_rot_ring=live , 150 cycles. (rotate_dt as a NUMBER = 20 x prod dt ~55 s; not auto.)
  W5c ON rt: 10 + rst + 10 with the W5b keys, must be bitwise (hydro_w / m1 bins + hst).
Report (caltech_results/V211/): CW, FATAL, s/cycle, GPU memory peak, dt range, series every 10 cycles of KE_int, L_top, L_mid, Lz_tot;
W5c bitwise yes/no. Order: after W4 (RSG/HG off vs oh) on the first free test node.

-- orion (for the USER)
