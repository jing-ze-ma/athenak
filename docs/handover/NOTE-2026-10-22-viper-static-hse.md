# NOTE 2026-10-22 viper -> caltech: static seed-free HSE on cs and sp: symmetric to round-off; no cs well-balance fix needed (agrees with 9c432f17)

Metric (scr/static.py, /viper/ptmp2/jinma/photinst_1022): per shell m=0 mean <v_r>, lateral rms of v_r - <v_r>, corr with
g = dOmega/<dOmega> - 1 (cs) or sin(theta) (sp).
(A) he_star_m1 RSG E, cs N=32, he_seed 0, no MLT scaffold (mlt_closure frozen), radiation coupled (he_star_m1 refuses inert radiation:
    force_reference = wb_arad needs coupling + gas_feedback), f55be371, to 4 d: lateral rms of v_r EXACTLY 0 on every shell 0.3-0.95 R_ph
    (one isolated 1e-6 cm/s); only an m=0 radial drift (the envelope readjusts without the scaffold): <v_r> +6e1 (0.3), -8e1 (0.5),
    +1e2 (0.7), +4e2 (0.8), -4.8e3 (0.9), +1.0e4 cm/s (0.95 R_ph) at 4 d.
(B) the same on a spherical-polar wedge (theta 60-120 deg reflect, phi 0-60 deg periodic, 24x24 = 2.5 deg cells), to 15 d: also exactly
    symmetric (lateral rms 0; 2e-3 cm/s at 0.95 R_ph at 15 d), and the SAME m=0 drift as cs to 3 digits (0.95 R_ph at 3 d: +5.72e3 sp vs
    +5.75e3 cs; 0.9 R_ph -2.85e3 vs -2.86e3).
(C) pure-hydro cs_test iprob 13 (isothermal HSE shell, constant g, PLM+HLLC), t 0-20 (~1750 cycles), CPU build: N=32 and N=64 lateral
    v_r rms 1e-16..5e-16 (round-off), m=0 drift 3e-6..2e-5 decaying; exact_panel_ghosts = 1 bitwise the same as default;
    <hydro>/cs_wellbalanced_src = true RAISES the lateral noise to 1e-12..1.5e-11 with |corr(v_r, dOmega)| 0.4-0.7 (the face-sum source
    leaves a small cell-geometry residual) -> not a fix; we do not propose it.
    Side note: the MI300A GPU build of cs_test segfaults at setup for iprob 13 (hydro and zero-field MHD) — GPU-only bug in that test pgen,
    not pursued.
=> A spherically symmetric star stays symmetric on cs to round-off, pure hydro and full RSG E physics; cs and sp behave identically.
The panel imprint needs lateral perturbations (seed, convection, cliff runaway) and comes from how their growth depends on cell size/shape.
Still running on viper: RSG E harmonic restart arm (40 -> 60 d), base / kirchhoff continuations to 80 d (scored by corr with dOmega).
