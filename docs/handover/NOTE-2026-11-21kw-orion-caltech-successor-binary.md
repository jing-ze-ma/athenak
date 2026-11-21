# NOTE 2026-11-21kw orion -> caltech: successor holders (after 17:46 PST) — new production binary B_p1 = bf4f917d4

Branch caltech-prod-bf4f917d4 (gated: radshock exact, stars hst identical with old keys, CPU + A100 restarts bitwise).
= the kp binary f496d606b's fwb_stage + gas Newton switch + harmonic/opac fixes + iterate-restart crash fix; xpat default OFF.
Build B_p1 with CUDA 12.9 + g++13 (never 12.2), run the tiny BSG smoke gate (pass 1 4.020e-06), record md5.

On the successor holders restart all four runs from their newest rst with B_p1 and the kp keys PLUS (named explicitly,
restarts otherwise keep the recorded values):
  rad_m1/clean_gas_newton_switch=4  rad_m1/clean_gas_newton_switch_min=6
  (keep: clean_stag_fwb_stage=true, clean_stag_thinf_held=false, fwb+fcone, flim off, clean_conv_max=1e-4; NO clean_pass_damp,
   NO time2_k1_xpat)
  env: export VGD_XFE=0   (sp early send off: removes a shared-GPU race; identical bits)
cfl: t9l62, m35mms, i1_m35mms_f1k at 0.9; t9r at 0.6 with hydro/tfloor = 1000 K as in NOTE iz (USER).
The switch cures the stage-stall class behind today's t9l62/m35mms fallbacks (V214 D_c). Post first-50-cycle tables as before.
V214 D_x (xpat) stays a diagnostic only.

-- orion (for the USER)
