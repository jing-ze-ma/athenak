# NOTE 2026-11-21kl orion -> caltech: USER — the production switch (ki/kj) ALSO takes F_B (fwb + fcone, flim off)

USER (2026-11-21 ~14:20 PST): "let's do F_B" in the production switch. This REPLACES the binary/keys of ki (procedure of kj unchanged):
  binary   B_fdef (caltech-test-fdef 03b027502; build per kk if not yet built, smoke gate)
  keys     keys from each run's restart file
           + rad_m1/clean_conv_max=1.0e-4
           + rad_m1/clean_face_opac=harmonic
           + rad_m1/clean_stag_opac_update=true          (pass cap 3 = default)
           + rad_m1/clean_stag_flim=false rad_m1/clean_stag_fwb=true rad_m1/clean_stag_fcone=true
           + time/cfl_number=0.9
           NO clean_pass_damp; do NOT add thinf_held_gate; nothing else.
Fallback per run if CW / FATAL / dt_retry in the first 100 cycles: first drop to cfl 0.6 (same keys); if it persists, go back to
flim (rad_m1/clean_stag_flim=true, fwb/fcone false) at cfl 0.6. Report which run is on which setting.
The V213 F_A/F_B pair (kk) is then no longer needed for t9l62 — skip it unless nodes are idle; W4g/W4h still run.
Post per run: the restart command line + the first 50 cycles (dt, passes/stage, CW, FATAL, KE, Ltop/Lin, the CLEAN realisability line).

-- orion (for the USER)
