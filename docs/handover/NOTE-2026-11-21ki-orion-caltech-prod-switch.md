# NOTE 2026-11-21ki orion -> caltech: USER GO — switch the 4 running production runs to the new settings (cfl 0.9)

USER (2026-11-21 ~14:05 PST): switch t9l62, t9r, m35mms and i1_m35mms_f1k now, inside their current holders, to:
  binary   B_opac3 (caltech-test-opac 772ee5ac6, md5 35416663..., the one V212 ran)
  keys     keys from each run's restart file
           + rad_m1/clean_conv_max = 1.0e-4
           + rad_m1/clean_face_opac = harmonic
           + rad_m1/clean_stag_opac_update = true        (pass cap 3 = default)
           + time/cfl_number = 0.9
           NO clean_pass_damp; do NOT add thinf_held_gate; no other key changes.
Procedure per run (one at a time is fine; production nodes only, same node as now):
  1. graceful stop (STOP file in the run dir, per NOTE jk) -> confirm the final rst is written (name + cycle);
  2. keep a copy of that rst (backup);
  3. restart from it with the binary/keys above in the SAME run dir lineage (new log name, e.g. *_k9.log);
  4. post the command line + the first 50 cycles: dt, passes/stage, CW, FATAL, KE, Ltop/Lin.
If any run shows CW or FATAL in the first 100 cycles, or dt_retry, fall back to cfl 0.6 with the same keys for that run and say so.
Then keep the usual PROD table. Holders still end 17:46 PST (successors are a separate USER OTP decision).
Evidence: V212 t9l62 cfl 0.3/0.6/0.9 all 0 CW, cfl 0.9 2.57x faster per simulated second, dKE 8e-5, dLtop 3.8e-4 vs cfl 0.3.
Note: m35mms / i1_m35mms_f1k were not CFL-tested yet (V8 m35 restarts were clean at damp 0) - watch them closely.

-- orion (for the USER)
