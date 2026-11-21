# NOTE 2026-11-21kp orion -> caltech: USER GO — restart all 4 production runs with the fwb_stage fix (t9r back to cfl 0.9)

Cause of t9r's stage divergence (reproduced on Orion A100, t9h cfl 0.9 + harmonic + oh3): the production fwb evaluates its residual s
from the current iterate, which makes the stage-2 pass map diverge at large dt (the NaN is the diverged solution; fwb/fcone only print
it). Fix: rad_m1/clean_stag_fwb_stage=true (s = per-stage constant of E_old) - same A100 deck 400 cycles 0 failures, 2.15 solves/stage.
1. Build B_fwbs from branch caltech-test-fwbfast = f496d606b (= caltech-test-fdef 03b027502 lineage + 1 commit), CUDA 12.9 + g++13, smoke gate.
2. For EACH of t9l62, t9r, m35mms, i1_m35mms_f1k: stop after the latest complete periodic rst (back it up as before), restart with B_fwbs:
     keys from the rst + rad_m1/clean_conv_max=1.0e-4 + rad_m1/clean_face_opac=harmonic + rad_m1/clean_stag_opac_update=true
     + rad_m1/clean_stag_flim=false + rad_m1/clean_stag_fwb=true + rad_m1/clean_stag_fcone=true
     + rad_m1/clean_stag_fwb_stage=true                       (NEW; a restart file lacks it, the command line sets it, NOTE expected)
     + rad_m1/clean_stag_thinf_held=false                     (USER 2026-11-21: thinf_held dropped at the next restart)
     + time/cfl_number=0.9                                    (t9r back to 0.9)
   NO clean_pass_damp, no gate, nothing else.
3. Fallback per run on CW / FATAL / dt_retry in the first 100 cycles: cfl 0.6 same keys; if it persists, flim at cfl 0.6 (as kl).
4. Post per run the command line + the first 50 cycles (dt, passes, CW, FATAL, KE, Ltop/Lin, CLEAN FWB/FCONE lines).
Holders end 17:46 PST (successors are a separate USER OTP decision; they should carry these keys).

-- orion (for the USER)
