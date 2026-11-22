# NOTE 2026-11-21jk Caltech -> orion: lk DONE — t9l62 on head 462210d51 from t9l62.00045 (USER confirmed in our session)

- B_j3 run stopped 22:58 (it had recovered past c33461, CW 5 / dt_retry 3 / FATAL 0, at c33600). Restart file: t9l62.00045 (cycle ~33350, written 22:41, before the first CW at c33458); later B_j3 rsts moved to t9l62/rst_j3_late/.
- Head 462210d51 built in place on testH1 (CUDA 12.9 + g++ 13, 20 files vs bf4f917d4, ~5 min), md5 2144006b9576934a061f6b2b9808a67e. Gate: pass 1 4.934e-06, NOT converged after 50, FATAL 0.
- Command-line keys: time/cfl_number=0.9 rad_m1/clean_conv_max=1e-5 rad_m1/clean_stag_flim=false rad_m1/clean_stag_fwb=true rad_m1/clean_stag_fcone=true **rad_m1/clean_stag_fwb_stage=true** output5/dt=21600; env VGD_XFE=0.
  * flim had to be named false: the J4 deck/rst records clean_stag_flim = true, so the head default would not apply on restart.
  * fwb_stage had to be named true: the head echoed "clean_stag_fwb_stage = false (a restart file without the key: the old default)" on the first launch (killed after ~1 min, relaunched 23:06). fwb_stage is the USER-approved kp fix for the cfl-0.9 fwb pass-map divergence (t9r 26418 NaN). Tell us if you want it off.
  * thinf_held stays as recorded in the rst (true).
- Running on testH1 until END 00:20 (graceful stop available on the head). First cycles, CW/dt_retry and the jh noise metric follow.

-- Caltech
