# NOTE 2026-10-22 viper -> caltech: kirch-1022 (61580891) built on MI300A; key-off bitwise; key-ON reproducible (re TASK 543cd795)

## Build
fork/kirch-1022 61580891 applied to the f55be371 (+fsub) build tree, incremental, rocm 7.1, he_star_m1, no extra flags ->
/viper/ptmp2/jinma/photinst_1022/bin/athena_he_apu_kirch61580891_r71 (md5 d90fca6d).

## Key-off gate (default implicit_face_opac = arith): PASS
RSG E N=32 base deck, fresh, 50 cycles (viper apudev 12319128): rsg.hydro.hst and rsg.user.hst bitwise identical to f55be371.

## Item 2: reproducibility with implicit_face_opac = kirchhoff: PASS
Restart from our growing RSG E state (base rst at t = 40.0 d, cycle 3456 -> 3506, rst every 25 cycles), 1 node x 2 MI300A:
- C vs C2 (rerun): both hst files and the rst at cycle 3506 bitwise identical.
- R (restart from C's rst at cycle 3481 -> 3506) vs C: hst rows identical; rst payload bitwise identical (the rst differs only in
  8 bytes of column padding in the embedded input text, output5 block; payload compared after the header).
- The key on a restart is accepted from the command line ("absent from the restart's input: added from the command line").
  Production-string smoke (-r rst00012 rad_m1/implicit_face_opac=kirchhoff) clean, no FATAL.

## Item 1: running
RSG E N=32 restart 40 -> 60 d (then 80 d) with implicit_face_opac = kirchhoff vs our base continuation (base: L(1.2-1.4) panel
contrast 1.1e-3 at 55 d -> 2.2e-3 at 57 d, first time2_reject retries at 56.1 d). Scored by panel contrast (primary), corr with
dOmega, cliff-shell sd, dt/retries, s/cycle. Results to follow.
Side note: our mlt_fsub_smooth = 4 restart arm (40 -> 49 d) cuts the cliff-cell e/rho sd 3-5x but RAISES the L(1.2-1.4) panel
contrast to +4.7e-3 at 49 d (base +6.6e-4 at 48 d; L sd 1.5 % vs 0.4 %); possibly a restart transient from changing F_sub, being followed.
