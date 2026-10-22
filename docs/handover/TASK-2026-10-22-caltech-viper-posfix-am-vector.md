# TASK 2026-10-22 caltech -> viper (USER): wire <rad_m1>/time2_posfix for <hydro>/am_vector (re f0bd7ea1)

Today time2_posfix > 0 with am_vector (or MHD) is FATAL (rad_m1_time2.cpp ~l.461: "not wired for MHD or <hydro>/am_vector"):
Time2PosFix puts the gas back to the stage-1 start (t2gsv) but not the am_vector spin reservoir's increment.
USER asks: wire it (save/restore the spin reservoir alongside t2gsv in the re-solve; the dropped part re-entering stage 2 must carry the
reservoir consistently), default unchanged (bitwise when posfix is off or never fires). Gates: bitwise key-off; C2==C and R==C with
am_vector + posfix on; L (all 3 components) budget over a posfix-active stretch. Then a validation run: repeat Caltech arm E2 with
<hydro>/am_vector = true (RSG E rsg.00025 at 125 d -> 129.92 d; E2 config in f0bd7ea1) -- Caltech can run it if you send the commit.
Meanwhile Caltech relaunches RSG E + HG with am_transport (L_z) + posfix (option 1, USER); am_vector only after this passes.
