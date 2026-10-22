# NOTE 2026-10-23 viper -> caltech (USER): MLT ramp rule v2 merged (f = 0.61); fixed-GM in n2seed_profile_1022.py; fin51-r-fix = d01ce628c

fin51-r-fix fast-forwarded 0429ff46b -> d01ce628c (= + n2seed-1022 a766d4f8: ramp rule v2 of TASK ca0602bb).
Rule (he_star_m1, fresh start, ramp keys auto = default): mlt_ramp_start = 0, mlt_ramp_time = t_sat = mlt_ramp_nefold (10.5) /
(mlt_ramp_sigma_frac (0.61) x sigma_max), sigma_max from the code's IC cells (code EOS + aT^4/3, Gamma1 from the first law, g = G m(r)/r^2).
Stored in the rst trailer, reused on restart; explicit numbers override; -1 = no ramp; an old rst without stored numbers gets no ramp.
Printed: RSG E 1/sigma_max 6.11 d -> t_sat 9.09e6 s = 105.2 d (matches your ~105 d); HG global 1/sigma_max 0.775 d -> 1.15e6 s = 13.3 d
(your deck: 7.3e5 s = 8.4 d -- see below).

WHY f = 0.61, not 0.74 -- your script docs/handover/n2seed_profile_1022.py uses g = GM/r^2 with GM FIXED at the column bottom:
 - RSG E: the envelope is massive; at the N^2 peak (0.995 R_ph) G m(r) = 1.84 x GM, so the script's sigma is low by sqrt(1.84) = 1.36.
   1/sigma_max: script (file points, fixed GM) 7.35 d; file points with G m(r) 5.42 d; code cells with G m(r) (what the run feels) 6.11 d.
   Your measured e-fold (~10 d at 0.995 R_ph) against the code's sigma_max gives f = 6.11/10 = 0.61 -> t_sat 105 d, same as you measured.
 - HG: G m(r)/GM = 1.013 there (gravity irrelevant); the difference is radial resolution: the N^2 peak (file dr/r 2.0e-4) is narrower
   than a code cell (dr/r 9.3e-4): file points 0.594 d, file N^2 averaged over code cells 0.686 d, code cells 0.775 d. The code-cell value
   is what the grid can grow. f = 0.61 is calibrated on RSG E only -- if you measure the HG e-fold, we can check it there.
 - Minor: radiation in Gamma1/P ~1 %, centred vs np.gradient < 0.1 %, balanced vs raw column < 1 %. With the code EOS your method
   reproduces your numbers exactly (EOS not the cause).
 Suggest fixing the script (g = G m(r)/r^2) for future calibrations.
Gates (prodmerge-1022e = d01ce628c vs 0429ff46b): old seed + ramp -1 and explicit ramp bitwise; old rst + auto -> no ramp, bitwise;
fresh n2 + auto C2==C, R==C (ramp read back); accretor D_B1v6 bitwise; RSG E N=32 stretch key-off and gpumem keys-on bitwise;
MI300A smoke clean. MI300A binaries: /viper/ptmp2/jinma/n2seed_1022/bin/athena_he_apu_d01ce628_r71 (2c068df6),
athena_apu_ry3d_d01ce628_r71 (c444bfe7). Your running RSG E (explicit ramp 43.0 -> 94.4 d) is unaffected.
