# NOTE 2026-11-21ka orion -> caltech: OPAC W4 = RSG + HG off vs oh (USER), BEFORE W3

USER wants the opacity keys checked on RSG and HG too. Replace W3 (t9r / m35 off/oh: t9l62 W1 already answered the T9 restart question)
by W4, to start as soon as the W2 T19 arms end:
  W4a RSG off  / W4b RSG oh : the RSG restart R2 used (rsg_E_1020/prodH/rst/rsg.01295.rst; R2 ran it clean at damp 0), NR=4, 150 cycles,
                               R2's keys (NOTE jt list) on B_opac (00ad333b7); oh = rad_m1/clean_stag_opac_update=true + rad_m1/clean_face_opac=harmonic.
  W4c HG off   / W4d HG oh  : the HG restart your C2F_hg arms used (say which), NR=4, 150 cycles, its production keys + conv_max 1e-4 + thinf_held_gate.
No clean_pass_damp. rst off. Report as V210 (series every 10 cycles of KE, L_top, L_in, L_bot; CW, FATAL, passes/stage mean/max, dt, s/cycle).
Then W3 only if nodes are idle.

-- orion (for the USER)
