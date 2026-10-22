# NOTE 2026-10-22 viper -> caltech (USER): <rad_m1>/time2_posfix wired for <hydro>/am_vector (re TASK 92c087ad, E2 config f0bd7ea1)

Branch fork/posfix-amvec-1022 = ed84cfed8 (on fin51-r-fix 1976bac68; NOT fast-forwarded):
 bd3974985 Time2FormStage(1) saves am_spin (t2ssv) with t2gsv; Time2PosFix puts it back and drops aminc in the marked cells;
           Time2FormStage(2) adds ccomp*amk1 back there (same stage-2 compensation as t2inc/csinc). posfix + am_vector no longer
           FATAL; MHD still FATAL. he_star_m1: problem/am_hist = true (named only) switches the user hst to the L budget
           RadM1AmHistory (Eg Er Lg Lr S Lf).
 ed84cfed8 (extra fix, found on the way) Time2Restore (stage rejection -> backward-Euler redo) also puts am_spin back to U^n
           (am_spin1); was an L leak per rejection. am_vector only; NOT exercised by the gates (0 BE steps in the stretch).
Gates, viper MI300A apudev, RSG E N=32 restart (rsg.00024, 78.4 d, cycle 7439):
 bitwise vs 1976bac68 (rst + hst): default IDENT; posfix=4 without am_vector IDENT (4 firings); am_vector without posfix IDENT.
 am_vector + posfix=4, 100 cycles: C2==C IDENT, R==C IDENT (restart at 7455, before the firings at 7477-7535); 18 re-solves on 17 cycles.
 L budget (gas + rad + spin, x/y/z): fix and no-fix arms give identical gas and radiation; without the fix the spin gains a step at
 each re-solve, total +1.57e44 / -8.0e43 / -6.8e43 after 18 re-solves (6 / 23 / 11 % of the spin reservoir); with the fix the
 per-cycle spin change on firing cycles matches the others. The stretch's total L change (~2e46, boundary-dominated, same in both
 arms) is not a closure test.
E2 rerun with am_vector (for Caltech to run): binary from ed84cfed8, restart prodF rsg.00025 (125 d), E2 keys plus am_vector:
   rad_m1/implicit_face_opac=harmonic rad_m1/time2_vet_col=extrap2 rad_m1/time2_accept_fallback=true
   rad_m1/time2_posfix=4 hydro/am_vector=true time/tlim=1.1225e7     (optional: problem/am_hist=true)
 am_vector needs rad_m1/implicit_am_cs = true (the RSG E decks have it). rsg.00025 has no spin: the run warns and starts the spin
 from 0. Compare against E2: dt record through 129.8 d, number of posfix re-solves, 0 FATAL, L budget (am_hist) across re-solves.
