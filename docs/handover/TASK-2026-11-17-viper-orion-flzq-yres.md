# TASK 2026-11-17 viper -> Orion (USER): push 17 is on clean-cc; which clean_yres is right on flz_q?

## Push 17 released
- **fork clean-cc = 5785ee51d** = b66ff2552 + one commit putting clean_yres back to **false** by default.
- So the push-17 defaults are exactly your bitwise A arm. The USER decided this after your NOTE 1af09d842.
- clean_yres = local stays selectable. Please rebase stag2 onto 5785ee51d when convenient.

## The question
Under clean_yres = local, flz_q (sp fv1, chi = 1, rho 1e2, T_gas 1e-4, the cold radiation-dominated rotating box) changes at O(1):
- max|E - E16| is 0.111 of max 0.124;
- m1_f moves by 0.8;
- the 40-step L_rz changes 3 % (1.050e-3 vs 1.086e-3);
- Lz and E stay at round-off.

Which answer is right? Please run on CPU:
1. flz_q with clean_yres = false / post / local, at the harness resolution and at 2x in r and theta (and 4x if cheap). Also with clean_tol and
   clean_tol_nl 100x tighter, to separate solver error from scheme difference.
2. Compare each with the deck's expected answer:
   - the analytic or steady L_rz;
   - the sign of the wall v_r (inward, as your earlier "wall v_r inward like cell" check);
   - E and F profiles if a reference exists, e.g. a fully converged run.
   - Find which arms converge toward each other with resolution and tolerance.
3. A short physical reading. In a radiation-dominated, cold, optically thin-ish box, the gas momentum the cell's local exact exchange gives
   (local) differs from "minus the radiation change" (false). Is local's drag-only share closer to the true radiation force on the gas?
   Look at the force balance per cell if you can.
- Reply as NOTE-2026-11-17-orion-viper-flzq-yres.md.
- If local is the better answer, viper will make it the default in push 18. If not, we will find what local gets wrong.
