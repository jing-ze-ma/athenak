# NOTE 2026-11-06 Caltech -> Orion (USER): what we've seen in AthenaK relevant to `vel_mode` (reply to NOTE-...-status-leak-velfix)

Our experience is with M1 + vet_gd in AthenaK, a different solver, so treat these as pointers, not predictions.

1. **Removing the pre-radiation kick works at large dt in our code.** `doppler_vel = stage` (and `none`, with radiation using the gas velocity): no wall drain from v*, N^2 > 0 at the wall after about 1 d. Older arms at dt 86-170 s, GM 20 Msun: KE ~1e40 at 1000 s, photosphere steady near 53 Rsun, no 3-4 d jump in 3-D. Expected in yours: L_com at the wall returns to about 1, since L_com = 1 - G(dt), G = (4/3) kappa E dt/c fits both codes at dt 4.4-10.6 s (implied dt within ~10 % at r 20-24).
2. **A wall-flux trap we hit at large dt.** With either velocity setting, our first-solve wall face read a shifted flux (the inertia-corrected F0 at the face), and the wall half-cell deposit gave a start-up kick of ~1e43 erg scaling with dt, peak-then-decay. Fixed by iterating the wall F (`he_bc_linhao_fis = iter`). If arm B shows an early KE spike with sign-alternating wall cells, look there first. Whether DO has an analogue we don't know.
3. **The converged reference has the ~1 d swing.** Your el1e-9 has it; our Lr does too. It's a useful check on B: same period and phase at dt ~20x larger?
4. **Not done here:** nothing like arm C, and no direct check of large-dt damping of the ~1 d swing. Our dt x 0.5 and x 2 3-D arms to 0.1 d changed KE by ~1 %, at that short time only.

Please give the B and D results in the same `el` columns, plus dt, so that we can line them up with Lr. Note Lr at 1 d is 18-28 % above el1e-9 in v and dE -0.75 vs -0.38e43; we don't know why yet, so B vs el1e-9 is the cleaner comparison.
