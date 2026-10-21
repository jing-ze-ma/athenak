# pgrst_1021 NOTE (PGRST worker)
- pgrst-merge-1021 @ 077dda23: merge set (pgen restore, k1_reset opt-in, booking A, finalize frees). Bins bin/*_077dda23*.
- pgrst-1021 (wt): + hesdirk2 hold booking B + K1 booking (887e4b74/a8d7b72b) -- NOT correct, off-branch.
- pgrst-vimp-1021 (wtv): implicit_vimp_ewf (b3f989f8; no effect on the accretor stream failures),
  implicit_cs2_wbsolve (455a5265) PARKED (USER 10-21): closed-form conservative2 write-back kick under time2_vstage.
  Caveat: the 3-pass loop gain om(1+d11)E/(rho c chat) only reaches ~1 when E/(rho c^2) ~ 1 (reduced c, the 1-D radwave
  P33 test); prod E/(rho c^2) <= 4.9e-3 with chat = c. It WOULD matter with RSLA chat << c. Untested, no gates.
- pgrst-comb-1021 (wtc) = pgrst-vimp-1021 33f3a46e (+ time2_k1_reset=pos) + gsdrag dacap (0484c323, 822d1e61) + bbd6b430 (ewf in the vimp Jacobian).
  Stream SS (sA..sE arms, runs/gss, runs/cpuss, runs/cpuss2): every cycle still NOT ADMISSIBLE with k1r/grow/pos/dacap/ewf(+Jacobian);
  cause = implicit_vimp positivity fallback (S2 down to -33 in 1e-8 back-sink cells beside the pgen-held nozzle, ~2/cycle) -> stage rejected
  via vimp_now != impl_vimp; with vimp off + matched: stream-front gas eint old vector (stale K1_EN), pos catches some (2.3e4 cells) but not all.
  Accretor production = BE + matched (deck decks/v12lzc_be_gm.athinput c6af572f) on pgrst-merge-1021 077dda23.

## DESIGN implicit_hold (pgrst-posfix-1021, 10-21/22): pgen holds as Dirichlet rows of the stage solves
Problem: ry_per_accretor overwrites held cells (nozzle E/F/gas, back sink, donor lobe E, sink gas) in every hydro stage
(RyPerSrcEnv) = an operator split the hesdirk2 stage algebra cannot see: stale K1, vimp kicks into held light, cs2 dA of
pgen-held F, stream-front jumps. Symptom patches (k1_reset, grow, pos, posfix, vimp_kill, accept_fallback) cap at ~half
the steps second order.
Design (opt-in <rad_m1>/implicit_hold = true, default false -> bitwise):
1. Declaration: RadiationM1::hold_m(m,k,j,i) (0/1). The pgen sets it in the same kernel that writes the held radiation
   (one source of truth, every hydro stage, zeroed at the stage start); the held value is the u0 E the pgen just wrote.
2. Solve (every stage solve and BE): a held cell is a Dirichlet row exactly like an M1_IBC_EFIX end cell, at every place
   that special-cases EFIX: row replaced (E' = EN), its old vector EN = the held E (no t2inc / K1), its gas old vector
   without t2inc, no work rows, no vimp coefficients and no vimp kick (P, dv = 0 before the halo exchange: neighbours do
   not see a held cell's kick), dc instead of plm, no od cross terms.
3. Write-back: a held cell keeps E = EN, its F and its gas untouched (the pgen owns them); its stage slopes K (E, gas,
   cell-F) are 0, so stage 2 / FSAL carry nothing for it.
4. Booking (step 2, after the admissibility test): the Dirichlet hold supplies whatever its faces carry: book
   dt_w x sum over held cells of the outward face fluxes (F0 + advective A E) with the solve's stage weights (as booking A);
   the per-stage reset sums then vanish (held E = target already). EDBG check K arm closes.
5. Gas holds: the pgen's gas resets stay in the hydro stages (after the RK combination: already consistent); the solve
   no longer moves held gas, so the reset only undoes hydro transport, as in BE.
Test: ST+SS v12lzc hesdirk2 + cs2 + dacap(24810dad) + key: stage vs BE steps, fast cells, budgets; default bitwise gates.
