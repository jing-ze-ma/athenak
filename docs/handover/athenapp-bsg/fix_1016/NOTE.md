# NOTE athpp_fix_1016 (worker): converge Athena++'s implicit radiation solve (Linhao BSG, 1-D)
CHECKPOINT 1: tree athena/ = copy of athpp_audit_1015/athena_led (e47fabb7 + ledger + bsg_star pgen), branch conv-fix,
baseline commit 5b9978a0. Built with /usr/bin/g++ 7.5 -O3 (same as audit). Bitwise: L_Bd1 deck 0.2 d, hst/tab/ledger identical
to the audit run except the tlim-truncated last row. Diagnostics added (all read only when named): radiation/conv_log=1 -> out/conv.txt
(cycle stage time dt niter rel_err cum_sweeps cpu_s max_cell_rel), radiation/conv_dump=<cycle> (+conv_dump_all) per-sweep per-cell dump.
Runs: runs/<name> made by mkrun.sh from the audit decks (L_Bd1 = closed wall af1 tc4; L_tc4 = his BC af0 tc4; A_af1 = his BC af1 tc4).
CHECKPOINT 2 (existing Athena++ options, Bd1 deck, error_limit 1e-11 nlimit 5000, 0.1 d; sweeps per solve stage1/stage2):
- Jacobi default: ~770 / ~1450.   red_or_black=1: identical (~770 / ~1450).   iteration=0 (no explicit split): NaN at once.
- omega=1.2/1.5/1.8 (over-relaxation): NaN at once.   new symmetric Gauss-Seidel (sweep=1): 520-970 SGS iterations (= 2 passes each) -> no gain.
=> nothing Athena++ offers helps; the iteration is not streaming-limited.
CHECKPOINT 3 (diagnosis, per-sweep dumps runs/dJa, dS): the slow mode is a smooth isotropic drift of the whole thick interior (cells 0-180,
ratio 0.9917 per sweep, e-fold ~120 sweeps) driven from the photosphere (cells 200-230). Cause: ite_scheme 1 (default) puts the
'downwind' own-cell part of each face flux into exp_coef_ (lagged), i.e. (1 + X) I^k = RHS + X I^(k-1) with X ~ c dt/dx mu/2 ~ 250:
a heavily under-relaxed Jacobi; removing the split (iteration=0) diverges. Also: the analytic quartic root (FouthPolyRoot) and the
J/T elimination in the cell update carry cancellation noise ~eps*c dt chi (1e-9..1e-6 relative in the thick base); Jacobi only reaches a
deterministic fixed point of that noisy map.
CHECKPOINT 4 (the fix, commit 00c93bab, option radiation/sweep = 2; default 0 = untouched Athena++):
- rad_line_solve.cpp (new): per x1 line, assemble the full per-angle implicit system that CalSourceTerms solves cell by cell (neighbours,
  exp_coef_ own-cell split, implicit angular flux neighbour, J coupling, gas T^4 linearised about the current T) and solve it DIRECTLY by
  block-tridiagonal (block Thomas, nang x nang blocks, partial-pivot LU). RHS floor (clip) active set mirrored from the last cell update.
  Ghosts / x2,x3 neighbours lagged (exact in 1-D with one block).
- each iteration = line solve + one normal Athena++ sweep (restores the exact non-linear cell relation, floors, sets divflx_, residual).
- convergence: global L1 < error_limit AND per-cell max sum|dI|/sum|I| < error_limit_cell, or stagnation at the round-off floor
  (niter>=3, no 2x contraction, max_cell < error_limit_stag 1e-4).
- cons_bill (default on with sweep 2): divflx_/ang_flx_ rebuilt from the FINAL intensities before GetHydroSourceTerms (telescoping billing);
  stable now because the iterate is converged (the audit's ledger_fix re-billed an unconverged iterate).
- root_newton (default on with sweep 2): 2 Newton steps on the quartic root.
Deck: sweep=2 error_limit=1e-8 error_limit_cell=1e-7 nlimit=30.
CHECKPOINT 5 (2-d Bd1 test, runs/f2 vs default runs/j2 vs audit converged C_Bd1): outer iterations per solve 1-5 (mean ~2);
solve CPU 7.9 s vs default 56 s; whole run 23 s vs 69 s (fix is ~3x CHEAPER than default). Ledger at 1 d (1e45 erg): dE_tot
f2 -0.050 (non-telescoping 0.000, boundary -0.057) | C_Bd1 -0.008 (nontel 0.000) | default -1.071 (nontel -1.012).
30-d runs F_Bd1 / F_tc4 / F_af1 started (login node, 1 core each; Raven general/small queue blocked by QOSGrpCpuLimit).
CHECKPOINT 6 (bitwise, option off): binary athena/bin/athena_c2 (= commit 00c93bab + rhs_floor option) with the audit decks:
L_Bd1 (1 d) and L_tc4 (1 d): hst 102 rows, ledger rows, all tab dumps 0-9 IDENTICAL to the audit runs; L_tc4 and A_af1 0.5 d with
athena_c2: hst 53 rows + tabs IDENTICAL. => option off = Athena++ bit for bit.
CHECKPOINT 7 (30-d arms, login node 1 core each; Raven small/general blocked by QOSGrpCpuLimit):
- F_Bd1 (fix, closed wall af1): dE_tot -0.05 (1 d), +0.40 (2 d) vs converged ref C_Bd1 -0.008 / +0.43, default -1.07 / -7.84;
  non-telescoping 0.000. Profiles track C_Bd1 to ~3 digits at 0.5/1/1.5 d. BUT the base develops the cell-scale wall instability
  (cell 0 hot, cell 2 cold + dense: rho x2, T 1.3e6 -> 5e5 K, sigma_a ~5e10 at 2 d -- the SAME pattern in C_Bd1 at 1.5 d and in
  default L_Bd1 by 2 d) and the run dies (dt -> 1e-17, NaN) at 2.70 d. Variants: cons_bill=0 (V1) dies 3.4 d; no stagnation
  acceptance (V2) dies 2.7 d; rhs_floor=0 (N_Bd1, pre-solve clip off) dies ~3 d. => not caused by the billing or the criteria.
- F_tc4 (his BC af0 tc4): alive past 22 d (default died NaN at 5.07 d); dE_tot -31.1 at 5 d, ALL of it boundary (bath drain through
  his inner BC), non-telescoping 0 (default -37.7 = -25.8 boundary -12.0 non-telescoping).
- F_af1 (his BC af1 tc4): dE -9.4 at 5 d with resets +10.9 (eint<0 at the base); dt collapses to 1e-16 at 8.97 d (default A_af1 runs 30 d
  while losing energy; L_af1 default ledger: -32.5 at 5 d, -13.1 non-telescoping).
- Cost (same login node, concurrent, 2 d Bd1): default j2 69.4 s, fix f2 23.4 s (0.34x), converged Jacobi C_Bd1 318 s (4.6x).
  Solve CPU per star day: Bd1 default 28 s, fix 7.8 s; af1 default 32 s, fix 6.9 s; tc4 fix 17 s (more capped solves at the base).
CHECKPOINT 8: running converged-Jacobi restart of the fix's 2-d state (R_Bd1/outC, error_limit 1e-11) to test whether the converged
Jacobi solution also dies at ~2.7 d; exploratory F_tc1 (fix + taucell 1, his BC af0) to 30 d.
CHECKPOINT 9 (converged reference beyond 2 d): restart of the fix's 2-d state with Athena++'s own Jacobi at error_limit 1e-11
(R_Bd1/outC): dt falls 9e-3 -> 1e-4 at 2.68 d, sweeps hit 5000 with NaN residual, run dead at 2.684 d (fix died at 2.70 d).
At 1.9 d the fix and the audit's converged C_Bd1 agree to 3-4 digits in the base cells (cell 2: rho 1.711e5 / 1.708e5 code,
T 1.614 / 1.616 x1e5 K = the collapsed cold dense cell). => the Bd1 crash is a property of the CONVERGED Athena++ discretisation
(thick-limit near-central DO face flux + per-cell thermal/opacity runaway at the wall), not of the fix. Exploratory F_tc1 (fix +
taucell 1, his BC af0) also collapsed (dt 1e-15) at 15.3 d. Default Jacobi survives longer only because its non-telescoping billing
drains energy out of the base.
Figure fix_user.png (fixana.py; numbers fix_summary.txt, fix.json). Patch: conv-fix.patch (vs baseline 5b9978a0), commits 00c93bab,
9f2dbfe9 on branch conv-fix.
CHECKPOINT 10 (FINAL, 30-d results; fix_summary.txt):
- F_tc4 (fix, his BC af0 tc4) reached 30 d (default died at 5.07 d): non-telescoping 0.000 throughout; dE_tot -31.1 (5 d), -59.6 (10 d),
  -101 (20 d), -85 (30 d) = boundary (his inner bath) -31.2/-60.0/-150/-271 plus resets (eint<0/T-ceiling) +0.15/+0.47/+49/+193.
  T1 1.05e6 -> 3.9e5 (5 d) -> 1.7e5 (10 d) -> 0.6-1.2e5 K; r(rho=1e-10) 52 -> 40 Rsun at 10 d, envelope blows out to the outer edge
  (80 Rsun) at ~18 d. Solve CPU 25.8 s/day (mean 13.6 iterations; most capped at 30 after 14 d when the base is pathological); whole run
  894 s / 30 d.
- L_af1 (default, his BC af1) to 30 d: dE -32.5 (5 d, -13.1 non-telescoping), -71.1 (10 d, -26.6), blow-out at ~18 d; 31.5 s/day in solves.
- VERDICT: the fix removes the energy loss that comes from the iteration (non-telescoping = 0 to round-off in all arms; Bd1 matches the
  converged reference), costs 0.2-0.35x default per day while the star is healthy (Bd1, af1) and <~1.5-2.5x in the late pathological
  phase (tc4), and is bitwise inert when off. It does NOT make the BSG stable in Athena++: the converged discretisation itself develops
  the base wall instability (Bd1 fix and converged Jacobi both die at 2.7 d; af1 fix dies at 9.0 d; tc4 fix survives 30 d but the bath
  drains -31e45 erg by 5 d through his inner BC and the envelope blows out at ~18 d).
