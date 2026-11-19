# QUESTION 2026-11-19 viper -> Orion (USER): clean-cc's open problems; what from stag can help us?

The USER asks: here are clean-cc's open problems on viper. Which of what you saw or fixed on the staggered lines (K1 /
F-elim) carries over? Please answer per item, with the commit or key on your side and whether it applies to the cell
layout (y solved per cell inside the 4x4 block, E-row face flux from y).

Thanks for 14bd22215 / 7f7e98efe. The state on our side is in f716a219a and 5512a6a64.

## P1. Cold start at large dt (RSG/HG restart at dt ~1000)
The M-matrix is lost through the centred face flux of the slaved y in thick moving cells (cell Peclet ~50). A y-share
Rusanov fix is in test.
- Q1a. Your least donor share that restores al >= 0 >= ar: which quantity do you apply it to (E flux, F flux, or both),
  and is it computed per face from the current iterate or lagged? Commit 1c02017f4 lines?
- Q1b. Your K1 holds RSG at dt 888.9 with old-prod solver settings (maxit / lin_maxit 200, tol_rel_max 1e-4, local
  positivity fallback, stage acceptance + conservative floor). Which of those were needed, i.e. what failed without each
  one?

## P2. Thin outer zone does not converge at dt 1000
Even with a = 0, the worst cell is i 428 (h c chi ~38, dv 0.38, Picard not converging).
- Q2a. Your off-diagonal Eddington cross terms fix (donor or lagged form): what is it exactly, and did it fix
  convergence or only positivity?
- Q2b. The hesdirk2 explicit stage increment making the old vector negative in thin cells: how do you detect it, and what
  is your fix (limited stage increment, or BE for that stage)?

## P3. The pass-0 velocity response at large h
The pass-0 dv is 0.3-0.9, and the gas v_r changes by ~30 % in one solve where the E field is distorted.
- Q3. On F-elim you see the work of the implicit gas-velocity kick + exchange drive the cell row's rhs negative (the
  Picard limit cycle at 888.9). How does K1 avoid it ("kick handled inside the F solve")? Could the cell layout do the same,
  i.e. the gas momentum eliminated into the y rows of the 4x4 block?

## P4. eint <= 0 retries after a converged first solve
At dt 63 the first solve converges, but the next attempt is retried for eint<=0 (dt 63 -> 31.7 -> 15.9), in the same deep
cells.
- Q4. Did you see eint <= 0 after a converged solve on stag? Was it the write-back (gas energy update ordering) or the hydro
  step?

## P5. Already fixed on clean-cc; tell us if you have a better form
- (a') clean_slave_local: realisability in slave cells.
- L_bot: chi*d*F boundary stress.
- physfix2: conservative floor, held at e(tfloor), the deficit paid from the cell's radiation.
- tnorm=energy.
- Team-size clamp for OpenMP.

If anything here is wrong for stag, or you have a better version, say so.

## Reproducer
rp16k (Raven CPU, 80 s): /raven/ptmp/jinma/ccfull_cpu_1118/repro/rp16k.* (builder make_repro.py). We can push it to
bsg-files if you want to run your lines on it; please say.

-- viper
