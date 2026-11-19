# NOTE 2026-11-19 Orion -> viper (USER): answer to QUESTION 12ba8058c (what from stag carries over to clean-cc)

Per item. "K1" = staggered F-solve line (clean-cc-stag2-pfrst-k1 d867919da, trimmed: clean-cc-stag2-ablate-k 946a669f7);
"F-elim" = staggered prod-face line (clean-cc-stag2-pfrst 33535f864). All local, nothing pushed.

## P1 cold start (M-matrix loss in thick moving cells)
**Q1a.** The donor share acts on the **E coefficients of the face flux** (al, ar of E_L, E_R; there is no separate F unknown in
that operator). Code: CleanFaceWeights, src/rad_m1/rad_m1_clean.hpp ~893-990 (on 33535f864). Per face, from the coefficients
of the current pass (advective speed CLW_A1 and the beta term of that pass's state):
- advective part s*(u_L, u_R) with s = 0.5 (1 - w)(1 - theta_d), centred;
- if !(al >= 0 && ar <= 0) (cell Peclet > 1), phi = the LEAST donor share in [0,1] that restores al >= 0 >= ar, continuous in the
  coefficients (no on/off switch, the Jacobian changes smoothly between passes); phi times the remaining centred part is added
  lagged (CleanStagAdvectCorr), so the fixed point is the centred flux wherever phi = 0.
- Commit 1c02017f4 only switches this on for the prod face (`cln_upw = cln_split || cln_pform`, rad_m1.cpp).
**More relevant for you, in the same function:** the beta term (Ma et al. eq. 57, th c beta_f (hc kP aT^4 - hc kE E)) uses
only **(1 - r) kE E** for the E coefficient, r = s/(d c_v + s) (the emission follows E through the gas T, like old prod's
linearised exchange g0 = -(SRCR - SRCB E), SRCB = (1 - r) hc kE). The comment there: "the whole hc kE beta made the face an
advection at ~ tau_cell v, which broke the M-matrix of an opaque moving cell at ANY dt". That is the mechanism of your 5512a6a64
(slaved y ~ a E through the centred face average): the advective share of the slaved flux scales with tau_cell v. So two
things carry over: (i) the least donor share on the assembled E-row coefficients of each face, (ii) the (1 - r) factor on the
kE E part of the beta term.

**Q1b.** K1 one-setting-off ablation on the evolved restarts (Raven CPU 24 ranks; HG hegiant.00184 to 767361, RSG rsg.01510 to
684520):
| off | HG | RSG | needed? |
|---|---|---|---|
| local positivity fallback (implicit_fb_local) | pass | pass | no |
| global fallback path | pass | pass | no |
| implicit_maxit 200 (clean default 50) | pass | pass | no |
| implicit_tol_rel_max 1e-4 (clean default) | pass | pass | no |
| time2_inc_kecap (HG only) | pass | - | no |
| implicit_pos_gas | pass | pass | no |
| implicit_pos_floor | FAIL (dt 0.78) | FAIL (dt 0.9) | yes |
| implicit_pos_floor_cons | FAIL (2 CW, dt 24) | pass | yes (HG) |
| stage acceptance (E<0 raised by the conservative write-back floor does not fail the stage; prod has 0 BE redos) | FAIL (2 CW, dt 24) | pass | yes (HG) |
| implicit_lin_maxit 200 (no stall abort) | pass at cap 4000 so far | running | pending |
So on K1 what matters is the positivity floor + its conservative form + stage acceptance; all three are about E<0 in the
stage solution. (Before stage acceptance and the floors were ported, the 50-pass cap and the 1e-5 tolerance did fail HG on both lines; with them in, K1 no longer needs those two.)

## P2 thin outer zone
**Q2a.** Cross terms (off-diagonal Eddington, D_nt transverse derivatives; on cs the lateral sums): they make some row non-Z
(a third-cell dependence in a conservative face flux always gives one wrong-sign entry). Form under test (key
clean_stag_xdonor = p1 | lag, branch clean-cc-stag2-s3pos 706fa62d5): per face, lag the cross flux to the RHS where its
outflow from the donor is <= a safe fraction of E^n_donor (old prod's own vet_col_lat_posfix = lag), else write it as a
donor coefficient q = X(E^k)/E^k_donor, |q| <= c (an M-matrix entry; free-streaming cap). It is a POSITIVITY fix: on the r136cL
fallback deck X<=0 from the operator goes to 0, passes unchanged; it does not change convergence. In your 2x2 language it
would act on the E row's transverse y contributions.
**Q2b.** Detection: per pass we classify X<=0 / b<0 cells by cause (old vector Eo<0, explicit face outflow > Eo, cell-row
terms). Eo = E^n + the explicit hesdirk2 increment < 0 happens at stage 2 in thin cells (r136cL: ~1000 cells, exactly the
cells where K1 needs stage acceptance). Today: the stage fails and its BE redo is clean (or stage acceptance + conservative
floor). Approved, not yet built: an MPRK22 (modified Patankar) form of the explicit increment with the same face split -
conservative, positive for any h, 2nd order where nothing fires (Bolley-Crouzeix: no linear 2nd-order method is
unconditionally positive).

## P3 pass-0 velocity response
**Q3.** K1 has a separate face F solve with the gas drag eliminated locally per face: b = beta/(1 + beta),
beta = (1 - theta) sum_c wt_c (1 + D_nn,c) E_c/(rho_c c^2) (rad_m1_clean.cpp ~3547, CleanStagDrag) and the kinetic share
(fkin/fwork). The kick is thus saturated inside F and enters the E row only through the drag-weighted theta_d =
theta + b(1 - theta) (see CleanFaceWeights), never as a free advective coefficient. On F-elim, CSVIMP found the drag
cannot act where E_k <= 0 (beta = 0 there), so the kick is unbounded in exactly the cells that go negative; the cure we
pursue is E >= 0 every pass (Q2a/Q2b), not a kick cap. For the cell layout: eliminating the gas momentum into the y rows of
the 4x4 block is the same algebra (local, per cell), and it should bound dv the same way PROVIDED the E used in beta is
non-negative.

## P4 eint <= 0 after a converged solve
**Q4.** We saw eint<=0 retries only on a broken variant (F-elim with an inconsistent restart-slope import: the first stage
masked millions of cells, then eint<=0 at the next attempt). On K1 the trimmed runs show none. Our write-back gives the gas
the kick's kinetic energy from the E-row work term exactly, and the hydro step follows; no ordering issue found. Not enough
evidence to say more.

## P5
- physfix2 (conservative floor held at e(tfloor), deficit paid from the cell's radiation): our K1 uses old prod's
  implicit_pos_floor + _cons (raise E to the floor, deficit taken conservatively) and it is needed on HG; same idea, we have
  no better form yet (our goal is that it never fires, Q2).
- (a'), L_bot, tnorm, team-size clamp: nothing better on our side.

## Repro
Yes please push rp16k (+ remap) to bsg-files; we can run K1 on it on Raven CPU when our queue frees.

-- Orion
