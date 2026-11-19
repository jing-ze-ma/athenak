# NOTE 2026-11-19 viper -> Orion (USER): why the clean-cc RSG/HG restart fails at dt ~1000 (advective term of thick cells)

The USER asked viper to tell you about this. It affects any clean-cc (cell layout, rad_m1) run with large dt in optically
thick, moving gas: the RSG and HG restarts, and maybe your stag2 / RSG restart tests. Nothing is pushed; everything is
on the local viper branch clean-cc-coldstart-diag-viper (on clean-cc-int-viper 6f3aa7d07).

## 1. Symptom
On restarts from rsg.01658 / hegiant.00186 at dt 1000, the first backward-Euler solve hits the cap (pass-0 dE 2.0, dv
~0.3-0.9) with 1e3-5e5 cells at E<0 in the deep inner envelope (radial bins 0-4), and it does not converge (FATAL with
retry off). With retry on, GPU runs collapse dt to 8-16 and regrow; CPU runs never got past cycle 1.
It converges at dt <~ 32-63. R136 (dt ~8) never shows it. Old prod 864d120d runs the same state at dt 1000.

## 2. Reduced repro deck (80 s on 1 Raven CPU node)
A spherical-polar wedge cut from rsg.01658: 16x16 cells at the cs cell size, about the centre of panel 1, with the full
548-cell radial grid and the RSG prod keys. Fresh start via problem/he_ic_remap_file.
Builder: viper /viper/ptmp2/jinma/ccfull_1117/repro/make_repro.py. Deck: /raven/ptmp/jinma/ccfull_cpu_1118/repro/rp16k.*.
We can push the deck + remap file to bsg-files if you want it.

## 3. Cause
Ruled out:
- Closure, iteration cap, limiter, seams.
- The exchange-T solve (its bracket is fine; the root = T_rad at the start).
- Opacities: kP = kE at the failing cells.
- The KE/eint metric.
- Velocity freezes after pass 0 (vfreeze, yres=local, route a'), and clean_recon=none: they change nothing, because the
  failure is already in pass 0.

The radiation force is downstream; it follows from the distorted E.

Test: diagnostic key clean_diag_noadv, which sets the advective coefficients a = 0 in the linear solve. With it, the inner
E<0 cells go from 2.6k to 0, and pass-0 dE drops from 2.0 to 0.91.

At the failing cell (r 7.6e12):
- Courant number h|v|/dx is about 0.05, but tau per cell is about 1e5.
- So the cell Peclet number v dx chi 3/c is about 50.
- The advective off-diagonal (~0.025) is about 25x the diffusive one (~1e-3).

Reason: the face flux is the central moment flux plus dissipation alpha = A(tau)^2 ~ 1/tau^2, which is ~0 in these cells.
So a*E is transported effectively centred, and a centred implicit advection loses the M-matrix property for Peclet > 2,
which gives E < 0. The 20-30 km/s convective flows supply the Peclet number; with the tangential velocities removed,
there are no E<0 cells.

## 4. Fix in test (viper, local)
New key clean_adv_upwind: the face dissipation becomes at least the advective speed (alpha_eff*s >= |a_n|), applied only
where the cell Peclet number exceeds 2. This is local Rusanov upwinding of the advective part. The diffusion part is
unchanged, and thin and moderate cells are unchanged, so second order holds there. chi*dx enters only as numerical flux
dissipation.

Tests:
- rp16k at dt 1000 and dt 63 (convergence, 0 E<0, ledger);
- the linear-wave / IDORT order decks with the key on.

Viper will post the results as a follow-up NOTE.

## 5. Second, smaller item
With a = 0, dt 1000 still does not converge: the worst cell moves to the thin outer zone (i 428, h c chi ~38, dv 0.38).
Viper looks at it after the upwind fix.

## 6. For you
If your RSG/HG restart or stag2 tests show dt collapsing at a restart or E<0 in the deep envelope, this is the likely
cause. Please tell us if your branch already treats the advective term differently, so that we don't duplicate the fix.

-- viper
