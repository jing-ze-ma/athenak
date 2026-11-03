# alfmach-1102: an Alfven-Mach fix for LHLLD (notes)

Branch `alfmach-1102` from `lowmach2-1102` 592bf7873 (worktree /resnick/groups/carnegie_poc/jingze/wt_alfmach).
Changes are confined to lhlld_mhd.hpp, eos.hpp/.cpp, utils/cb_sensor.hpp (optional component argument)
and two test pgens (gresho_mhd b_phi / cb_ut, wb_atm MHD port), so they cherry-pick onto prod-1102.

## 1. Where LHLLD dissipates a slow vortex (asymptotics)

Scales: flow speed U, sound speed c = U/M, Alfven speed v_A = U/M_A, fast speed
c_f = sqrt(c^2 + v_A^2) = U/M_f with M_f = (M^-2 + M_A^-2)^-1/2 (fast Mach). In the low-M, low-M_A
(incompressible-MHD) limit the physical total-pressure fluctuation is O(rho U^2): the Alfvenic field
perturbation is perpendicular to B0, so B0 . dB (the compressive part) is O(M_f^2)-suppressed and
p_T is the Lagrange multiplier of div u = 0, as p is in hydro.

Writing the HLLD face flux as central + dissipation, per channel:

| channel | where it acts | size of the dissipation | physical scale | ratio |
|---|---|---|---|---|
| outer fast waves S_L, S_R (HLL part) -> velocity-jump term of p_T* | normal momentum, energy (div u) | (rho c_f/2) du_n | rho U^2 | O(1/M_f) (Guillard-Viozat for MHD) |
| same term after the LHLLD fix, phi = chi(2-chi), chi = c_u/c_f | " | rho c_f phi du_n ~ rho max(\|u\|, c_a) du_n | rho U^2 | O(1/M_A) when M_A < 1 (the Alfven FLOOR of c_u) |
| outer fast waves on the transverse variables | v_t, B_t | ~ B_n B_t/(rho (c_f^2 - c_ax^2)) jumps; vanish for B_n = 0 or B_t = 0 | - | O(1) |
| Alfven fan S*_L, S*_R (rotational discontinuities) | transverse momentum and transverse B (Elsasser jumps) | ~ (c_ax*/2) dz, a viscosity/resistivity ~ c_ax dx | flow upwinding \|u\| dx | O(1/M_A) |
| contact S_M | everything advected | \|S_M\| dU | \|u\| dU | O(1) |

* B perpendicular to the vortex plane (B_z, 2-D): B_n = 0 on every face, so there is NO Alfven fan;
  the only O(>1) channel is the p_T* velocity-jump term. LHLLD's chi = c_u/c_f -> c_a/c_f =
  sqrt(2/(2 + gamma beta)): 0.11 at beta 100 (phi 0.21), 0.74 at beta 1 (phi 0.93 = HLLD). This
  is the observed residual 20 % HLLD dissipation at beta 100 and "LHLLD = HLLD" at beta 1.
* B in the plane (aligned B_phi or uniform B_x): both the floor term and the Alfven fan act; the fan
  term is Minoshima's "numerical resistivity max(\|u\|, c_ax) dx" (Minoshima 2024 notes LHLLD leaves it
  at the HLLD value). It is the upwinding of genuine Alfven waves (the characteristic speed of
  incompressible MHD), so removing it makes linear Alfven waves central (as LHLLC/LHLLD already make
  small-amplitude sound/fast waves central); 2dx transverse modes must be protected by a sensor.

What the Alfven floor protects: (1) in M&M 2021/Minoshima 2024 it keeps the pressure dissipation at
the super-Alfvenic (high-beta MRI) scale; (2) at LOW beta c_a/c_f -> 1 so the floor makes LHLLD = HLLD
there, i.e. it switches the low-Mach fix OFF in magnetically dominated gas -- a robustness insurance
for low-beta shocks/blasts; (3) in the lowmach2 tests the floor (vs the face-normal chi) was what
stopped a round-off 2dx checkerboard; with cb_cont that job is now done by the sensor.

## 2. The fix

**A: `<mhd>/lhlld_alfmach = true` (pressure channel).** Scale the Alfven floor of c_u by the local
Alfven Mach number, per side: chi_A = min(1, \|u\|/c_a),
  c_a^2 -> c_a^2 chi_A^2 = min(c_a^2, \|u\|^2),  c_ax^2 -> c_ax^2 chi_A^2.
Then c_u is unchanged (bitwise) where M_A >= 1, and c_u = \|u\| (1 + sqrt(1 - c_ax^2/c_a^2))^1/2
in [\|u\|, sqrt2 \|u\|] where M_A < 1, i.e. chi = O(\|u\|/c_f) = O(M_f) -- the flow-Mach fix of
hydro LHLLC with c -> c_f. cb_cont keeps the 2dx compressive mode at full HLLD dissipation.

**R: `<mhd>/lhlld_alfmach_rot = true` [+ `lhlld_alfmach_rot_min`] (Alfven channel).** On faces whose
Alfven fan straddles the face (S*_L < 0 < S*_R, ** states present), the transverse momentum and
field fluxes become F_ref + psi_A (F_HLLD - F_ref), F_ref = (F*_L + F*_R)/2 (the two outer-star
fluxes, i.e. HLLD without the rotational-discontinuity resolution), psi_A = chi_A(2 - chi_A),
chi_A = min(1, max\|u\|/max c_ax*), floored by the cb sensor applied to the two TRANSVERSE momenta
(full Alfven dissipation of 2dx transverse modes) and by psi_min. The energy flux gets the matching
work/Poynting terms v.dF_m + B.dF_B (as lhlld_emf_eta_pos). Conservative; consistent (D = 0 for
equal states).

**CT/UCT.** LHLLD runs with the GS07 corner EMF, built from the face EMFs (ey, ez = -F_by, F_bz) that
R already modifies, so the induction equation sees the same reduced resistivity. UCT-HLLD
(hlld_uct_mhd.hpp + CornerE_UCT) carries its own HLLD dissipation coefficients and is NOT paired with
LHLLD here (rsolver = lhlld uses CornerE/GS07); if LHLLD is ever run with emf_scheme = uct_hlld, the
UCT d-coefficients would need the same psi_A for consistency -- not done (out of scope).

RESULTS: see section 3 (filled after the runs).

## 3. Results so far (11-03 ~02:00; rest suite and GPU gates still running)
* Gresho, B perp (B_z), M 1e-2, KE kept after 1 turnover, N32/64/128:
  beta 100: HLLD 0.45/0.70/0.94, LHLLD+cb_cont 0.61/0.90/0.98, A 0.86/0.975/0.995.
  beta 1:   HLLD 0.45/0.63/0.91, LHLLD+cb_cont 0.45/0.64/0.92 (= HLLD), A 0.84/0.973/0.995.
* Aligned B_phi Gresho: not a dissipation test; KE GROWS with resolution even for HLLD (2.2x at N128, beta 100):
  physical instability of the current-carrying equilibrium. Uniform in-plane field: field winds up for all solvers (KE 0.28-0.34).
* Linear waves (lw2dm, amp 1e-4): A = 2nd order, errors equal to or below LHLLD/HLLD (fast, Alfven, slow).
  R and A+R: Alfven/fast/slow errors STALL at ~1e-5 (no convergence) -> R FAILS; also blows up on B_phi Gresho beta 1. R rejected.
* Checkerboards (normal un 30 deg, transverse ut along B), beta 100/1, M 1e-2/1e-3: all variants damped to round-off.
* Kriel+26 magnetised Balsara vortex (M 0.01, 1 orbit; 3 orbits too slow on the shared CPU): ME/KE/total kept per orbit:
  64^2 PLM: Luct 88.3/84.7/86.0, A(+R) uct same; ppm4 Luct 88.4/96.1/93.2; GS07 Lgs 90.6/80.2, HLLD GS07 79.0/27.9.
  128^2 PLM: Luct 98.5/96.8/97.4, A uct 98.5/96.7; GS07 (H, Lgs, Lsm, ARsm) 81-98/35-41: the GS07 corner EMF DISRUPTS the
  advected vortex at 128^2 (grid-scale stripes, fig/kq_maps.png); uct_hlld keeps it clean. Kriel: 92.0 (64^2) / 99.6 (128^2).
* CPAW: built-in cpaw pgen crashes at setup (Kokkos view free in the pgen lambda), pre-existing, dropped.

## 4. 11-03 (after the login4 reboot)
* Rest batch (jobs_rest.txt, 93 runs) had FINISHED before the reboot (go_rest.out ALLDONE); only edel/iso_b1 (H, L, A, AR)
  failed and they fail for every solver incl. plain HLLD (non-finite at j = 66, t 29 H/L, t 12 A/AR): the beta-1 isentropic
  wb_atm setup is not a valid test. ana.py now also reduces runs/kb (Kriel beta_k 100 / 1e4). Takeaways:
  - B_z Gresho M 1e-3: A 0.67/0.95 (N32/64) vs L 0.47/0.56, H 0.49/0.46 (beta 100); beta 1: A 0.62/0.94 vs L 0.50/0.47.
    M 1e-4 N32: all ~0.47-0.50 (round-off/eint dominated at N32), A no worse.
  - low-beta Gresho glb_z: A 0.97 vs H/L 0.49; glb_p (B_phi): R KE x59, AR x2.8 (R rejection confirmed), A 0.21 vs H 0.15.
  - Brio-Wu (also beta 2e-3), Orszag-Tang: A = HLLD (dmin/pmin same, no NaN). edel iso_b100: A well balanced (max M 4.4e-7 =
    HLLD), vertical 2dx share 0.06 vs L 0.09, H 0.12.
  - Kriel vortex weak field (beta_k 100 / 1e4): A = L (84.3 / 60.4 % ME per orbit at 64^2); only A+R gains (92.5 / 76.5) --
    in-plane field = Alfven-fan channel, which A does not touch by design.
* Branch: merged prod-1102-r1m 0d85a5539 (r1mspeed) -> aea989762, GPU build athena_gpu_deep_hot_jupiter_rt_alfmach_aea98976_fmadoff
  (md5 b7f098e8, 235 files, 10 min, -j14 on hpc-sm-02-03).
* GPU gates (gpu/gates.txt, hpc-sm-01-11 GPU 3) vs prod-1102 0d85a553: hydro W121 40, r1m v6 40 and 300, v6 with
  lhlld_alfmach=false 40, restart 300 vs 151+149: all rst IDENTICAL.
* r1m try (gpu/runs/try_off, try_on; same new binary, rst trc_1101/keep/off_0p50.rst 8774 -> 13774, 5000 cycles ~0.18 rot,
  hpc-sm-01-11 GPU 3): ms/cycle 47.54 (off) vs 47.58 (on, +0.1 %), dt mean 3.898 vs 3.889 s, RKG stages 15.3 both.
  hst totals on/off-1 at end: KE +1.0e-4 (max 2.3e-4), ME +2.2e-4 (max 4.5e-4), E -2e-10, M 2e-13.
  Pressure-binned (reduce_mhd, on-schedule dumps t 60588/66093/71603): KE, ME differ <= 2e-3 in 1e-6..1e-3 bar,
  <= 3e-4 in 1e-3..1 bar, <= 1.5e-4 at 1..1e3 bar; radial 2dx v_r in 1e-6..1e-3 bar -15 % / +2 % / +8 % (sign flips = chaotic
  divergence, not systematic). -> lhlld_alfmach is essentially inactive in r1m (flow is super-Alfvenic where it matters).
  Recommendation: do NOT add to r1m (no gain, opt-in key for slow strong-field problems only).
