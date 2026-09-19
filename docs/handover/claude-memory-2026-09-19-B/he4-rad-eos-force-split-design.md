---
name: he4-rad-eos-force-split-design
description: design agreed with the user 2026-09-19 night - split the radiation taper's two jobs: EOS weight density-only with a very low window, radiative FORCE gated by optical depth inside the two-stream
metadata:
  type: project
---

**Problem.** The EOS radiation weight w(rho,T) = max(w_rho, w_T) (rad_taper::WeightGated; window rho 4.97e-10..1.93e-9 = tau 0.3..3 of the IC, T gate 4.52e4..6.33e4 K) does two jobs: (a) thermodynamics (e = e_gas + w a T^4, p = p_gas + w a T^4/3), (b) deciding where the force is -grad(P_rad) (w=1) versus the explicit (1-w) rho kappa F/c + P_rad grad w in two_stream_rt.hpp (~line 4776). On the He4 star the density at tau~1 is only 3x below interior densities (CZ 2.3e-9, thins to 6e-10 in the pulsation), so a density window cannot do both; the T gate was the patch, and it wrongly radiation-loads shock-heated transparent gas above the photosphere (tall-domain arm b_tall: c_s hits cs_max, dt 0.6 s). Density-only with the present window fails thermodynamically: w drops at fixed e -> T jumps (e_rad ~50 e_gas) -> evacuation loop (tests_r6).

**Design (user: "we probably need to do this").**
1. EOS: density-only weight, window far below anything stellar: ~1e-11..1e-10 (eos_rad_rho_lo/hi), T gate OFF (eos_rad_t_hi = 0). LTE P_rad = aT^4/3 is fine down to tau~0.1 (Eddington approximation).
2. FORCE, in the two-stream where tau per cell exists: new default-off switch (e.g. problem/rt_force_tau_gate, rt_force_tau_lo=0.3, rt_force_tau_hi=3): s(tau) = 1 for tau<=lo, 0 for tau>=hi, cubic in log tau;
   f1 = s*[rho kappa F/c + d/dr(w a T^4/3)] + (1-s)*[(1-w) rho kappa F/c + P_rad dw/dr]
   (s=0 reproduces today's formula exactly -> bitwise when off); transverse: f2,3 = s*d_t(w P_rad) + (1-s)*P_rad d_t w. Needs w at i+-1, j+-1, k+-1 (WeightGated on rhoN, T_g: T_g is ghost-filled one angular cell wide since de1d41e4) and tau-to-top in the force kernel (a `tautop` exists in the cell-report code: check it is available there).
3. IC: the file stores (rho, eint) with the OLD w. Regenerate eint above ~0.97 R with the NEW w at the same (rho, T) (offline EOS in tests_r12/relax_ic.py: temp_of with old taper_w, then e_gas + w_new a T^4). Gas-pressure balance is unchanged at t=0 IF s(r) = 1 - w_old(r): old net radiative force = -w grad P + (1-w) kappa rho F/c, new = -(1-s) grad P + s kappa rho F/c. The w thresholds were read off tau = 3 and 0.3, so s(tau) with the same two levels is close; verify with tests_r12/fbud.py (per-cell force at t=0, top 5 % of the radius) on the 1-D sp column (tests_r13/r13_arm.sh).
4. Check the EOS table inversion brackets (eos_table.cpp ~300-400 assume "gate => hot cell has w=1"; with the gate off that branch is unused) and that e(T) stays monotone.
5. Then the tall domain (tests_r13: fitgrid.py, mk_tall_ic.py, floor atmosphere rho 3e-13) should no longer need cs_max to tame the gas above the photosphere.

**Why:** makes the region above the photosphere (wind launching, Gamma ~0.83-1) physically sane so x1max can go to ~1.25 R and the pulsation stops draining the envelope through a top at 1.014 R. See [[he4-r11-floor-arms-verdict]] for the evidence.
**How to apply:** implement default-off, bitwise check vs athena_v19, 1-D validation before any 3-D.
