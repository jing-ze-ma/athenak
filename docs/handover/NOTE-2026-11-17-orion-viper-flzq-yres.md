# NOTE 2026-11-17 Orion -> viper (USER): which clean_yres is right on flz_q?

Reply to TASK-2026-11-17-viper-orion-flzq-yres.md.

## Verdict
- **Neither answer is "the right one" at the harness resolution, and the gap is not an error in local's force law.**
- In one step from the same state, false, post and local give the same gas momentum to 2e-14 (backward Euler), and to 7e-7 with
  hesdirk2. That holds with the harness tolerance and with a tolerance 100x tighter.
- The O(1) end-state gap comes from a 5e-9 seed that the run amplifies by about 10x every 2-3 steps. The amplifiers are the discrete
  branches of the scheme:
  - the positivity face marking;
  - the hesdirk2 stage fallback to backward Euler.
- Other tiny changes are amplified the same way:
  - the tighter tolerance alone flips false at 2x (L_rz 17 %);
  - gas_rot x (1 + 1e-6) gives 1.4e-2 of max E at 1x;
  - halving cfl changes the 1x answer by 0.5-0.8 of max E.
- So the flz_q end state is not converged in dt or dx at the 1e-1 level. The 3 % L_rz and 0.9 max|dE| of the regression are inside
  that spread.
- **Nothing here argues against local as the push-18 default.** Its results conserve L_z and E to round-off, as false does, and
  restart bitwise.
  - Weak point in local's favour: across 4 matched perturbation pairs, the 1x local runs sit closer to both 2x solutions (EL1 0.05-0.07
    for local, 0.12-0.13 for false; 0.10 vs 0.18 at cfl 0.15).
  - Weak point against local: at 1x it is more sensitive to clean_tol. Going from 1e-14 to 1e-12 moves local by 0.68 of max E and false
    by 4e-4. The likely cause is local's pass-skip rule `inl == 0 || nit != 0`. At 2x, false is just as sensitive (see above).
- **flz_q's E/m1_f should not be a bitwise or tight gate for yres-type changes.** Gate it on L_z and E conservation, the iteration
  counts and the L_rz/L_eq ratio.

## Setup
- Binary: build_p17c = b66ff2552 (md5 90608574d0593dbaa47033327e950cc6). It differs from fork clean-cc 5785ee51d only in the clean_yres
  default (`git diff --stat b66ff2552 5785ee51d`: rad_m1.cpp/.hpp, 4+/7-), so I wrote clean_yres explicitly into every deck.
- Deck: the harness flz_q (fb.py flz, sp fv1, 32x8x8, chi = 1, rho 1e2, T_gas 1e-4, gas_rot 2e-3, vet_gd GL16x32, hesdirk2).
  - Changes from the harness: `tlim = 32.8` (= the 40-step time) with nlim open, so all arms stop at the same t.
  - Keys in the decks: `clean_yres`, `clean_tol` / `clean_tol_nl` (default 1e-14/1e-13, tight 1e-16/1e-15). I removed the deck's own
    duplicate tol lines.
  - Resolutions: 1x 32x8x8, 2x 64x16x8 (r and theta).
- 4x (128x32x8) is too slow for this window: t = 6.6 of 32.8 after 1.3 h at 40 threads. Up to t = 6.6, false and local agree to
  2.2e-5 in L_rz. The other 4x arms were stopped.
- The 2x tight arms reached t = 30.2 (false) and 25 (local) before I stopped them. post_t_x2 got few CPUs and is not used.
- All runs: srun --overlap + pick.sh in 230761 / 230775. 0 CLEAN WARNING in every finished run. Restart dumps deleted.
- Files: /orion/ptmp/jinma/clean_1115/flzq_yres (ana.py, ens.py, force.py, lrz.py; runs under reg/fb/runsf/flz and ex/, f1/, rst/).

## 1. Arms at t = 32.8 (E, F vs the 2x solutions, volume-restricted onto 1x)
| arm | steps | stage fb | L_rz | L_rz / L_eq | Emax / EL1 vs false_x2 | Emax / EL1 vs local_x2 | inner-wall v_r |
|---|---|---|---|---|---|---|---|
| false_d x1 (harness) | 40 | 3 | 1.0862e-3 | 0.998 | 1.04 / 0.132 | 0.91 / 0.118 | -4.2e-4 |
| false_t x1 | 40 | 3 | 1.0862e-3 | 0.998 | 1.04 / 0.132 | 0.91 / 0.118 | -4.2e-4 |
| post_d / post_t x1 | 40 | 3 | 1.0862e-3 | 0.998 | 1.04 / 0.132 | 0.91 / 0.118 | -4.2e-4 |
| local_d x1 | 40 | 2 | 1.0507e-3 | 0.987 | 0.65 / 0.069 | 0.59 / 0.054 | -4.2e-4 |
| local_t x1 | 40 | 2 | 1.0507e-3 | 0.987 | 0.66 / 0.069 | 0.60 / 0.055 | -4.2e-4 |
| false_d x2 | 89 | 15 | 1.0684e-3 | 0.979 | 0 | 0.74 / 0.050 | +1.2e-4 |
| local_d x2 | 90 | 13 | 1.0756e-3 | 0.981 | 0.51 / 0.049 | 0 | +1.2e-4 |
| post_d x2 | 89 | 15 | 1.0684e-3 | 0.979 | 2.5e-6 / 1e-7 | 0.74 / 0.050 | +1.2e-4 |

Same-grid pairs and tolerance:
| pair | Emax | EL1 | FL1 |
|---|---|---|---|
| post vs false, x1 | 2.7e-8 | 6.6e-9 | |
| post vs false, x2 | 2.5e-6 | 1.1e-7 | |
| local vs false, x1 | 0.87 | 0.155 | 0.88 |
| local vs false, x2 | 0.51 | 0.049 | 0.75 |
| false d vs t, x1 | 6.8e-6 | | |
| local d vs t, x1 | 8.1e-3 | | |

- **No analytic L_rz in the deck.** As the expected answer I use the drag equilibrium L_eq = sum V R (4/3) E v_phi / c^2 (P = E/3),
  since t_drag = 1/(c chi) = 0.1 << dt.
  - Every arm sits at it to 1-2 % (above). The 3 % L_rz gap is false and local at equilibrium with slightly different E fields.
  - L_rz is not smooth in t. It jumps 5-10 % in single steps in every arm (positivity marking events). At t = 2.5 it is 0.79e-3 at 1x
    and 1.14e-3 at 2x.
- **Wall v_r: false and local give the same sign at each resolution.**
  - Inner wall: inward at 1x (-4.2e-4), outward at 2x (+1.2e-4), outward at 1x with cfl 0.15 (+8.6e-4).
  - So "inward" is not a converged property of this deck, and it does not separate the yres variants.
  - Outer wall: +1.5e-3 (false) and +1.2e-3 (local) at 1x; +9.8e-4 and +7.6e-4 at 2x.
- Where they differ: E differs O(1) only in the outer-wall cells (i = 31). There, E has a theta-sloshing pattern with opposite phase in
  the two arms. The bulk mean E differs by up to 8 %, with total E conserved in both.
- Convergence toward each other:
  - Tolerance does not bring false and local together. local_t vs false_t is 0.877, the same as at the default tolerance.
  - Resolution narrows the gap a little: EL1 0.155 at 1x, 0.049 at 2x.
  - post equals false at every resolution (1e-8 at 1x, 1e-6 at 2x).

## 2. Where the split starts (1x, L_rz and Er rel. difference local vs false)
| step | 1-6 | 7 | 10 | 14 | 17 | 22 | 25 | 28 | 29 | 31 |
|---|---|---|---|---|---|---|---|---|---|---|
| L_rz | <2e-12 | 5.5e-9 | 1.3e-8 | 1e-7 | 6e-6 | 4.6e-5 | 1.3e-4 | 3.7e-3 | 3.5e-2 | 7.7e-2 |

- The CLEAN lines are identical in pass/mark/neg until cycle 27, where the marked faces differ (32 vs 24).
- Cycle 32: false's stage 3 ends with neg 40 after maxpos and falls back to BE, then takes an auto-BE step. local accepts the stage.
- At cycle 31 (before that fallback) E already differs by 0.94 of max E.
- At 2x the split (> 1e-3 in L_rz) starts at t = 18.0 for local vs false, and at t = 26.6 for false tight vs false default.
- Other 1x perturbations at t = 32.8, vs the harness false run (max of E):

| perturbation | Emax |
|---|---|
| tol 1e-12 | 4.3e-4 (false), 0.68 (local) |
| gas_rot x (1 + 1e-6) | 1.4e-2 |
| cfl 0.15 | 0.52 (false), 0.83 (local) |

## 3. Force balance per cell (one step, cycle 20 -> 21, all arms restarted from the same local 1x restart)
- dm_rad = m(variant) - m(gas_feedback = false), with the same step. The reference is the tight-tolerance solve.

| time_scheme | arm | vs tight solve: max / L1 (of max\|dm_rad\| = 2.3e-2) | vs h chi (F - (4/3) E v)/c: max / L1 |
|---|---|---|---|
| be | false_d / post_d / local_d | 2.0e-14 / 1.2e-14, 2.0e-14 / 1.1e-14, 8.3e-14 / 2.4e-14 | 0.098 / 0.046 (all) |
| be | false_t vs local_t | 8.2e-14 | |
| hesdirk2 | false_d / post_d | 2.3e-14 | (stage combination, not one BE drag) |
| hesdirk2 | local_d and local_t | 7.3e-7 / 5.2e-7 | |

- At solver tolerance, the linear residual that false gives to the gas and local does not is about 1e-14 per step. So **local's
  drag-only share and false's minus-the-radiation-change are the same force to round-off.** Neither is closer to h chi (F - A)/c.
- The 4.6 % L1 against the Eddington estimate is the P = E/3 approximation, and it is the same for every variant to 3 digits.
- The one systematic difference is 7e-7 per hesdirk2 step, independent of the tolerance. It is local's stage handling (likely the pass
  skip `inl == 0 || nit != 0` in rad_m1_clean.cpp:726-727, or the stage weights on cln_dy), not drag physics.
- The chaotic growth turns this 1e-7 to 1e-9 seed into the O(1) end-state difference.

## 4. Restart (local, 1x, harness deck with clean_yres = local written in)
- 40 cycles straight vs 20 + restart + 20: hydro_u and m1 bins at cycle 40 have max diff 0.0, and the user.hst last 20 rows are
  identical. **Bitwise.**
- The "no implicit_predictor increment" restart warning is the known stale text.

## Suggestions
- If viper wants a converged flz_q reference, run 4x (and 2x tight) on GPU to t = 32.8. On CPU it takes more than 6 h.
- A shorter tlim (t <= 10, before the marking-driven growth) would make flz_q a meaningful E/F regression gate.
