# m1-positivity (09-30): positivity of the implicit M1 solve in the He-star wedge plume

Branch `m1-positivity` (off `he-presn-m1` abf547b5), worktree /viper/ptmp2/jinma/wt_m1pos,
run tree /viper/ptmp2/jinma/m1pos_0930 (diag*/ iso*/ smoke/ gate_off/ runs/ hebox/).
Binary (He star, HIP ROCm 7.2): /viper/ptmp2/jinma/m1pos_0930/bin/athena_gpu_he_star_m1_7192a176.

## 1. What fires first (instrumented, read-only: dbg_t2_admiss now also reports the
implicit_vimp fallback and every solve that ends with E <= e_floor or eint <= 0)

ad3d_128 rst 00011 (t = 25850 s, 5.5 turnovers), keys of ad3d_128/r2.sh, 13 cycles (diag0/diag1):
- every hesdirk2 stage (13/13) is redone with backward Euler. The stage is CONVERGED and its
  final E, T are positive (min(E,T) = 0.2177); it is rejected only because the implicit_vimp
  positivity fallback fired on the FIRST Picard pass (vimp_now != impl_vimp).
- that pass-0 linear solve has E < 0 in 281 + 91 cells (min -2.28e8 against E ~ 7e6), in dense
  plume cells (rho/<rho> 1-10, r = 1.46-1.65e11, v_r 1-3e7). Their rows are Z-matrices (0 positive
  off-diagonals, vimp block ~1e-2 of the diagonal); the right-hand side is the problem:
  x1 face term -c dt v_f g0_f = -1.9e8 (base 2.3e7), opacity-Newton rhs +-1e7..9e7.
- g0 = rho (kappa_E E - kappa_P a T^4) of the ITERATE: in these radiation-dominated, stiffly
  coupled cells (c dt rho kappa_P ~ 5e5, G0 ~ 30, chat dt g0 ~ 4e12 >> E) it is the difference
  of two numbers ~1e5 x the converged exchange; at pass 0 it is garbage.
- then (g0 fixed): the opacity-Newton NEIGHBOUR entries keep the right sign but break row
  dominance (sum|off|/diag = 2.4 at r = 1.46e11, 2.0 at the photosphere r = 2.39e11): E < 0.
- then (both fixed, vimp kept): Picard limit cycle at res ~1.5e-7 > tol 1e-8 (iso5/E plog),
  6/8 solves NON-CONVERGED -> stage rejected anyway.
- gas eint <= 0 appears in the OLD vector (after hydro) in 37 near-void cells (rho/<rho> down
  to 2e-5); the radiation write-back makes them positive in all but 1 of 26 solves.
- The transverse operator is NOT the trigger: advection is already upwinded in x2/x3
  (M1SphTransRow / m1_impl_tcell) and every failing row had 0 positive off-diagonals.
  Item 1 of the brief (ap_hll on x2/x3 faces) was therefore not built.

## 2. New keys (all read only when named: absent = bitwise the old code, restart echo too)
- `implicit_g0_exchange = true`: from pass 1 on g0 = -(SRCR - SRCB E^k)/(chat dt), the
  exchange of the previous pass's linearised source row (identical to the pointwise g0 at the
  Picard fixed point).
- `implicit_g0_limit = w` (1.0): pass 0 only, |chat dt g0| <= w (max(E^k, E_old) + (chat/c)
  max(e_gas,0)). Clipping every pass was tried and REJECTED (1.7e7 cell-passes clipped, 5/6 NC).
- `implicit_opac_newton_guard_mode = 6` (bit 4 new): a row's remaining Newton face terms are
  dropped when they raise sum|off|/diag above max(1, value without them) x (1 + guard).
- `implicit_pos_gas = true` (+ `implicit_pos_gas_frac`, 1e-3): written-back gas eint below
  max(e(rho, tfloor), frac x max(e_gas old, 0)) is raised with energy taken from E (down to
  e_floor); e_gas + (c/chat) E per cell exact; counted (cells, energy).
- `implicit_pos_floor = true`: E raised to e_floor is paid by the gas (down to its limit);
  the uncovered part is counted as "energy created".
Counters: one `<rad_m1> m1-positivity:` line in the end-of-run report.

## 3. Short tests from rst 00011 (4 cycles, apudev, binary pos4 = 7192a176)
| arm | keys | Picard mean / NC | stage fallbacks | vimp fallbacks | E<=floor / eint<=0 |
|---|---|---|---|---|---|
| base (diag0, 13 cyc) | none | 19.5 / 2 of 26 | 13/13 | 26 (min E -9.7e8) | 31 / 1 |
| iso4 A | exchange | 16.7 / 0 | 3/3 | 6 (-2.3e8, pass 0) | 0 / 0 |
| iso5 G | all + guard 6, vimp on | 156 / 6 of 8 | 4/4 (NC) | 0 (min E +2086) | 0 / 0 |
| iso6 H | vimp off only | 49 / 2 of 9 | 3 | - | 18 / 4 |
| iso6 J | all + guard 6, vimp off | 28 / 0 of 8 | 1 of 4 | - | 2 / 3 (both limited, 0 created) |

## 4. Gates
