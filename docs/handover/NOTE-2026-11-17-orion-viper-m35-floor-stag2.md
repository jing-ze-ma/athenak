# NOTE 2026-11-17 Orion -> viper: M35MMS floor runaway: would stag2 avoid it? (part 1 of TASK-2026-11-17-viper-orion-m35-floor-and-cs-corner)

Short answer: **yes, by construction for the gas momentum and the gas energy.** The staggered layout gives the gas the
physical force and the exchange evaluated on the solution. The linear-solve residual is put into E, never divided by rho.
A 1-D sp floor test agrees. Staggered results are identical at clean_tol 1e-10 and 1e-5. Production conservative2 pumps
the floor and collapses dt at loose tolerance; matched does not.

Two caveats for the cell layout (not stag2):
- With clean_gas_exact = false, the residual goes into eint. The test reproduces this: floor T x8, dt /15.
- The cell momentum is still the conservative "minus the radiation change" of the BiCGStab y unless clean_yres = true
  (default false). This route was not excited in the test.

Code: clean-cc-stag2-gpu2, e9227103a (Orion CPU build). Line numbers are for that tree; the head 0e959f97c only adds 7
lines to rad_m1_clean_solver.cpp.

## 1. How the gas gets its momentum and energy

| | staggered (stag2 defaults) | cell layout | production conservative2 |
|---|---|---|---|
| unknowns of the linear solve | E only; face F is an explicit function of E, plus its own exact line (Thomas) F solve `CleanStagFSolve` rad_m1_clean.cpp:1697, stop at 1e-14 of max\|F\| (:1867) | E and y = F/c per cell (4x4 blocks, BiCGStab) | E and F per cell |
| gas momentum | `CleanMomentumChange` rad_m1_clean.hpp:945: dm = -wt [(F - G)/c^2 + h d_n P - h diss]. Since F satisfies its face equation to round-off, this **is** h chi_f (F - A)/c (+ beta terms) shared in proportion to chi_c V_c (:925-:932). So dv = h kappa (F - A)/c, with no 1/rho, and \|F\| <= c E_donor (R2, al = 1 on violators). The consistency bracket is on the energy flux only, not the momentum F (rad_m1_clean.hpp:772, `!x.hi`), so it gives the gas no force | rad_m1_clean_cell.cpp:1860 (write-back) and :1514 (iterate): dm = -(y - y_old + h div fY)/c from the BiCGStab y. This is conservative, **the y-row residual reaches the gas as residual/rho**. clean_yres = true (:1853-:1857) first removes R_y/(1 + h c chi); default false (rad_m1.cpp:250) | dm = -(F'_c - F^old_c)/(c chat) - sum Z_f/2: residual and mismatch divided by rho |
| implicit drag | `CleanStagDrag` :2185, b = beta/(1+beta) (:2216), beta ~ (1-theta)(1+D)E/(rho c^2): bounded for any E/(rho c^2) | Newton on y | vimp prediction (cs2_vimp) |
| gas energy | write-back :2318-:2338. The conservative den is formed first. Then rd = (den - EP + ex - wk)/(1 + h c kE (1-r)) is **moved to E** (:2333-:2335). The gas keeps -ex(T_k, E) + beta terms + W: the exchange on the solution at the iterate's T. In the source comment: "the linear solve's residual goes to E, not eint" | clean_gas_exact = true (default, rad_m1.cpp:249): eint is the root of the cell's exact local exchange with the face fluxes fixed (rad_m1_clean_cell.cpp:1891-:1911), and E = C - eint takes the rest. Off: den = -(dE) - div ... goes to eint unchanged (only a flag if eint < 0, :1912) | conservative "minus radiation change" |
| W (KE gained) | wk = (p + dm/2).dm/rho (:2302) of the physical dm, paid by E | same form (:1876), but of the conservative dm, so a residual kick in dm becomes a large E debit | - |

## 2. Other routes, checked

| route | staggered | cell |
|---|---|---|
| F-solve residual | direct tridiagonal per line, iterated only for the lagged MeshBlock line ends, to 1e-14 relative (:1867; f_max floored by 1e-8 c\|E\|) -> dm error ~1e-14 \|F\|/c^2 | the y rows are in the global BiCGStab: \|\|r\|\| <= clean_tol \|\|b\|\| (rad_m1_clean_solver.cpp:718), a global norm, so a floor cell can hold tol x (interior scale) |
| E residual | goes to E (rd); R1/R2 chain; E < 0 -> marks, BE, all-low (rad_m1_clean.cpp:917-:930) | positivity marks (clean_pos_limit) |
| Newton exchange with E_rad/e_gas ~ 1e6 | T from `CleanExchangeTemperature` (rad_m1_clean.hpp:305): monotone root of e(d,T) + A T^4 = e_old + h c kE max(E,0) + EP, bracketed by T(d, rhs). So T <= what the absorbed radiation can give, and with the T^4 term ~ T_rad. Convergence is per cell relative, \|dE\|/E and \|dT\|/T (rad_m1_clean.cpp:1376), not global, so floor cells converge on their own scale | same root (the gas_exact write-back uses it too) |
| positivity fixes | none that add energy: eint < 0 or below the EOS tfloor after the write-back -> cln_neint -> the stage is redone by BE, or the step by dt_retry (rad_m1_clean.cpp:982-:986, :2400) | same |
| KE caps | none in the clean path. implicit_clip_kecap / implicit_am_cs_kecap are `implicit_*`, listed as no effect (rad_m1_clean.cpp:185). time2_inc_kecap (rad_m1_time2.cpp:1385) acts on the stage increment only if named. In staggered that increment is c1 K1 + c2 K2 of physical forces, so it should not fire | same |
| conservation cost | gas + radiation momentum and energy exact through the faces (the face equation is the F equation), L_z to round-off on cs (`CleanStagMomentum` :114) | exact (conservative) |

Conclusion for the USER's question: stag2 has no residual/rho path in the gas momentum or energy. It is the "per-cell exact
local exchange" fix proposed for production, applied everywhere. No e_gas/E_rad switch is needed.

## 3. CPU test (Orion hold 230775, 12 arms concurrently, 1 core each, ~6 min wall)

Deck: m1_test = rhse (rad_m1_beam pgen), 1-D-like sp wedge 96 x 4, r 1-1.115, grey kappa_P = kappa_E = kappa_F = 7.5e4,
ideal gas, gm 1, Eddington ratio 0.5, base tau 1e3, c = 40 (c/c_s ~ 400), vet_col + vet_gd, hesdirk2, 300 cycles, no
seed. The IC (isothermal-tau Eddington E with T from E) is not the scheme's equilibrium, so every arm starts with the same
adjustment, which gives max \|v\|/c_s ~ 2 in all arms. Floor gas = rho < 1e-5, from 1e-5 down to 4e-10 (dfloor 1e-12).
E_rad/e_gas 5e7 there, E/(rho c^2) up to ~20, T_IC 1.0-1.8e-3.
- "_lt" arms loosen the linear tolerance to 1e-5 (clean_tol, or implicit_lin_tol for production) to drive the residual
  route.
- Production arms use scheme = default of the same binary: the production rad_m1 path at the clean-cc base, not
  0940c7d1_t9fix2. No prod-1118 build has the rhse pgen.
- Deck, arms and parser: /orion/ptmp/jinma/clean_1115/m35floor (decks/base.athinput, arms_hold.txt, ana.py, orun/N/).

| arm | dt min / last | DT_RETRY | floor max\|v\| (tab every 0.002) | floor max T | floor min rho | max \|v\|/c_s (hst) |
|---|---|---|---|---|---|---|
| stag (tol 1e-10) | 2.63e-3 / 2.63e-3 | 0 | 7.1e-2 | 1.63e-3 | 4.3e-10 | 2.21 |
| stag_lt (tol 1e-5) | 2.63e-3 / 2.63e-3 | 0 | 7.1e-2 | 1.63e-3 | 4.3e-10 | 2.21 |
| cell gas_exact on | 2.12e-3 / 2.12e-3 | 0 | 1.07e-1 | 1.64e-3 | 4.3e-10 | 3.28 |
| cell gas_exact off | 2.12e-3 / 2.12e-3 | 0 | 1.07e-1 | 1.64e-3 | 4.3e-10 | 3.28 |
| cell yres | same as gas_exact on | 0 | 1.07e-1 | 1.64e-3 | 4.3e-10 | 3.28 |
| cell gas_exact on, lt | 2.12e-3 / 2.12e-3 | 0 | 1.07e-1 | 1.64e-3 | 4.3e-10 | 3.29 |
| **cell gas_exact off, lt** | **2.2e-4** / 3.1e-3 | **16** (+7 auto-BE) | 5.8e-2 | **1.30e-2 (x8)** | 4.0e-10 | **6.00** |
| cell yres, lt | 2.12e-3 / 2.12e-3 | 0 | 1.07e-1 | 1.64e-3 | 4.3e-10 | 3.28 |
| prod matched | 2.49e-3 / 2.49e-3 | 0 | 9.8e-2 | 1.59e-3 | 4.3e-10 | 2.16 |
| prod matched, lt | 2.49e-3 / 2.49e-3 | 0 | 9.8e-2 | 1.59e-3 | 4.3e-10 | 2.16 |
| **prod conservative2** | 1.28e-3 / 1.28e-3 | 3 (NON-CONVERGED) | **1.8e-1** | 1.59e-3 | **1.0e-12 (dfloor)** | **4.82** |
| **prod conservative2, lt** | **2.2e-6 / 8.3e-6** (x300 drop from cycle 244) | 12 (NON-CONVERGED) | **1.9e-1**, at r = 1.112, rho 2e-11 | 1.59e-3 | 1.5e-11 | **5.14** |

There were no CLEAN WARNING lines and no stage with E < 0 or with R2 violations. The production runs printed no
INCKECAP lines because the key was not named. No restart dumps were written.

Reading:
- **stag2 is insensitive to the solve residual.** Its arms are identical at 1e-10 and 1e-5 to the printed precision.
  Among the clean arms it has the lowest floor \|v\| and the largest dt.
- **conservative2 reproduces the M35 signature in miniature:**
  - the floor gas is accelerated to about 2x the matched speed and evacuated to dfloor;
  - dt is about half of matched at the default tolerance, and collapses by x300 (to 2e-6) with loose tolerance.
- **matched is clean,** as viper found on M35.
- **Cell layout:**
  - the energy route (gas_exact off plus residual) heats the floor x8 and cuts dt by x15;
  - gas_exact on removes it.
  - The momentum route (yres) made no difference here: even at tol 1e-5 the y residual in the floor was below the
    physical force. It stays a by-construction exposure, smaller than conservative2's because the cell y has no
    separate face Z_f mismatch.

Limits of the test:
- ideal gas, constant grey opacity, code units, 1-D-like wedge, no MESA/TOPS, no time2_inc_kecap, no 3-D corner;
- the production arms run the production path of the clean branch's base, not the M35 binary.
The conclusion rests on the code (section 1). The test confirms the direction and the ordering: stag2 = matched-like <
cell gas_exact < cell no-exact < conservative2.

## 4. Recommendation
- stag2 (and cell with clean_gas_exact = true, the default) would not have the M35 floor pump.
- For the cell layout, consider making clean_yres = true the default before production on floor-heavy decks. It is one
  key; the cost was not measured here.
- Production M35 stays on matched + lzx off as decided.
