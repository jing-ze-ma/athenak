CHECKPOINT 1: cloned athena master e47fabb7dd5efb22e1e2c7ff193148d723e3e952
# NOTE athpp_audit_1015 (worker): Athena++ radiation audit re Linhao BSG collapse
- Athena++ clone athena/ (origin/main == "master"): e47fabb7dd5e (2026-10-14). Radiation last touched 66adc6bd/42346cc2 (2026-05-14, pole/psi fixes).
CHECKPOINT 2: taucell = <radiation>/taucell default 5 (rad_integrators.cpp:63); multiplies cell tau in HLL tau-factor f=1/(taucell(tauL+tauR)) (implicit_first_order.cpp:75-81, GetTaufactor rad_integrators.cpp:615) and adv factor.
CHECKPOINT 3: analytic: thick-limit HLL numerical diffusion carries eps = Sum w|mu|/(4 taucell Sum w mu^2) ~ 3/(8 taucell) of L (7.5% @5, 9.4% @4) WITHOUT radiation force (force from cell-centred F = phys part only). Candidate: radiative acceleration deficit Gamma/(1+eps). Testing with 1D slab.
CHECKPOINT 4: slab test (pgen athena/src/pgen/thickslab.cpp, runs slab/): thick-region face flux/cell F = 1.0759 (eps=3/(8 taucell) confirmed). NEW: sharp tau jump -> per-angle RHS<0 (cells 149/150) -> TINY floor in CalSourceTerms creates radiation (L_out 11x L_in); with source_flag=1 the created energy is billed to the GAS (tot-E falls at ~ excess L). Scan slab/scan (taucell 4/5 x hcells 0..4) running.
CHECKPOINT 5: scan done: taucell5 face/cell 1.0761, taucell4 1.0951 (eps confirmed); spurious L_out/L_in = 10.5/5.0/2.7/1.0/1.0 for tau e-fold per 0,0.5,1,2,4 cells (c dt/dx 310). Force check: rad momentum source ∝ F_cell (not face flux) -> radiative accel deficit 1/(1+eps). bsg_star.cpp does NOT compile on public Athena++ (any branch/history): uses Hydro::rad_source_term -> Linhao runs a locally MODIFIED Athena++. Patched copy (output only) in athena_bsg/.
CHECKPOINT 6: 1D Athena++ run of bsg_star.cpp (his pgen, output patched) in bsg1d/: units rho 1e-9, T 1e5, l Rsun, mu .62, GM 20 Msun, 20-80 Rsun nx1 384, nmu 4, cfl .3; arms tc4 / tc40 / tc400 (nlimit 100, err 1e-6) and tc4n (nlimit 400, err 1e-9), 30 d. Note: Jacobi hits nlimit 100 at err ~3e-5 at start.
- tc40/tc400 diverged at cycle 1-2 (Jacobi rel err 0.6 -> NaN): near-central DO flux unstable; removed. New arms tc1/tc10/tc20 (run2.sh).
CHECKPOINT 7: 1D Athena++ (his pgen) results: ALL arms drain energy & collapse at the INNER WALL: cell1 T 1.05e6->2.8e5 K, rho x3, cell Fr1 -0.1..-0.3 code (= -100..-300 L, INTO the bath) within 1-3 d; Egas halves in 16 d (tc1), r(rho=1e-10) 52->44-48 Rsun; tc4 NaN at 5.5 d, tc4n at 3.5 d, tc10/20 crash, tc1 survives 16 d. Floor diagnostic (athena_bsg rad_source.cpp DIAG): RHS<0 floor fires at cell i=1 (r 20.23, tau_cell 1.4e6) with negE ~0.7 E0 per stage. Next: no-floor variant to test causality.
- nf4/nf1 (pre-solve RHS floor removed, athena_nf): wall odd-even/inward-flux pattern STILL develops by 1 d (cell0 T 12.4, F -0.29) -> floor is not the only driver. Testing no-Doppler-boost bath (athena_nb, nb4) + ghost-1 velocity print (gv4).
CHECKPOINT 8 (final): v_ghost1 (HLLC zero-contact) ~ +4.2e4 cm/s early -> bath Doppler boost ~ +26 L. nb4 (no boost) and nf4/nf1 (no RHS floor) STILL show the wall pattern (cell0 hot/cell1 cold, Fr ~ -0.2..-0.4 code = -200..-400 L into the bath) and the energy drain (nf: Egas -8% in 5 d vs tc4 -19%). => the floor is an amplifier (~half the gas drain), not the sole cause; the core is his heat-bath ghost + Athena++ thick-limit (near-central) DO face flux at the wall. tc1 survives to 21 d with Egas 4.1e9 -> 0.9e9 and r(rho=1e-10) 52 -> 44 Rsun, then blows out (22 d). nb4 still running at report time.

# HEAD-TO-HEAD 1D (Athena++ vs AthenaK, same setup)
CHECKPOINT H1: AthenaK binary h2h/bin/athena_he_cpu_d67d2599 (fin51-linhao-bc d67d2599 = fin51-r-fix b16b66db + Linhao IC/BC; md5 a68e623b; CPU OpenMP+MPI gcc13).
- 1-D: AthenaK rad_m1 on a 1-D mesh FATALs (implicit_vimp needs multi-D; 1-D uses a different column-solve path); nx2=nx3=2 FATAL (mesh >= 4);
  USED 384 x 8 x 8 (theta pi/3..2pi/3, phi 0..pi/3, periodic, no perturbation -> stays laterally uniform), 1 block.
- Deck: production Linhao deck bsgLI_full_fin3_bc.athinput (FIN3 keys, hesdirk2, vet_gd, he_bc_linhao) with: nx1 384 uniform 20-80 Rsun
  (stretch off), lateral 8x8, he_ic_linhao_amp 0, tlim 30 d, outputs 0.25 d. Arms (h2h/decks): K1 prod keys cfl 0.3 (dt ~175 s);
  K2 cfl 0.364 (dt ~ Athena++ 212 s); K3 = K1 with implicit_pos_floor/_solve/pos_gas false; K0 = K1 with he_bc_linhao false (our fixed-flux BC).
  Athena++ bsg1d dt 212 s (cfl 0.3, gas sound speed). Batch job h2hK (h2h/run_k.sbatch, general, 4 arms x 18 OMP).
CHECKPOINT H2 (RESULTS, HEAD-TO-HEAD 1D). Figures h2h/h2h_user.png (user), h2h/h2h_timeseries.png, h2h/h2h_profiles.png; numbers h2h/h2h_summary.txt, h2h/h2h.json; script h2h/h2h.py.
- dt: Athena++ 212 s (cfl 0.3, gas cs); AthenaK K1 175 s mean (cfl 0.3, rad_signal_speed), K2 218 s (cfl 0.364).
- Athena++ arms (his bsg_star.cpp; runs h2h/runs/A_*): A_tc4 (his taucell 4, angle_flag 0) and A_nf4 (no RHS clip) die at 5 d;
  A_af1 (angle_flag 1, nzeta 4: spherical angular flux) dies at ~24 d. ALL: wall cell 1 T 1.05e6 -> 2-3e5 K, rho x3-4, wall F
  -100..-400 L (into the bath), dE_tot -3.4e46 (5 d, af0) / -7.1e46 (10 d, af1) erg full sphere (~100-300 x L t);
  A_af1: r(rho=1e-10) 54 -> 49 Rsun at 10 d, 45 at 17 d (envelope deflates), Fe-bump <T> 1.77e5 -> 1.0e5 K by 16 d, then blow-out.
- AthenaK arms K1 (prod FIN5.1 keys, his BC), K2 (Athena++ dt), K3 (no pos floors: BITWISE = K1 -> floors never fire), K0 (our flux BC):
  all 30 d, 0 FATAL. Wall T1 1.049e6 -> 1.057e6 K (+0.8 %), rho1 -1.8 %; bath L_in mean 2.5 L (K1) / 2.8 L (K2) INTO the star;
  dE_tot +1.6e45 (K1) vs own int(L_in - L_top) +2.3e45 (30 % hst-diagnostic gap); K0 dE_tot +2e43 (int 8.8e43).
  r(rho=1e-10) oscillates 50.5-54 Rsun first 8 d, settles 51.5-52.5. Fe bump: <T> 1.75e5 steady, max|v| <= 3.7 km/s, Gamma 0.7-0.96.
- Terms: (1) WALL: same BC, opposite sign: Athena++ bath DRAINS 100-400 L, AthenaK bath FEEDS +1.5-2 L; this is the collapse driver.
  (2) FLOORS: Athena++ RHS clip amplifies (tc4 vs nf4: wall T1 3.2e5 vs 5.3e5 at 5 d) but nf4 still collapses; AthenaK floors inert.
  (3) angle_flag 0 (Athena++ default) in spherical 1D: global inward acceleration ~50-70 cm/s^2 (8-10 % of g) at t < 0.1 d,
  INDEPENDENT of taucell (1,2,4,5 identical) -> missing spherical angular-redistribution term; angle_flag 1 removes it (<12 cm/s^2,
  ~AthenaK). (4) force deficit eps (taucell): not visible at 0.1 d (IC F set exactly; would build up on diffusion time); subdominant to (1).

# EMULATION (Athena++ ingredients inside AthenaK, h2h 384x8x8 harness)
CHECKPOINT E1 (interim): branch athpp-emul-1015 (worktree /raven/u/jinma/ATHENAK/wt-athppemul) from d67d2599; commit 1c210bb0
(CPU bin h2h/bin/athena_he_cpu_1c210bb0, md5 036c0def): rad_m1 keys emul_wall_do (R1), emul_clip (R2, pgen user source), emul_af0 (R3),
emul_hll_diff (R4), emul_taucell, emul_absmu; each read only when named. Decks h2h/decks/{E0,R5a,R5b,R5c,R1,R1t1,R2,R1R2,R3,R1R3,R4}.
- Smokes (900 s): R3, R4 run (0 FATAL). R5a FATAL at setup: implicit_flux = central is rejected with implicit_flux_faces = all
  (needs berthon|blend) -> R5a needs implicit_flux_faces changed too (to do). R2 v1 (unclipped P1 intensities) collapsed hydro dt
  175 -> 2.5 s in 4 cycles (P1 I < 0 in the streaming atmosphere = spurious clips) -> fixed: P1 clipped at 0 (uncommitted, rebuild).
- R1: runs, but dt shrinks 175 -> 153 -> 97 -> 28 s by cycle 5 and Picard NON-CONVERGED at dt ~0.03 s by cycle 26 (dt_retry halving);
  ~20 s wall per cycle. Under investigation (wall face explicit 0.5(F_g+F_0^n) with the boosted ghost F; cause of the dt drop not yet known).
- Bitwise key-off check (E0 = K1 deck on 1c210bb0 vs d67d2599 K1, 30 d): not yet run.
CHECKPOINT E2: 9a68a5e3 = R1 acts on the ENERGY row only (v1 1c210bb0 also put the DO flux into the face's M1 flux f0 -> the half
ghost boost (4/3 v_g E_g ~ 23 F) entered the wall momentum deposit: v0 +5.8e5 -> -3.4e6 -> +2.4e7 -> -1e8 cm/s (vceil) in 4 steps);
R1 smoke now dt 175/150 s normal. R2 v2 (clipped P1) still created E in the streaming atmosphere (hydro dt -> 1e-3 s) ->
d6f557ba restricts R2 to tau_cell >= 1 cells (Athena++'s actual clips were at the wall, tau 1.4e6). R5a needs, besides
implicit_flux = central: implicit_flux_faces x1, implicit_thin_switch none, implicit_flux_beam closure, implicit_blend tau
(setup FATALs otherwise); smoke: hydro dt 175 -> 31 s at cycle 2 (dt COLLAPSE).
- FINAL binary h2h/bin/athena_he_cpu_d6f557ba (md5 309a6253). 30-d jobs: emulE0 31207682 (E0 R1 R1t1 R4), emulR2 31207683 (R2 R1R2 R3 R1R3),
  emulR5a 31207684 (R5a R5b R5c). R5c mapping (Athena++ nlimit 100 Jacobi / L1-change 1e-6): implicit_tol 1e-6 (Picard |dE|/E),
  implicit_lin_tol 1e-6, implicit_lin_maxit 20, implicit_maxit 5.
CHECKPOINT E3 (d6f557ba, h2h/emul_summary.txt, emul_timeseries.png, emul_user.png, emul_profiles.png):
- BITWISE key-off: E0 (K1 deck, d6f557ba) vs K1 (d67d2599), 30 d: user.hst + hydro.hst IDENTICAL; final rst bsg3d.00006.rst cmp IDENTICAL.
- 30-d arms, all stable, NONE drains/collapses (K1: wall T1 1.057e6 K, L_in +2.5 L, dE_tot +1.6e45, r(rho=1e-10) 52.1 Rsun at 30 d):
  R1 (DO wall face, taucell 4) / R1t1 (taucell 1): IDENTICAL to each other to 4 digits; L_in +3.2 L (bath FEEDS), wall T1 1.061e6,
  dE_tot +2.5e45, rsurf 52.1. R4 (HLL eps diffusion, energy only): = K1 (+0.01 %), L_top 1.02 L. R5b (GM x1.0078): L_in 2.9 L,
  wall cell F -16 L, rsurf 52.1, dE +2.2e45. R5c (loose solve, at 21.7 d, NC every step, dt 17-175 s): = K1 within ~10 %.
  R3 / R1R3 (angle_flag-0 curvature): cell-0 F -81 L (cell-centred, the extra 2P/r term), L_in +4.4 L, dE_tot +4.6e45, envelope
  EXPANDS to rsurf 68 Rsun at 5 d then back to 50.9 Rsun at 30 d; wall T1 +1.8 %.
- Failed arms: R2 (per-angle clip, tau >= 1 cells) died 0.73 d (dt_retry 9x, dt 5.6e-5 s): the P1 clip fires in the outer envelope
  (Fe-bump Gamma -> 35, max|v| 1000 km/s at the top, blow-OUT, wall T unchanged) -> emulation not faithful (Athena++'s clips were at the
  wall only); R1R2 same behaviour, cancelled at 1.6 d (dt 6-8 s, 1594 retries). R5a (central flux) died at t = 486 s: hydro dt COLLAPSE
  175 -> 8e-3 s within 8 cycles (dt_retry FATAL) = central M1 flux without the HLL/blend part is unstable in this setup (a result).
- note: R3/R1R3 r(rho=1e-10) oscillates 44-71 Rsun (minima 44-47 at ~8, 21, 28 d) while E_tot RISES; wall T1 unchanged.
CHECKPOINT E4 (FINAL): figures h2h/emul_user.png (user), emul_timeseries.png, emul_profiles.png; numbers h2h/emul_summary.txt, emul.json.
Branch athpp-emul-1015 (LOCAL, not pushed) @ d6f557ba (1c210bb0 keys -> 9a68a5e3 R1 energy-only + R2 P1 clip -> d6f557ba R2 thick only);
CPU bin h2h/bin/athena_he_cpu_d6f557ba (md5 309a6253). Key-off BITWISE vs d67d2599 (E0 vs K1, 30 d hst + final rst).
VERDICT: none of R1 (taucell 4, 1), R3, R1R3, R4, R5b, R5c drains the wall or loses energy; all keep the bath feeding (+2.4..4.4 L) and
E_tot rising (+1.2..4.6e45). R3 (angle_flag-0 term) is the only one that changes the envelope strongly (oscillation 44-71 Rsun, no secular
deflation). R2/R1R2 (clip) blow OUT the envelope (not faithful), R5a (central) is numerically unstable at t < 0.01 d.
=> The Athena++ collapse is not reproduced by any single term transplanted into the moment scheme; it needs the per-angle DO
odd-even wall mode itself (Athena++ nf4/af1 still collapse without the clip / with angle_flag 1).
R5c left running in job 31207684 (at 21.75 d, ~9 h to 30 d at ~1 d/h); summary uses it to 21.7 d.
- R5c finished 30 d (exit 0); summary/figures regenerated with full R5c.

# BC SWAP (inside Athena++: his bsg_star.cpp 1-D, only the inner radiation BC swapped)
CHECKPOINT S1: tree athena_bc (= athena_bsg + pgen option problem/bc_swap, read only when named; 0 = his BC); athena_bcnf = athena_bc +
the nf patch (pre-solve RHS clip removed). Mapping (code units: E = sum w I, F = sum w mu I, sum w mu^2 = 1/3, P1 ghost I = E + 3 mu F):
- Ba (bc_swap 1, OUR fixed-flux BC): ghost I = E_0 + 3 mu F_in at every ghost layer, E_0 = cell-is E (continuous), F_in(r_g) =
  Fr(ghost 1) r_g1^2/r_g^2 (= L/4 pi r^2), no Doppler boost. Hydro ghosts = his (bd_data + HLLC zero-contact v).
- Bb (bc_swap 2, diffusion-limit bath): ghost I = T_g^4 + 3 mu F_b, F_b = -(E_0 - T_g^4)/(3 chi_f (r_is - r_g)), chi_f = mean(chi_g, chi_is),
  T_g = his bd_data T (profile), no boost.  Bc = Bb on athena_bcnf.
- Default check: Bchk (A_tc4 deck, athena_bc, bc_swap absent) vs A_tc4: all tab dumps 0..0.9 d and the 100 hst rows to 1 d IDENTICAL; the
  1-d dump differs only because tlim truncates the last step (t 4.547922 vs 4.555754 code). => default = his BC bit-for-bit.
- Arms Ba0/Ba1/Bb0/Bb1/Bc0/Bc1 (0/1 = angle_flag 0 / 1, decks = A_tc4 / A_af1 + bc_swap), batch job bcswap 31212475 (out dirs runs/B*/outb).
CHECKPOINT S2 (RESULTS; h2h/bcswap_user.png, bcswap_timeseries.png, bcswap_profiles.png, bcswap_summary.txt, bcswap.json):
Added two diagnostic arms beyond the brief: Bd (bc_swap 3) = CLOSED fixed-flux radiation wall (ghost I(mu) = I_is(-mu) + 6 mu F_in: mirror
+ imposed flux, DO face flux = F_in exactly, nothing advected out), Be (bc_swap 4) = Bd + REFLECTING hydro ghosts (no zero-contact loop,
no bd_data). Both af0 and af1 run for every arm (cheap: ~12 min for 30 d on 1 core).
| arm | t_end d | T1 5d K | dEtot 5d erg | rsurf min 5-20 d |   (K1 AthenaK: T1 1.05e6, dE +8e44, rsurf 50.5)
| B0 A_tc4/A_af1 | 5.0 / 23.8 | 3.2e5 / 3.8e5 | -3.8e46 / -3.2e46 | 54.8 / 44.9 |
| Ba0/Ba1 our flux BC | 7.2 / 30 | 3.6e5 / 4.3e5 | -5.3e46 / -4.9e46 | 44.0 / 45.1 |
| Bb0/Bb1 diffusion bath | 5.0 / 30 | 4.1e5 / 4.1e5 | -5.2e46 / -4.0e46 | 57.0 / 43.4 |
| Bc0/Bc1 Bb no clip | 5.0 / 30 | 6.3e5 / 6.3e5 | -4.2e46 / -3.3e46 | 55.7 / 45.5 |
| Bd0/Bd1 closed flux wall | 5.0 / 30 | 5.4e5 / 3.8e5 | -4.7e46 / -4.7e46 | 55.9 / 44.9 |
| Be0/Be1 closed + reflect hydro | 5.0 / 30 | 5.7e5 / 4.9e5 | -4.2e46 / -4.5e46 | 55.5 / 46.0 |
(af0 arms end at 5-7 d = NaN crash; af1 arms reach 30 d but blow out the top after ~18 d like A_af1.)
VERDICT: swapping ONLY the inner boundary does NOT stop the collapse in Athena++ -- not with our fixed-flux BC, a diffusion bath, a closed
fixed-flux radiation wall, nor a fully reflecting wall. => his inner BC is NOT the cause in his own code. The clip only softens it (Bc vs Bb:
T1 6.3e5 vs 4.1e5 K at 5 d). Energy budget (Bd1, A_af1 at 5 d): ALL of dE_tot (-4.7e46 / -3.2e46) is lost inside r < 30 Rsun (cumulative
dE(<r) flat beyond 30 Rsun), with 0.02-0.03 Msun moved inward into the first ~1 Rsun; the wall face carries only F_in in Bd/Be, so the
energy disappears in the interior of the Athena++ solution near the thick base (non-conservation in the radiation-hydro update there),
not through the boundary. Cell-centred lab flux at 0.1 d is ~-50 L throughout r < 22 Rsun in every arm (gas infall ~1 km/s advecting
(4/3) v E), i.e. the "wall flux" metric is advective, not a bath leak.
Next diagnostic (not done): per-term energy ledger in Athena++ at the base (gas eint<0 / t_ceiling resets in AddSourceTerms, the pre-solve
clip, the point-mass energy source vs face mass flux, rad-source billing), to name the sink.

# LEDGER (per-term energy ledger inside Athena++, 1-D his bsg_star.cpp)
CHECKPOINT G1: tree athena_led (= athena_bc + src/led.hpp ledger, key problem/ledger = 1; diagnostics only). Accumulated per cell, LAST
stage only (vl2 discards stage-1 writes), x cell volume: hydflux (AddFluxDivergence IEN), gravmass (Phi_i x mass-flux divergence),
pointmass (IEN source), radT_bill (transport the gas billing uses: (E_new - E_ini) + rad_source0/prat, i.e. divflx_ of the LAST Jacobi
sweep = neighbours' previous iterate), radT_cons (same with every face rebuilt from the FINAL intensities: telescopes = boundary flux),
clip (gross per-angle RHS deficit; NOT energy created -- the solve divides it by (1 + cc + c dt chi)), reset (applied IEN - billed in
AddSourceTerms), floorE/floorM (EOS floors), radD, gasbill. E_tot = IEN + prat Er + rho Phi.
- Bitwise: ledger off (A_tc4 deck, 0.9 d) = A_tc4 (30 dumps identical); ledger ON (L_tc4) = A_tc4 to 5 d (150 dumps identical).
- Closure: measured dE_tot(region) = sum of terms to <= 1e-11 relative (all regions, all times).
RESULTS (full sphere, 1e45 erg, change since t = 0):
| run | region | t | dE_tot | non-telescoping rad transport | consistent rad transport (boundary) | hydro | gravity | resets+floors |
| L_Bd1 (closed wall, af1) | all | 5 d | -47.4 | -47.8 | +0.3 | +0.0 | -0.0 | +0.1 |
| L_Bd1 | all | 10 d | -98.1 | -98.6 | +0.3 | +0.0 | -0.2 | +0.3 |
| L_Bd1 | r<30 | 5 d | -52.7 | -48.2 | +1.9 | +0.7 | -7.2 | 0 |
| L_Bd1 | cells 0-10 | 5 d | -129.0 | -49.2 | +18.1 | +18.1 | -116.0 (mass x Phi redistribution, globally 0) | 0 |
| L_tc4 (his BC, af0) | all | 5 d | -37.7 | -12.0 | -25.8 (bath drain through the inner face) | 0 | 0 | +0.1 |
| L_tc4 | r<30 | 5 d | -44.5 | -12.3 | -23.8 | +0.9 | -9.3 | 0 |
Where: Bd1 at 5 d the non-telescoping loss is cells 0/1/2 = -22.3/-12.2/-8.0 (89 % in cells 0-2; 103 % in cells 0-10); tc4 cells 0-10 -11.4.
When: Bd1 non-telescoping starts at ~0.05 d and grows ~linearly (-0.1 at 0.17 d, -1.0 at 1 d, -7.8 at 2 d, -48 at 5 d, -99 at 10 d).
Mechanism: the implicit Jacobi iteration stops at nlimit = 100 with relative change ~1e-5 > error_limit 1e-6 (17/18 sampled solves);
GetHydroSourceTerms bills the gas with divflx_ of the last sweep (neighbour intensities from the previous iterate), so the face
fluxes do not telescope and the mismatch at the thick wall (odd-even mode, slow Jacobi convergence) is destroyed energy.
TEST (converged iterations, C_Bd1: error_limit 1e-11, nlimit 5000; converges in 1000-1700 sweeps): non-telescoping = 0.0000 and dE_tot
stays at -0.025 (initial transient) from 0.02 to 0.1 d, while L_Bd1 loses -0.0925 by 0.097 d and the rate keeps growing. C_Bd1 is ~1 d
per 13 h of CPU, so 30 d is not affordable here; only the early-time contrast is shown.
FIX VARIANT tried (problem/ledger_fix = 1: rebuild divflx_ from the final intensities before the billing = conservative billing of an
UNconverged iterate): FAILS -- hydro dt -> 2.5e-10 s at t = 200 s (the per-cell transport mismatch of an unconverged iterate is huge;
billing it to the gas destroys the base). => the actual cure is converging the iteration, not re-billing.
Figure: h2h/ledger_user.png (left: whole-domain energy by term vs t for L_Bd1 / L_tc4; middle/right: cumulative by term vs r at 5 d).
- C_Bd1 update at 0.53 d: dE_tot -0.030, non-telescoping -0.0000; L_Bd1 at 0.53 d: dE_tot -0.551, non-telescoping -0.491 (18x larger loss). C_Bd1 continues in the interactive srun (1.5 h limit).
- C_Bd1 FINISHED 2 d (318 s CPU after the slow start): dE_tot +0.43e45, non-telescoping 0.0000; L_Bd1 at 1.99 d: dE_tot -7.61e45, non-telescoping -7.62e45. => with a converged iteration the Athena++ energy loss disappears (to 2 d).

# THINCLOS (BSG closure-switch test; work dir /raven/ptmp/jinma/thinclos_1016)
CHECKPOINT T1: branch thinclos-1016 (worktree /raven/u/jinma/ATHENAK/wt-thinclos) from b16b66db; 02ddf595 (key) + 34b1b5f8 (include fix).
Key <rad_m1>/thin_closure_switch = m1 | eddington (| off), thin_closure_tau (default 1), thin_closure_start (default cycle 1), read only
when named. RadiationM1::ThinClosureOverride (rad_m1_vetcol.cpp, called at the end of every VetColBuild, i.e. after vet_col + vet_gd
(VetLatBuild) of every build): per column (one MeshBlock along x1 required), tau_R from the top summed in the block (kappa_T = M1_OP_T,
half own cell), every cell with tau < thin_closure_tau gets tau_ten slots 0..3 = (chi, n) of M1 (chi = M1Chi(|F|/cE), n = F/|F|, from
u0) or Eddington (1/3, n = F/|F|), lateral slots LAT0..LAT0+5 zeroed. NOT changed: the formal-solution J/H/half-range moments
(implicit_flux blend/halfrange beam, vet_source), vet_col surface q (deck: Marshak q fixed).
Builds: bin/athena_he_{cpu,a100}_34b1b5f8 (md5 f878851b / a4bd0fc3). Jobs: gate 31224483 (base/m1/eddington, 5 d), E1 key-off
31224484 (gpudev, 46 cycles vs overlap_1015 XB), CPU key-off 31224485 (BSG 6 cycles vs 0da548c5 cpu = b16b66db object code).
CHECKPOINT T2 (RESULTS; thinclos_1016/ana/thinclos.png, thinclos.json):
- Key-off BITWISE: A100 E1 production deck 46 cycles vs overlap_1015 XB (b16b66db): rst BITWISE + 94 hst rows IDENTICAL; A100 BSG gate
  base arm (5 d) vs icstd bsg_P10 (b16b66db, same deck/keys): first 1440 rows of hydro.hst and user.hst IDENTICAL; CPU BSG gate deck
  6 cycles vs 0da548c5 cpu (= b16b66db object code): hydro.hst + user.hst IDENTICAL. (CPU E1 not run.)
- Switch active: 26368 cells/rank switched (tau_R < 1: ~103 of 256 radial cells) from cycle 1; arms differ early (0.5-1.5 d).
- Gate (a), 0.99-1.1 R_ph (R_ph 51.57 Rsun), max|<v_r>| [km/s] / rms / max|rho/rho0 - 1|:
| t [d] | base (vet_gd) | m1 | eddington |
| 0.5 | 0.47/0.18/0.43 | 0.30/0.15/0.43 | 0.66/0.24/0.43 |
| 1.0 | 2.46/0.73/0.56 | 1.49/0.49/0.56 | 3.17/0.95/0.62 |
| 1.5 | 4.97/1.48/0.62 | 3.54/1.30/0.62 | 7.02/1.93/0.62 |
| 2.0 | 3.21/1.46/5.0 | 3.23/1.50/4.4 | 3.35/1.49/5.9 |
| 3.0 | 20.2/9.5/206 | 20.3/9.5/200 | 20.1/9.5/210 |
| 4.0 | 27.4/10.0/105 | 27.1/9.9/104 | 27.4/10.0/103 |
| 5.0 | 27.4/8.1/1.9 | 27.5/8.2/1.9 | 27.1/8.0/1.9 |
| max 0-5 d | 54.3 (4.8 d) / 25.6 / 222 | 54.6 / 25.5 / 219 | 54.1 / 25.3 / 226 |
  FFT periods (shell-mean v_r, 0.5-5 d), all three arms identical: 2.24 d at 50.0-51.1 Rsun, 1.49 d at 51.6-52.8 Rsun, 1.12 d at 53.6 Rsun.
VERDICT: the oscillation, its periods (1.49 d at R_ph = the photospheric acoustic cut-off of RESULTS 7b, 1.35 d) and its saturation
(54 km/s at 4.8 d) are the same to a few % with the M1 or the Eddington closure in the tau_R < 1 layer: it is NOT driven by the
vet_gd thin-layer closure; it is a hydrodynamic atmospheric resonance (physical for this setup). Caveat: the switch replaces the
closure tensor only; the formal-solution half-range moments used by implicit_flux = blend/halfrange and vet_source stay VET.
- Branch pushed: origin thinclos-1016 @ 34b1b5f8. No production touched.

# BSG ATMOS (damped-start release test + old deck-P comparison; /raven/ptmp/jinma/thinclos_1016)
CHECKPOINT A1: thinclos-1016 @ 88557a88 (pushed) adds DEBUG problem/he_damp_rate, he_damp_r0, he_damp_r1, he_damp_tend (read only when
named): v_r *= exp(-rate dt) in r0 < r < r1 while t < tend, KE removed from the total energy (as the top sponge). No existing key could
do a time-limited band sponge (he_sponge_* is the top sponge r > r0, always on). Key off BITWISE: base7 (88557a88, no key) = base
(34b1b5f8) on the first 1440 rows of hydro.hst and user.hst (A100). bin/athena_he_a100_88557a88 md5 a9c4aaf4. Job 31226537 (1 A100).
Runs (gate (a) deck/IC, P grid, seed 0, 7 d): base7 (undamped), r15 / r30 = rate 1e-3 /s in 0.98-1.2 R_ph (50.54-61.88 Rsun) until 1.5 / 3 d.
Band 0.99-1.1 R_ph, max|<v_r>| km/s / rms / max|rho/rho0-1|:
| t [d] | base7 | r15 (release 1.5 d) | r30 (release 3 d) |
| 1.0 | 2.46/0.73/0.56 | 0.08/0.05/0.56 | 0.08/0.05/0.56 |
| 2.0 | 3.21/1.46/5.0 | 4.86/2.09/0.66 | 0.17/0.06/0.66 |
| 3.0 | 20.2/9.5/206 | 21.7/10.7/218 | 0.10/0.05/0.72 |
| 4.0 | 27.4/10.0/105 | 26.1/10.1/154 | 6.64/2.43/7.2 |
| 5.0 | 27.4/8.1/1.9 | 28.6/8.6/2.0 | 6.51/2.62/10.2 |
| 6.0 | 17.1/4.4/1.0 | 15.6/4.1/1.1 | 4.11/1.50/4.9 |
| 7.0 | 8.3/4.0/16 | 8.7/3.9/15 | 10.3/4.0/28 |
| max | 54.3 at 4.8 d | 58.5 at 4.85 d | 24.2 at 5.35 d |
FFT (shell-mean v_r, 0.5-7 d): base7 and r15 identical, 1.29 d at 51.6-52.8 Rsun, 2.15 d at 50-51 Rsun; r30 (only 4 d after release)
1.07 d at 52.3-52.8 Rsun, 2.15 d at 50 Rsun.
=> SELF-EXCITED: from 0.08-0.1 km/s at release the band regrows by ~x70 within ~1 d (e-folding ~0.25 d) to the same 5-50 km/s state:
released at 1.5 d it reaches the undamped saturation (58 vs 54 km/s, 4.85 vs 4.8 d); released at 3 d it is at 24 km/s / rms 4 km/s
(= base7's late level) by 5.4-7 d. Not a decaying startup transient: a linear instability of this atmosphere (growth ~0.25 d e-fold,
period ~1.1-1.5 d at the photospheric acoustic cut-off), saturating nonlinearly in shocks. Physical vs numerical: the THINCLOS test
excluded the vet_gd thin closure; resolution/grid was not varied here.
OLD deck-P production (bsgprod_1012/prod/P, 43 dumps, 0.5-d cadence, R_ph 53.03 Rsun of its IC, same band, ana/oldP.txt):
max|<v_r>| 41 (0.5 d), 83 (1 d), 116 (1.5 d), 49 (2 d), 57 (3 d), 42 (4 d), 72 (5 d), 123 (6 d), 143 (7 d), 18 (8 d), 38 (10 d),
80 (15 d), 23 (20 d) km/s; max 165 km/s (9.5 d); rms up to 125 km/s; rho swings to x2.4e4. Periods: not resolvable (0.5-d dumps,
Nyquist 1 d; FFT gives 3.5-10 d, i.e. the slow infall/refill cycles).
=> the old setup is NOT quiet: from t = 0 its atmosphere falls in (the known -83 km/s at 1 d, -148 km/s and rho x9e4 by 6.5 d) with
2-5x larger velocities than the new setup's resonance -- a different, large-amplitude non-equilibrium (infall/refill) motion of its
IC atmosphere (radiative-only extension + isothermal atmosphere with F_r = 0 above r_ph, total-pressure scale height, ~1 cell per H at
R_ph on dr 0.156 Rsun, no sponge), which would mask a 1.3-1.5 d resonance; with 0.5-d dumps one cannot tell whether the resonance is
also present underneath. Linhao-IC productions: Raven prod/LI and Viper prod/LI never ran (no dumps); Caltech's is remote: skipped.
Figure: thinclos_1016/ana/bsg_atmos.png (left: band max|<v_r>| base7/r15/r30; middle: rho swing; right: old P vs new setup).
CHECKPOINT A2 (all other old 3-D BSG runs with dumps; thinclos_1016/atm/*.txt|json, script atm/atm_runs.py, meta.py):
R_ph = innermost r with <E> <= 1.1 x 4 L/(4 pi r^2 c) in dump 0 (checked: P 53.09 = its IC r_ph 53.03). Cells/H = H_P/dr at R_ph,
H_P = ((2/3) e + E/3)/(rho g) (ideal-gas approx; sk67 has the general EOS). ALL old runs dumped every 0.5 d (cFR 0.33 d) => no run
can resolve a ~1.3 d period (Nyquist 1 d; cFR 0.66 d gives 1.9 d at 1.02 R_ph, inconclusive). Band 0.99-1.1 R_ph: max|<v_r>| / median
(t >= 1 d) / median rms [km/s]; infall = fraction of dumps with band-mean v_r < -1 km/s; interior 0.5-0.9 R_ph median |<v_r>|.
| run (where) | code | IC / outer BC / sponge | grid at R_ph | span | atm max/med/med-rms | infall | interior | verdict |
| P deck-P prod (Raven bsgprod_1012/prod/P) | 9972a56d (FIN3) | Ma2026 IC (ic_ma2026_f0, isothermal F=0 atm), copy top, no sponge, ideal | R_ph 53.1, dr .156, H .51 -> 3.3 cells/H | 20.6 d | 165/49/33 | 0.72 | 6.9 | infall-refill |
| repro (Viper merged_repro; Raven repro_sf_4n+cont) | fsc_d6164a6a | Ma2026 IC, noinflow top, no sponge, seed 1e-2 | 53.0, .153, .44 -> 2.9 | 57.4 d | 165/19/12 | 0.59 | 1.7 | infall-refill |
| repro_4n (Raven bsg_raven_1002) | gd_12b78e4f | Ma2026 IC, noinflow, no sponge, seed 0 | 53.0, .153, .44 -> 2.9 | 9 d | 155/72/51 | 0.74 | 4.8 | infall-refill |
| R bsgrec (Raven bsgrec_1011/prod/R) | 0197f020 | Ma2026 IC f0, copy top, no sponge | 53.1, .156, .51 -> 3.3 | 6.5 d | 138/70/51 | 0.71 | 6.1 | infall-refill |
| sk67 (Viper merged_sk67 = Raven sk67_4n + Viper) | fsc_d6164a6a | own SK67 IC (general EOS), hse top, sponge 1e-3 above 53.6 Rsun, seed 1e-2 | 46.0, .133, .70 -> 5.2 | 57.4 d | 60/14/7.6 | 0.43 | 0.55 | resonance-like (quiet interior) |
| arm2r2b_4n (Raven) | e8cf2d64 | own arm2 IC, hse top, sponge 1e-3 above 59.2 Rsun, seed 1e-2 | 50.7, .153, .36 -> 2.4 | 57.4 d | 66/14/7.3 | 0.44 | 0.55 | resonance-like |
| arm2seed2_2n / arm2seed2 (Raven) | bcbaef52 | arm2 IC, same | 50.9, .153, .36 -> 2.4 | 9.5 / 3.5 d | 52/15/8.8; 52/9/4.0 | 0.25; 0.0 | 3.4; 1.2 | resonance-like |
| cFR gd_phys (Raven) | gd_12b78e4f | arm2 IC, same | 50.7, .153, .36 -> 2.4 | 5.9 d (0.33 d dumps) | 55/17/7.9 | 0.26 | 2.6 | resonance-like (FFT 1.9 d at 1.02 R_ph, unresolved) |
| new setup gate (a) base7 (thinclos_1016) | 88557a88 | IC-tool BSG, outflow top + sponge, general EOS | 51.6, .156, .85 -> 5.5 | 7 d (0.05 d dumps) | 54/16/5.7 | - | <0.3 | self-excited resonance 1.1-1.5 d |
(rho swings: in the arm2/sk67 runs the band starts near the floor and fills, so max|rho/rho0 - 1| = 1e5-1e7 there is not an amplitude.)
READING: two families. (1) The Ma2026-IC paper-reproduction runs (P, repro, repro_4n, R): no sponge, isothermal F=0 atmosphere,
~3 cells per H: the atmosphere falls in and refills (band mean < -1 km/s in 59-74 % of dumps, 50-130 km/s infall), 2-5x the
new setup's amplitude, and even the interior moves (median 2-7 km/s). (2) Our own-IC runs with the hse top + sponge (arm2*, cFR,
sk67): a resonance-like atmosphere of the SAME size as the new setup's (max 52-66 km/s at 2.5-3 d, median 14-17, rms 7-9 km/s,
balanced in/out motion) and a QUIET interior (median 0.55 km/s in sk67 and arm2r2b over 57 d) -- this is the "very stable in 3-D"
run the user remembers: the star is stable, the atmosphere above R_ph is not. None of them is quiet in 0.99-1.1 R_ph; the 0.5-d
dump cadence hides whether their atmospheric motion is the same 1.1-1.5 d resonance.
- Caltech TASK prepared (NOT pushed): thinclos_1016/caltech_task/TASK-2026-10-16-caltech-bsgatmos.md + bundle bsgatmos-1016/
  (atm_runs.py, gate_atm.py, bin_convert.py, meta.py, MD5SUMS); intended repo path docs/handover/bsgatmos-1016/.
