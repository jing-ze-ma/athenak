---
name: fecz-box-projects
description: 2026-09-13 B-star and He-star iron-zone (FeCZ) Cartesian box projects staged in bench/bstar_fecz and bench/hestar_fecz (columns integrated with the code's Rosseland table, gas+radiation EOS, MLT alpha 1.5; MESA is NOT on viper). B star 15 Msun (log L 4.28, log g 4.16, Teff 29360 derived): FeCZ T 1.3-2.1e5 K, tau 580-3300, 1.4 Hp thick, base 8.5e9 cm (0.023 R) below tau=2/3, v_mlt 2.4 km/s Mach 0.034, F_conv/F 8e-4, turnover 2.8 h, thermal 1.7 h, Prad/Pgas 0.47, Gamma_e 0.033; box 13 Hp (2 Hp below base to tau 1e-2) x 4 Hp wide; 32/min-Hp with refined top = 2.4e8 cells, 36 GPU-h per turnover, 10 turnovers/24 h on 8 GPUs. He stars (Woosley 2019 He/2) 3/5/8 Msun Gamma_e 0.04/0.09/0.15: FeCZ 1.1-2.1e5 K, 2 Hp thick, 1e7-3.5e7 cm below tau=2/3, F_conv/F 5e-9..3e-7, Mach 3e-3..1.7e-2, local Eddington factor 0.99 (5) and 1.09 (8 Msun = 1D invalid, the 3D case); 8 Msun box 8.2 Hp, refined 9.5e7 cells, 9.6 turnovers/24 h on 8 GPUs; launch 5 Msun first
metadata:
  type: project
---
Opacity: a REAL X=0 helium table was built from the Farag+2024 OPLIB set (Zenodo 15277019; recipe data/stellar_opac/PROVENANCE.md +
tools/stellar_opac/merge_rosseland.py, md5-verified) -> bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt (+Z0.014, X0.2 arms); the
solar-minus-kappa_es approximation would have been 25-46 % too opaque. Code: table EOS with eos_radiation=true is radiation-aware
in c2p, Gamma_1, lhllc and the WB background (no gamma-law assumption); box pgen has NO RT (Newton cooling lid) -> the red-giant
two-stream is being ported to plane-parallel (bench/wt_rgbox, problem/rt_two_stream) for the tau 1e-2 top; static WB has no hook in
the box pgen -> dynamic polytropic WB chosen. Uncertainties: Woosley vs Gotberg surface params (0.3 dex log g moves the FeCZ 5x
deeper), CNO-ash metals vs scaled-solar Z, super-Eddington 8 Msun column. See [[fecz-low-lum-gap]], [[rg-box-sweep]].

B-STAR SMOKE (job 11679292, apudev, 2.4e6 cells, 2 turnovers, no RT, deep Newton thermostat lid, bc_mode 0 both walls): numerics
CLEAN (dt flat 2.40 s, empty event log, zero FOFC, E balanced to 1.3 % of F A t). The FeCZ convects: rms v1 grows on a 0.5-0.7
turnover e-fold to 1.4-1.6 v_MLT (3.9e5 cm/s), v1/v_h 3.6 (vertically elongated; the 100 %-capped transverse radiative operator,
max x_i 1e3 vs cap 0.5, is the suspect), superadiabaticity held to 3 % of the MLT 4.7e-2. LID FAILS: thermostat -> hot rarefied
coherent wind (top T +47 %, rho -46 %, 0.35-0.56 c_s) + z-independent through-flow 5.6e-3 g/cm2/s = 2.3 % of the box mass per 2
turnovers, refilled through the bottom column ghosts. Fix = RT-on top (smoke_rt) + bc_mode 2 (zero mass flux) at both walls;
run >= 4-5 turnovers to saturate. Analysis: bstar_fecz/smoke/analysis/smoke_profiles.py.

HE-STAR (5 Msun) no-RT smoke staged in hestar_fecz/smoke but NOT submitted: tau=10 lies INSIDE its FeCZ (tau 26 -> 1.6, top
2.2e7 cm below the photosphere) and the Newton lid cannot work (box holds 4.8 s of F; dry run dt 0.36 -> 0.053 s in 100 cycles,
Mach 1.6 1-D wind). X=0 is safe in the EOS (eos_metal_ionization is dead code at X=0, set false). The He-star arm must be
RT-on from the start: lid tau 1e-2, Lz 8 Hp, ~2e6 cells at 8/Hp, 15-min job. rad_cap_ang caps 100 % of cells (x_i 5e2-1e3).
HE-STAR RT smoke (11683058): FAILED, dt collapse to 6e-5 s by t=159 s, E x3 then back, mass +1.8 % with closed walls; same handover bug, harder (FeCZ inside the blend).
B-STAR RT smoke RERUN with the fix (smoke_rt2, job 11685738): energy leak GONE (E/E0 1.000 flat), mass +1e-3 then flat, ringing
10x weaker, KE_h > KE_r by 1500 s (3-D motion starting) BUT dt fell 1.49 -> 0.32 s by 1658 s and floors/FOFC/tclamp fire after
cycle 0 (cycle 4129: efloor 1.2e5, fofc 3e4, tclamp 2.8e4) - location under analysis (top tau 1e-2 cells, rho 8e-11, suspected).
smoke_rt2 VERDICT (analysis/rt2_profiles.py): the fix works cleanly - retention -0.2 % of F A t (was +25.6 %), handover deposit
0.32 % (was 12.8 %), diffusion flux 1.02->0.99->0.90 F straight through the blend, T/T_Edd 0.96 over 40 cells above tau 10,
emergent flux 0.89 F oscillating about F, 1-D ringing gone (coherence 0.03-0.4), FeCZ nabla-nabla_ad +4.85e-2 intact (4 % of MLT),
interior drift < 1 % T / 2 % rho. Remaining: a bounded photospheric FOFC/floor transient at 600-1000 s (tau 1.6..4e-3, peak tau 0.1,
zero by 1659 s) and a slow rarefaction of the tau 1e-2 top cell (min rho 8e-11 -> 1.2e-12, 12x above dfloor; c_s ~ rho^-1/2 there
because Prad-dominated) that erodes dt (0.32 s, velocity-limited in the top cell); a growing subsonic 3-D motion in the tau 10-100
region (h-KE doubling per 600 s) of unknown nature. Recommended: dfloor 1e-12 (done in rt2turn/), then the 2-turnover continuation
(~9 h on 2 GPUs, 5 x 2-h restart links, staged in bstar_fecz/rt2turn, NOT submitted); deeper x1max for production.
HE-STAR 5 Msun RT rerun with the fix (smoke_rt2, 11686862) COLLAPSED AGAIN for PHYSICS (smoke_rt2/he_rt2_DIAGNOSIS.txt): the
handover is fine (T/T_Edd 1.001, emergent 1.02 F, E -0.5 % of F A t) but Gamma = kappa F/(c g) = 0.98 at tau 12 with the code's
table (gas pressure carries 1.8 % of gravity), the tau 3-45 Fe-bump layer evacuates into a shocked pile under the closed lid,
voids at 2-3 dfloor set dt via the radiation sound speed (c_s 0.38 c) -> dt 1e-4 s by 85 s; bc_mode 2 walls leak +3.4 % once the
interior leaves the column. dfloor only buys sqrt. DECISION: first He target = 3 Msun (Gamma_max 0.75 at tau 30, gas carries 25 % of
g); 5 Msun needs an OPEN top and is a marginal-stability problem; 8 Msun (Gamma 1.09) must not run in a closed box.
HE 3 Msun RT smoke staged + submitted (smoke_rt_3msun, job 11689071): IC Gamma max 0.747, dt flat 0.1 s over 100 CPU cycles, no
density collapse (min rho 1100 dfloor with dfloor = 1e-3 x top rho = 2.3e-12), 50/134 cells at full diffusion weight. Remaining
defect = bc_mode 2 BOTTOM-WALL MASS LEAK +6.5e-4/s (bottom cell rho +16 %) whenever the interior departs from the column; a
conservative bc_mode 3 (WB-consistent hydrostatic ghosts + red-giant-style exact face-flux cancellation) is being built.
bc_mode 3 DONE (4e0054e7, wt_rgbox, unpushed): ghosts = the wall cell's own (rho,e) walked across the wall by WBAdvance with the
WB closure (WB-consistent), v1 mirrored, + problem/wall_noflux exact cancellation of the face mass and energy flux after RKUpdate
(red-giant style; only dm*h removed where conduction shares the energy channel). 1-D column 1e4 cycles: mass/E drift 1e-15 (mode 2:
2e-4/8e-4; mode 0 with a 10 % T bump: -1e-2); B-star RT input 100 cycles: mass flat 1e-13, no start-up transient. Modes 0/1/2
bitwise. WHY mode 2 leaked: ghost rescaled by the INITIAL column ratio -> not the hydrostatic continuation of the evolved wall cell
-> the WB stencil and the wall Riemann problem are asymmetric -> rectified one-signed mass flux. 3 Msun He smoke with bc_mode 2 ALSO
collapsed (112 s, dt 1e-3, mass +1.2 %, a persistent void; diagnosis pending); rerun staged in smoke_rt_3msun_bc3.
3 Msun He RT smoke (bc_mode 2) collapse DIAGNOSED (smoke_rt_3msun/he3_DIAGNOSIS.txt): NOT the wall (leak +1.2 % saturates at 75 s),
NOT Eddington (Gamma mean peak 0.67, 2 % of cells > 1); a NEW instability: transonic horizontal turbulence born at the rad_tau_lo=10
edge (rms v_h 6.5e6 at tau 6 = 350 x v_MLT, coherent downdraft, filamentary voids at 1 dfloor with 0.45 % fill at tau 6-10, the
tau 14-105 layer drained 4.5x the leak), dt plateau 1e-3 s (dfloor works). The B-star RT box shows the same growth (h-KE doubling
per 600 s, 78 % in tau 10-100) slower. Suspects: the blend edge and the capped transverse radiative operator (rad_cap_ang 0.5
throttles horizontal damping 1e3x, so columns decouple). Tests submitted 09-13 ~13:00 (apudev): bc3 baseline 11689919,
bc3 + rad_tau_lo/hi 100/300 11690266, bc3 + rad_cap_ang 50 11690267 (caveat: cap 50 violates the explicit transverse diffusion
CFL, x_i > 0.5, so a blow-up there is numerical, not a refutation). Real fix candidates: implicit transverse diffusion, or the
two-stream everywhere with a higher lid.
3 Msun test matrix results (09-13 13:00): bc3 baseline (10/100) collapses like bc_mode 2 (dt 1.1e-3 at t=102 s) -> WALL CLEARED as
the cause; cap50 -> dt 4e-16 at cycle 7 (explicit transverse CFL, as predicted; not a test); tau100 (blend 100/300) -> dt 1.3e-2 at
t=111 s (10x better than 10/100 at the same time, instability weaker/slower) but died at t=113 s from NaN in the OUTERMOST TOP
GHOST cell (bc_mode 3's WB walk under Prad-dominated low-density top state; a Newton e<=0 rescue preceded it). Fix of the ghost
walk + a 10/100 vs 100/300 comparison in progress.
bc3 baseline (11689919) final: t=106 s, dt 9e-4, dfloor 4.6e7 fofc 1.3e6 = same collapse as bc_mode 2; wall cleared.
MECHANISM of the He-star box instability (smoke_rt_3msun_bc3/he3_compare.txt): moving the blend edge 10/100 -> 100/300 REMOVES it
(h-KE 745x smaller, dt 0.1 s held to 50 s, box quiet; worst capped cell moves from tau 6-30 to the deep box, max x_i 700 -> 120),
moving the wall does nothing -> it is the CAPPED TRANSVERSE radiative diffusion (rad_cap_ang 0.5 throttles it ~1400x at tau<30,
x_i 700-1400) acting where the diffusion operator has weight; not a physical photospheric instability. The cap cannot be raised
explicitly (cap50 died at cycle 7), so the practical setting is 100/300 (diffusion only deep) and the real fix is an IMPLICIT
transverse operator. Second lid bug (exposed at 100/300): the two-stream Newton step leaves e<=0 in the optically thin
Prad-dominated top cell and the rt_rescue_eq 'radiative equilibrium' rescue sets T = 2.6e13-3.6e13 K -> dt collapse + NaN ghosts
(fix in progress). rt2turn (B star) re-staged with 100/300 + bc_mode 3, still unsubmitted.
bc_mode 3 ghost NaN FIXED (39ce5a5c): the polytropic WBAdvance walk with a ~1e13 K top cell (from the two-stream rescue bug) gave
d -> 0 by exp underflow and e = +inf (aT^4/rho), and +inf > 0 passed the positivity guard. Now: guard finite + positive + within
problem/wall_walk_maxfac (100) of the column ratio, else fall back to the rescaled mirror; ghosts floored. Conservation results
bit-identical to 4e0054e7; sick 1-D reproducer clean 2000 cycles. NEXT after the Newton-rescue fix: rebuild HIP, rerun 3 Msun at
100/300 + bc_mode 3, then the B-star rt2turn.
TWO-STREAM NEWTON FIXED (bcb1fb5f): the semi-implicit refinement's bare Newton on F(de) overshot to e<=0 in optically thin
Prad-dominated cells (iterates 3-4 decades past e; RG_v4/out.txt shows it twice, dhj never); now a step that would cross zero
energy bisects against the true lower bracket -(1-1e-3)e; bitwise where it does not fire (dhj ck + cs red giant byte-identical),
rad tests 4/4. The rt_rescue_eq rescue only COOLS (to 1e-3 e) and was not the source of T=2.6e13 K (that is the EOS inversion of an
over-pressured cell above the table ceiling). ALSO FOUND: hydro/tfloor UNSET in these inputs -> rt_use_cons clamps a non-positive e
to e(rho, FLT_MIN) = 1e-45; set tfloor = 5e3 K in the He v2 and B-star rt2turn inputs. New rt_efix_rec diagnostic prints the cell
state with the warning. Two agents edited the same worktree without conflict (disjoint files) - avoid in future.
HE 3 Msun v2 (bc3 + 100/300 + tfloor + all four fixes; job 11692296): FIRST CLEAN He-star RT run - t = 672 s (0.18 turnover) in
13 min, dt 0.109 -> 0.024 s (no collapse), mass exact, E +4e-4, no Newton warnings, event log 4 efloor + 164 fofc total. Profile
analysis pending (dt limiter, residual h-KE 9e30, FeCZ). 2-turnover continuation staged in hestar_fecz/he3_2turn (~75-90 min on
2 GPUs, apu chain), NOT submitted (user decision), like bstar_fecz/rt2turn.
v2 ANALYSIS (smoke_rt_3msun_v2/he3_v2.txt): HEALTHY. tau 6-10 turbulence gone (rms v_h 9.8e5 vs 8.4e6), h-KE saturated at
9e30 with no trend, no voids (min rho 91 dfloor), FOFC photospheric and decaying (783 -> 4 cells), top-cell tau 0.07 (vacuum BC
valid), T/T_Edd 0.99 above tau 10, emergent 0.94 F, E -0.6 % over 88 box energies. FeCZ convects: rms v1 1.5e6, v_h 9.4e5,
nabla-nabla_ad 55-75 % of the column retained. dt set by the top 3 cells (rho 8.7x below the column -> c_s up 2.5x; bounded
>= 2.5e-3 s); trimming x1max to tau ~0.3 buys 1.6x. CAVEAT for any result: a box-wide velocity floor of 50-90 v_MLT (rms 9e5,
Mach 0.12) even in the stable zone at tau>100, origin unknown (start-up transient? the still-capped deep transverse operator, max
x_i 120-240?) - the FeCZ signal is only 1.7x above it; needs a control before v_rms or F_conv/F is quoted. Cost: 27.5 cycles/s on
2 MI300A -> 2 turnovers ~3 h (he3_2turn staged, unsubmitted).
