---
name: solar-convection-test
description: Goal + status of tuning the solar_convection pgen to look like solar convection
metadata: 
  node_type: memory
  type: project
  originSessionId: 0df72474-f27b-49d2-9277-124d645c5193
  modified: 2026-08-05T15:51:30.136Z
---

Working on `src/pgen/solar_convection.cpp` (+ `run/sun/sun.athinput`). It is a 3D Newtonian
hydro convection box with a hand-rolled gray two-stream radiative transfer applied as a
`user_srcs` source term (`two_stream_RT`, a vertical column sweep). x1 = vertical/depth.

**Decision (user, firm):** do NOT use AthenaK's built-in radiation module (`prad`) — it is
GR-only, full speed of light (no RSLA), a bad fit for subsonic stellar convection. Keep the
two-stream RT framework untouched.

**Goal:** produce solar-like convection — sinking plumes + broad upflows in the convective
envelope, with a stable, ~hydrostatic atmosphere above. Reference: Stein & Nordlund.

**Diagnosed problem (from run/sun dumps):** original profile was strongly superadiabatic
(nabla=0.60 vs nabla_ad=0.40) with only a marginally stable transition; convective overshoot
pumped waves into the atmosphere that amplified as rho^-1/2 and reflected off the top BC,
giving vx_rms ~1.1 km/s at the top (velocity in the atmosphere). Mass IS conserved (the
history-file mass drop was a nan artifact at t=199, not real outflow).

**Changes made (2026-07-20):** (1) `get_wb_Tp` -> two-zone profile: mildly superadiabatic CZ
(nabla_cz=0.45) for p>=0.2*p0, isothermal stable atmosphere above (photosphere ~44% up the
box). (2) seed perturbation in the CZ (p>p_ph) at 1% amp. (3) fixed a var-shadowing bug in
SourceFunc cartesian branch (z was uninitialized) and added a top sponge layer (top ~22%,
tau=500s) that drains KE from total energy. Two-stream RT untouched.

**Run mechanics:** build per [[freya-build-procedure]] with PROBLEM=solar_convection. Test in
`run/sun_test/` (tlim shortened to 20000) via `sbatch submit_orion.sh` (partition p.exclusive,
16 MPI x 4 OMP). Analyze dumps with vis/python/athena_read.py; hydro_w vars = dens,velx(vertical),
vely,velz,eint(=p/gm1). Key check: horizontally-averaged vx_rms(z) should decay into the
atmosphere, not rise at the top. `/u/jinma` == `/orion/u/jinma` (same tree).

Status (2026-07-20): FIRST ITERATION WORKS. Job 179578 completed clean (tlim=20000, mass
conserved, no nan). Result: envelope convects at vx_rms~0.78 km/s; atmosphere quiet at
~0.066 km/s (z=9e7) and 0.12 km/s at top (vs 1.1 km/s in old run) — the atmospheric-velocity
pathology is gone. vz slices show plumes + overshoot in the envelope and a quiet atmosphere;
convection turns turbulent by t=20000. Weak solar-like up/down asymmetry (downflows faster,
smaller filling) emerging but limited by 64^3 resolution.

Iteration 2 (random white-noise seed, job 179582): convection more turbulent/irregular (good),
but atmosphere noisier — vx_rms rose from 0.066 (smooth) to ~0.19 (9e7) and ~0.47 km/s at top.

Diagnostic assessment of atmospheric velocity (assess.py): the atmospheric motion is internal
GRAVITY WAVES from convective overshoot, NOT convection (atmosphere strongly stable, N^2>0,
oscillatory sign-flips). Wave action v*sqrt(rho) DECAYS upward (dissipation) — not runaway
rho^-1/2. RT is a net DAMPER not a driver. Verdict: primary lever = top BOUNDARY (reflecting
-> resonant g-mode cavity); RT exonerated; resolution secondary. Much of the wave motion is
physical (real solar sims have it, handled by transmissive top / damping layer).

RT audit (two_stream_RT): numerical scheme CORRECT (short-char weights, quadrature, flux
divergence, indexing, units all verified). But found: (1) Tint4=0 => ZERO driving flux from
base, box only cools (secular), convection not sustained. (2) source was fully EXPLICIT
(implicit Newton block disabled) -> stiffness risk. Portability: NN=1040 stack arrays break on
GPU. Opacity kappa~rho^0.5 T^9 is fine (H- scaling).

Iteration 3 (job 179586, running): applied 3 changes user requested: (1) Tint4=SQR(SQR(5778))
turns on stellar driving flux F=sigma*Teff^4 at base; (2) re-enabled implicit RT update
(verified algebra: e_new-e_old=4pi kappa rho (J - B(T_new)) bdt, correct split); (3)
transmissive top BC = fmax(0, vertical mom) in outer_x1 ghost (outflow-only, less reflecting).
CAVEAT: initial photosphere ~4847K < Teff 5778K, so box heats/restructures over thermal time
~1e4 s; run shows relaxation toward balance, not instant steady state. Sponge still ON (confound
- do a sponge-off follow-up to isolate transmissive-top effect). Knob: Teff (5778 vs 4847 for
IC-consistent near-balance start).

RESULTS of driven runs (jobs 179586 implicit, 179592 explicit): BOTH bad. Implicit RT ->
box COOLED -10.5%, atmosphere ran away to ~560K (too cold). Explicit RT -> box HEATED +25%,
violent convection 1-1.5 km/s everywhere incl atmosphere. So RT source flips totE by ~35%
depending on time-integration = RT source far from balance, numerically fragile. No nan.

KEY PHYSICS INSIGHT (user, 2026-07-21) -- reason the driving was wrong: in EFFICIENT
convection the deep CZ carries energy by CONVECTION, not radiation. Box bottom is in the CZ,
so the upward RADIATIVE flux there should be tiny. Quantified: F_rad_bottom ~ (16 sigma T^3/
3 kappa rho)|dT/dz| ~ 1.1e8 vs sigma*Teff^4 = 6.3e10 => radiation carries only ~0.2% of total
flux at the base. So injecting sigma*Teff^4 as RADIATION at the base (Tint4=5778^4) over-supplied
by ~600x -> caused the pathologies. Tint4=0 (original) was actually correct for the radiation.
The ~sigma*Teff^4 energy must enter as CONVECTIVE ENTHALPY FLUX at the base, not radiation.
(Also: current bottom hydro BC is a frozen hydrostatic Dirichlet = fixed-T closed bottom, no
inflow; earlier Tint4=0 runs convected nicely but slowly cooled b/c no sustained energy source.)

PLAN FOR NEXT SESSION (agreed):
1. Revert Tint4 back to 0 (currently code has Tint4=SQR(SQR(5778)) - WRONG, revert first).
2. Add base energy input as CONVECTIVE flux, not radiation. Two options:
   (A) open fixed-entropy bottom (true Stein-Nordlund): upflows enter with fixed deep-adiabat
       entropy at prescribed pressure; downflows leave. Most physical, moderate effort.
   (B) [recommended first] bottom HEATING layer: keep closed bottom, add volumetric heating
       source in bottom few cells totaling sigma*Teff^4 per unit area (= the convective flux
       entering the box). Simple, one block in SourceFunc. Test via energy conservation:
       totE should PLATEAU, not swing +/-10-25%.
   Open Q for (B): spread heating over bottom ~5-10% of box vs just bottom couple cells.
3. RT explicit-vs-implicit matters much less once driving is balanced (RT source becomes small
   except near surface); implicit is safer choice for surface layers.

CURRENT CODE STATE (src/pgen/solar_convection.cpp): two-zone profile (nabla_cz=0.45, photosphere
p_ph=0.2 p0), random white-noise seed (amp 1% in CZ), top sponge ON (top 22%, tau=500),
explicit RT, transmissive top ON (fmax(0,mom) in outer_x1). Build: incremental make -j after edit.

SESSION 2026-07-21 (implemented the agreed plan): (1) Tint4 REVERTED to 0 in two_stream_RT
(radiation only cools; deep base carries ~0% radiative flux). (2) Added BOTTOM CONVECTIVE-FLUX
HEATING LAYER in SourceFunc: Fbase=sigma*Teff^4 (~6.3e10, Teff=5778) deposited as uniform
volumetric heating Qheat=Fbase/d_heat over the bottom f_heat=10% of the box (d_heat=0.1*(r1-r0)=
1.3e7 cm ~6.4 cells; Qheat~4863 erg/cm^3/s). Kernel: `if (x1v < z_heat_top) u0(IEN)+=Qheat*bdt`
(uses x1v not z so it's coord-agnostic). By construction column heating = Fbase = surface
radiative loss, so total energy should PLATEAU (key diagnostic vs the +/-10-25% swings of the
radiation-driven runs). Sponge unchanged (controlled single change from baseline).

BUILD NOTE: orion02/orion003 login nodes ARE Freya (share MPCDF soft + module system). Must run
the full module+cmake+make procedure in ONE bash call (shell state doesn't persist between calls;
default node compiler lacks -march=sapphirerapids). `module load openmpi/4.1` only resolves inside
the combined `module load gcc/13 openmpi/4.1 cmake/4.0` command, not step-by-step. Build succeeded.

RUN: job 179610 (submit_orion.sh, p.shared, 16 MPI x 4 OMP, tlim=20000) COMPLETED in run/sun_test/.
Previous bad run archived to run/sun_test/prev_179592_explicitRT_driven/.

RESULT of job 179610 (bottom-heating change): totE PLATEAUS at +25% (2.46e32) - SAME plateau as
the old explicit radiation-driven run 179592. Atmosphere still NOT quiet (vx_rms flat ~1.3-1.4
km/s most of the box). KEY LESSON: moving energy radiation->convective-heating changed NOTHING,
because with kappa~rho^0.5 T^9 the deep base is optically thick (tau~2e5) so the old radiative
Iint was absorbed within a few cells of the base anyway = already acted as near-base heating.
The injection CHANNEL was never the lever.

DECISIVE DIAGNOSIS (2026-07-21, scratchpad stability.py + tau.py, analyze last dump):
- tau=1 PHOTOSPHERE at x1=1.005e8 = 77% of box height, T=4688K. Below: optically thick (tau->2e5
  at base). Above (top ~23%): optically thin. RT IS making a photosphere + stable cap (physics OK).
- Schwarzschild: whole column CZ+lower-atm has grad~=grad_ad=0.40 (marginally superadiabatic,
  N2<0) => the box is a single near-adiabatic CONVECTING column. The flat vx_rms is REAL
  convection, not waves. Only the top ~11 cells (above tau=1, x1>1.05e8) are stable (N2>0).
- In that thin stable cap vx_rms DOES decay (1.18->0.84 km/s) - normal overshoot (Mach~0.15),
  NOT anomalously strong. Problem = only ~11 stable cells, overshoot enters at 1.2 km/s and
  can't decay to quiet before the top.
- ROOT CAUSE: photosphere MIGRATED UP from intended ~58% to 77% because the box OVER-HEATED
  +25% (too-cold IC, surface 4847K, relaxed toward Teff=5778K). Hotter box -> tau=1 pushed up
  -> atmosphere squeezed from ~27 cells to ~11.
- CONFOUND: top sponge (top 22%, starts x1=1.014e8~78%) OVERLAPS the entire stable atmosphere
  (tau=1 at 77%) -> sponge drains the whole stable layer; natural decay not cleanly visible.

ITERATION (2026-07-21, job 179689 RUNNING): EQUILIBRIUM IC (the real fix for the runaway).
Job 179617 (taller box) FAILED at equilibrium: totE never plateaued (+52% & climbing at t=20000),
CZ grew upward (photosphere 45%->64%), atmosphere re-filled with convection (21/53 atm cells
superadiabatic, vx_rms up to 1.9). ROOT CAUSE confirmed: NO thermal equilibrium - T(tau=1)=4854 <
Teff=5778, so surface radiates < base heating Fbase -> perpetual heating; box height only delays it.
FIX: rewrote get_wb_eos_arr to build the IC in RADIATIVE-CONVECTIVE EQUILIBRIUM (top-down integration
using the code's OWN opacity get_kapr): atmosphere = gray Eddington T^4=(3/4)Teff^4(tau+2/3) [so
T(tau=2/3)=Teff, emergent flux=sigma*Teff^4=Fbase]; interior = adiabat once Schwarzschild-unstable.
Shooting loop on p_top places photosphere at 50% box. Now get_wb_eos_arr/get_wb_eos take a 3rd array
Tarr(z) (tabulated T); get_wb_Tp is now UNUSED. Seed perturbations by height (x1v<0.45*box) not p_ph.
Offline profile_check.py verified: photosphere z=1.0e8 (50%), T=5776K; base T=13718K, rho=7.9e-7;
rho_ph=2.2e-7 (~solar!); skin T=4860K; base cs=17.8km/s -> dt~1.2s. Box starts BALANCED so it should
NOT heat/restructure. Prev taller-box run archived run/sun_test/prev_179617_tallbox/. CHECK: totE
should PLATEAU near 0% now (was +25/+52%); photosphere stays ~50%; atmosphere quiet+stable.
KEY OPEN Q if still churned: opacity norm C=100 vs Teff=5778 - rho_ph=2.2e-7 gives Mach~0.57
photospheric convection (realistic). Reusable: scratchpad/profile_check.py mirrors the C++ integration.

RESULT of job 179689 (equilibrium IC): IC CORRECT at t=0 (Tbase=13470, Ttop=4860 skin, balanced -
the construction works!) but the run CRASHED to nan at t~3380 (dt exploded 0.25s->doubling->15000s
once solution went nan). TWO failure modes seen in dumps t=0..3000:
(1) BASE OVERHEATS at startup: Tbase 13470->15185K in 1000s, Tmax->2e4, totE +15.5% by t=3000.
    Convection doesn't develop fast enough to carry Fbase, so the bottom heating layer overheats
    locally. (Low base density 7.9e-7 = low heat capacity -> heats ~6x faster than the old dense IC.)
(2) TOP ATMOSPHERE EVACUATES (the actual nan killer): rho_min 5e-9 -> 1.2e-9 -> 2.8e-12 (t=2000)
    while Ttop crashes 4860->2055->932K. The tall box puts ~4 scale heights of ultra-low-density
    atmosphere on top; RT cooling + the outflow-ONLY transmissive top BC (fmax(0,mom)) drain it to
    ~0 density -> negative density -> nan. |vx|max hit 9-14 km/s.
FIX OPTIONS for next iteration (NOT yet done - user to steer):
 (a) robustify low-density top: density/pressure floor, OR lower the box top (fewer scale heights
     of ultra-thin atmosphere), OR change top BC from outflow-only to hydrostatic (stop the drain).
 (b) base startup: spread heating over a thicker layer (raise f_heat) or RAMP Fbase in time so
     convection develops before the base overheats; bigger/faster seed.
 (c) two-stream RT != exact Eddington, so RT applied to the Eddington IC gives small net cooling in
     the optically-thin top -> may need to build IC from the code's own 2-stream closure, or damp RT
     in the thin top.
NOTE: initial totE now ~3e31 (10x lower than old 2e32) because equilibrium IC is low-density.

RESULT of job 179694 (equilibrium IC, HEATING OFF): STILL CRASHES the same way -> proves the crash
is NOT the base heating, it's the TOP ATMOSPHERE COLLAPSE. t=0 balanced (13470/4860). Heating off ->
base COOLS (13470->9289, totE -32%, expected). But top evacuates: rho_min 5e-9 -> 2.3e-10 -> 5.4e-16
(t=2000), |vx|max -> 35 km/s, nan at t=3787. Ttop RT-cooled 4860->1240K. Mechanism: ultra-thin top
(~4 scale heights, rho~5e-9) is over-cooled by 2-stream RT -> pressure drops -> gas falls + the
OUTFLOW-ONLY top BC (fmax(0,mom)) drains mass -> rho->0 -> negative density -> nan.
NEXT FIX (top-atmosphere robustness, heating OFF first to confirm survival, then re-enable heating):
 (1) TOP BC: outflow-only -> HYDROSTATIC extrapolation (stop the drain) -- likely the biggest lever.
 (2) density/pressure FLOOR to keep the thin top alive.
 (3) limit RT cooling in optically-thin top and/or lower box top (less ultra-thin fragile gas).
Heating currently DISABLED via `if (false && x1v < z_heat_top)` in SourceFunc. Reusable diag:
scratchpad/inspectdump.py (per-dump Tbase/Ttop/rho_min/vxmax/nan).

S&N BOUNDARY ITERATIONS (2026-07-21). User gave exact S&N (1998) spec: horiz periodic; vertical
transmitting. TOP (at temp minimum): larger zone, density gradient HYDROSTATIC, velocity gradient
ZERO, internal energy HELD CONSTANT in time+space at top fiducial layer. BOTTOM: outgoing fluid
leaves with its properties; incoming fluid - adjust PRESSURE so NET MASS FLUX through bottom
VANISHES (uniform over horiz plane, varies in time; = no boundary work on oscillation modes), damp
incoming-fluid velocity fluctuations (long time const), adjust density+energy of incoming fluid at
constant pressure to FIX its ENTROPY (space+time).
 - job 179695 (first S&N try, bottom pressure extrap from interior): RUNAWAY HEAT (totE +32000%,
   Tbase->95000K) - interior-tied extrap = positive feedback. CANCELLED.
 - job 179697 (fixed-reference bottom pressure p_ref, inflow density from entropy): MASS RUNAWAY
   (+110% by t=4000, Tbase fine ~14000) - fixing inflow DENSITY pumps mass in. CANCELLED.
 - job 179701 (RUNNING, faithful-ish): BOTTOM density+velocity ZERO-GRADIENT both ways (keeps net
   mass flux ~0 - the tractable stand-in for S&N's pressure-adjust-for-zero-flux), inflow (vx>=0)
   ENERGY set to fix entropy K_bot at continuous density: eint=K_bot*rho_is^gamma*igm1 (buoyant
   push when base cooled below deep adiabat = energy source, self-regulating). Outflow zero-grad.
   TOP: rho hydrostatic from interior (rho*exp(-dphi/(Rgas*T_top))), velocity zero-grad, specific
   internal energy CONSTANT at T_top_fix=4860K (temp min). Heating still OFF.
 DEVIATION from exact S&N (note for later): not doing the global horiz-uniform bottom pressure
   adjusted for exactly-zero net flux (needs global horizontal reduction across the 64 meshblocks
   tiling z=0; hard in a pgen BC), nor the incoming-velocity damping. Using zero-gradient density
   as the mass-conserving stand-in. If mass still drifts, implement the pressure controller.
 - job 179701 (zero-GRADIENT density both ways): opposite failure - whole box DRAINED to dfloor
   (rho=1e-11 everywhere, mass->0) in <1000s. No mass source -> downflows drain the open bottom.
   Proved BOTH extremes fail: fixed-density inflow pumps mass IN, zero-gradient drains OUT. The
   essential S&N ingredient is the ZERO-NET-MASS-FLUX condition. CANCELLED.
 - job 179702 (RUNNING/promising - FIRST S&N run that WORKS): bottom removes the horizontal-MEAN
   vertical momentum per meshblock (botmom(m) precomputed in a par_for over m; ghost IM1 =
   u0(IM1,is)-botmom -> net flux ~0), density zero-grad, inflow(mx_g>=0) entropy fixed to K_bot,
   top fixed-T. RESULT t=0..5000: ALIVE, CONVECTING, NO crash/runaway/drain. Top holds ~4900K
   (fixed-T works). Base cools 13470->9500K & LEVELS. |vx|max 10-17 km/s (transonic). Mass loss
   DECELERATES to plateau ~-41% (not collapse). totE ~-60% leveling. rho_max stable ~9e-7.
   job 179702 ran to t=20000 stable (archived prev_179702_SN_meanmom). But mean-momentum subtraction
   was NOT faithful S&N.
 KEY EMPIRICAL FINDING (jobs 179704 global-instant, 179706 global-relaxed): the PRESSURE/density
   approach (incoming rho_in from p_bot=K_bot*rho_in^gamma, zero-gradient velocity) DRAINS the whole
   box to dfloor in <1000s (total collapse, rho=1e-11 everywhere) REGARDLESS of how rho_in is set
   (instant algebraic sum_out/sum_in OR slow-relaxed). Zero-gradient VELOCITY cannot stop a coherent
   net downflow. Only removing the net vertical MOMENTUM stops the drain. So faithful "adjust
   pressure for zero net flux" alone fails here; must zero the net momentum directly.
 - job 179711 (RUNNING): FAITHFUL CO5BOLD (Freytag+2012) OPEN BOTTOM, per the exact recipe (found
   via websearch of arxiv 1110.6844 / ar5iv). Implemented as a SOURCE-TERM relaxation on the lowest
   active layer (is) at the END of SourceFunc using the stage dt bdt (user: "just use bdt for dt" -
   BC func has no bdt, so it lives in SourceFunc not the BC). 5 steps with GLOBAL horizontal means
   (Kokkos::parallel_reduce over bottom plane + MPI_Allreduce): (2) upflow(vx>0) entropy relax to
   s_inflow=cv*ln(K_bot): dr=rlx*(-r^2 T(g-1)/(P g))*ss, des=rlx*T*(1/g)*ss, rlx=CsChange*bdt*cs/dz,
   CsChange=0.1; (3) pressure damp to <P>: dr=rlxp*(1/cs2)*(<P>-P), des=rlxp*(1/(g r))*(<P>-P),
   rlxp=CPChange*bdt*cs/dz, CPChange=0.3; (4) drho4=(<rho0>-<rho2>) uniform density restore (mass
   cons/closed to plane waves); (5) v3-=<rho3 v3>/<rho0> (zero net mass flux). 3 passes: reduce
   <rho0>,<P>,N -> modify(2,3) -> reduce <rho2>,<rho2 v>,<v> -> modify(4,5). BC now just does
   isothermal-hydrostatic ghost extrapolation from the relaxed is layer. TOP = hydrostatic open
   (user reverted from fixed-T; NOTE it drained+crashed in 179708 - watch if it drains again here).
   Built OK. Prev: prev_179708_SN_globalmom_hydtop.
   RESULT 179711: CO5BOLD BOTTOM WORKS! Base ROCK-STABLE (rho 7.7e-7->7.9e-7, T 13470->12573 at
   t=4000), MASS CONSERVED (+0.3->+1.9% over 4000s, vs -41/-100% before!), totE stable ~-5%,
   convecting (|vx|max ~10 km/s). But CRASHES nan at t~5424 - failure is ENTIRELY THE TOP: upper
   atmosphere (76%+) cools (Ttop 4860->2660K) & drains (rho_min->2e-11) out the hydrostatic-open
   top, then the isothermal extrapolation goes singular -> nan (same as 179708). CONCLUSION: bottom
   SOLVED; the hydrostatic-open top (interior-tied energy) drains/crashes. FIX = revert top to
   CONSTANT-INTERNAL-ENERGY (fixed T_top_fix=4860) form = what BOTH S&N & CO5BOLD actually specify
   for the top ("hold internal energy at top fiducial layer constant") AND was stable to t=20000 in
   job 179702. Awaiting user OK to switch top back to constant-internal-energy.

*** WORKING CONFIG (job 179713, 2026-07-21): CO5BOLD bottom + constant-internal-energy top ***
User approved switching top back to fixed-T (T_top_fix=4860, hydrostatic density, zero-grad vel).
RAN FULL t=20000, NO CRASH (first S&N/CO5BOLD config to survive!). RESULTS:
- Convection: SOLAR-LIKE. vz slices show broad red upflows + narrow COHERENT blue downflow lanes
  from the corrugated surface (z~0.9e8=45%) down to the base. tmp slice: hot well-mixed CZ (>6000K)
  + sharp corrugated photosphere at ~0.85-1.0e8 + cooler stratified atmosphere (2000-4600K).
- Thermodynamics TEXTBOOK: deep CZ ~adiabatic (grad~0.45), SUPERADIABATIC PEAK grad-grad_ad up to
  +0.24 just below photosphere (z~0.65e8), sharp photosphere transition (grad crosses 0.4 at ~0.85e8,
  T~5300-5600K), then STABLE atmosphere: 0/53 atm cells superadiabatic (was 21-23!), N2>0 everywhere
  (mean 5.5e-4). THE STABLE ATMOSPHERE PROBLEM IS SOLVED.
- Stable, mass ~conserved. Photosphere stayed at ~45% (no upward migration).
REMAINING (refinements, not blockers):
- Mass grows slowly ~LINEAR +0.7%/1000s (+12.7% by t=18000). Small boundary imbalance, likely the
  fixed-T top: rho_ghost=rho_i*exp(-dphi/(Rgas*T_top_fix)) with T_top_fix=4860 > actual cooled top
  T~2600-3300 -> over-dense ghost -> slow mass influx. Fix: lower T_top_fix to the real atm top T,
  or tie it to interior T minimum.
- Atmospheric vx_rms HIGH ~1.6-2.7 km/s, does NOT decay, GROWS toward top (rho^-1/2 wave amplif +
  fixed-T top reflecting waves). Mach~0.2-0.3 (subsonic, ~realistic for chromosphere but strong).
  Fix: add a wave-damping layer near top, or a more absorbing top. Convective velocities: CZ 3.2,
  photosphere 2.9 km/s (Mach~0.18-0.26, realistic solar granulation).
Archived: prev_179711_CO5BOLD_hydtop. This run's outputs in run/sun_test/ (job 179713).

REFINEMENTS (2026-07-21, all on top of the working CO5BOLD-bottom + fixed-T-top config):
- Strengthened top wave-damping SPONGE (usrsource in SourceFunc): top 12%->20%, tau 500->50s.
  Cut reflected-wave component: top vx_rms 2.48->1.84 (-26%). Residual = intrinsic rho^-1/2 waves.
- T_top_fix 4860->3500 to reduce mass growth (only ~18% help; mass source mostly elsewhere).
- KEY PHYSICS on atmospheric waves: they are NOT over-driving (Teff=5790K=solar via tau=1 temp).
  They are convection-driven acoustic/gravity waves amplifying as rho^-1/2 in the thin atmosphere.
  Under-damped because kappa~T^9 COLLAPSES in the cool atmosphere -> RT decouples -> no radiative
  wave damping. Also: our EOS Rgas=1.38e8 => mu=0.6 (IONIZED); real photosphere is NEUTRAL mu=1.3,
  so our surface gas is ~2x too light (rho too low, P too high) - a fixed-mu ideal EOS can't match
  both ionized base & neutral surface (real codes use ionization EOS). rho_ph ~ (C*Rgas)^(-2/3).
- *** BIG WIN (job 179717, C=20 + OPACITY FLOOR): ***
  (1) get_kapr C: 100->20 (kapr=2.0e1*sqrt(rhom6)*Tp9) -> denser photosphere rho_ph~6.2e-7 ->
      v_conv~ rho^-1/3 -> CZ velocity 3.2 -> 2.2 km/s (~solar). Base state updated to C=20
      equilibrium: p_base=4.38e6, rho_base=2.31e-6 (in BOTH the BC constants AND the CO5BOLD
      s_inflow=cv*ln(K_bot) block). Verified via scratchpad/profile_check.py.
  (2) OPACITY FLOOR kapr=fmax(kapr,1.0e-2) -> keeps cool atmosphere radiatively coupled -> RT
      DAMPS the waves. Atmospheric vx_rms COLLAPSED: top 2.48->0.64 km/s (-75%!), now DECAYS into
      atmosphere (min 0.36) instead of rho^-1/2 runaway. Atmosphere still fully stable (N2>0).
  Floor barely shifts base state (photospheric kappa~0.11 >> 0.01 floor). Teff stays solar.
  Archived prev_179717_C20floor_64.
- RESOLUTION (job 179718 RUNNING): horizontal was UNDER-RESOLVED - dx_horiz=62.5km (nx2=nx3=64)
  vs dz=20.8km, 3x anisotropic; downflow lanes (~1 H_p=291km, ~100-300km) only 2-5 cells. Bumped
  nx2=nx3 64->192 (dx=20.8km CUBIC, ~solar-sim-like; matches Stein-Nordlund ~10-25km). meshblock
  nx2=nx3 8->16 (144 blocks). 3.5M cells (9x), ~9x runtime (dt already set by dz). No code change,
  athinput only.
  RESULT (job 179718, 192^3, t=20000, ~1.5-2hr): WELL-RESOLVED SOLAR GRANULATION. Horizontal
  photosphere slices (plot_photosphere.py, new script) show broad upflow granules (~1000-1500km)
  + sharp intergranular downflow-lane network reaching -4 to -6 km/s. Subsonic (M~0.6-0.7 peak).
  CZ vx_rms 2.2-2.4 km/s (~solar). Atmosphere STABLE (N2>0), waves DAMPED (vx_rms 1.7->0.55 decay,
  mild rise to 0.85 at top). FLAGS: Teff(tau=1)~5110K (~10% below solar) + mass +6.8%/totE +5.8%
  growth => box NOT fully equilibrated / slight energy imbalance (residual boundary mass-influx,
  never fully closed). Atmospheric v higher at 192^3 than 64^3 (more resolved waves).

VIZ TOOLS: plot_photosphere.py (run/) makes HORIZONTAL (x2-x3) granulation slices at the photosphere
(auto-located where horiz-mean T~Teff): vz (blue=up/red=down, seismic_r, +-6km/s) + T (gist_heat,
5000-6400K to show dark intergranule lanes). Saves to sun_test/plot/photosphere/. plotsun.py has
vz/tmp/tau/mach y-z vertical slices (vz +-2km/s, mach inferno 0-1.5).

OPACITY (get_kapr, user-updated 2026-07-21): composite = fmin(H- [20*rho6^0.5*(T/1e4)^9],
Kramers [3e3*rho5*(T/1e5)^-3.5]) + electron-scatter [0.34 if T>1e4] + floor[1e-2]. kapsc decl fixed
(was missing Real). INACTIVE in current box: H-/Kramers cross at ~24500K > base 13718K, so H- wins
everywhere; e-scatter negligible where H- large. Would only matter for a deeper/hotter box.
INFLOW ENTROPY CHECK (user asked re deep box): s_inflow=cv*ln(K_bot) set by PHOTOSPHERIC adiabat
(K=P/rho^gamma constant along convective adiabat) -> INDEPENDENT of box depth. Deep-box test
(scratchpad/deepbox.py, 1e9 box): K_bot/K_ph=1.0000, deep CZ stays fully convective (grad_rad>>
grad_ad even with Kramers), s_inflow=7.64e9=unchanged. So inflow entropy consistent, no change
needed. (User decided to KEEP original 2e8 box, not go deep.)

GRANULATION APPEARANCE (2026-07-21): user noted the tau=2/3-surface TEMPERATURE looks too turbulent/
filamentary vs real solar granules (smooth interiors + thin lanes). Diagnosis: NOT mainly a slicing
artifact (tau=2/3 surface, smoothed by linear-tau interp in plot_tau_surface.py, still filamentary).
Real causes ranked: (1) ideal-gas EOS has NO IONIZATION -> no H recombination latent-heat buffering
of T fluctuations (biggest). (2) two_stream_RT is column-by-column 1D (2 angles just scale dtau/mu,
no HORIZONTAL radiative transport) -> at optically-thin tau~1 real 3D RT smooths horizontal T
contrasts, ours can't (strong #2, acts right at the imaged surface). (3) GRAY (no opacity binning)
-> diffuse not sharp surface cooling. Also possible: explicit RT + the kappa floor may add thermal
noise (cheap to test via implicit-RT/no-floor). Velocity/morphology are solar-like; T smoothness is
the toy limitation.

SESSION 2026-08-05 -- COOLING-DRIVEN-CONVECTION EXPERIMENT SET UP.
(1) FOUND+FIXED A REAL BUG in solar_convection.cpp: get_wb_eos_arr had `grad_ad =
    0.8*(gamma-1)/gamma` = 0.32 (comment said 0.4). That 0.8 factor made the IC interior
    SUB-adiabatic AND, because the hard-coded bottom-BC base state (p_base=4.38e6,
    rho_base=2.31e-6, T=13718, K=1.085e16) had been computed with grad_ad=0.4, left the
    boundary injecting Delta_s/cv = ln(1.087e16/7.99e15) = 0.31 of EXCESS ENTROPY at the
    bottom -- comparable to the entropy contrast that drives solar convection, i.e. the
    bottom was driving the convection itself. Likely source of the +6.8% mass / +5.8% totE
    growth in job 179718. Reverted to `grad_ad = (gamma-1)/gamma`; base state then matches
    the hard-coded constants to <1.5% (Delta_s/cv = 0.0065, negligible) so NO BC edit needed.
    LESSON: the IC base state and the bottom-BC constants must be recomputed together.
(2) NEW PGEN src/pgen/cooling_convection.cpp (copy of solar_convection.cpp) for the
    experiment. Changes vs the original:
    - get_wb_eos_arr now takes gamma and grad_int as ARGUMENTS (gamma no longer hardcoded
      5/3 inside; user preference: take gamma from the EOS/problem setup, not constexpr).
    - The two uses of grad_ad are SPLIT: `grad_switch = (gamma-1)/gamma` is the Schwarzschild
      threshold (fixes the photosphere depth), `grad_int` is the gradient actually imposed
      below it. grad_int < grad_switch => sub-adiabatic stable interior WITHOUT moving the
      photosphere. Runtime knob `<problem>/grad_int`, defaults to (gamma-1)/gamma (= exactly
      reproduces solar_convection).
    - Base state now COMPUTED from the IC integration into file-scope p_base_ic/rho_base_ic/
      T_base_ic (anonymous namespace) and read by BOTH the bottom BC and the CO5BOLD
      s_inflow=cv*ln(K_bot) relaxation -- kills the stale-constant bug class permanently.
      Computed BEFORE `if (restart) return;` (restart returns before gamma is known), using
      a local `MeshBlockPack *pmbp_ic = pmy_mesh_->pmb_pack` (pmbp not yet in scope there).
      Needs #include "globals.hpp" for global_variable::my_rank in the printout.
    - Prints "cooling_convection: grad_int=... base state: T=... rho=... p=... K=..." at
      startup -- use it to confirm IC/BC consistency.
    VERIFIED: builds clean; grad_int=0.4 -> T=13717.6 rho=2.28205e-6 p=4.31998e6 K=1.09214e16
    (matches the offline mirror exactly); 0.32 -> T=12134.6 rho=3.01056e-6 K=8.032e15;
    0.30 -> K=7.384e15; 0.36 -> K=9.416e15.
    GOTCHA: AthenaK requires a cmdline-overridden param to EXIST in the athinput, so
    `grad_int` must be added to the <problem> block before `problem/grad_int=...` works.
    NOT YET DONE: no athinput/run dir for this pgen yet; seed perturbation still 1% below
    0.45*box (consider weakening so the result is unambiguously cooling-driven); bottom
    heating layer still disabled (keep it that way).
    build/ currently holds the cooling_convection executable (clean rebuild wiped the
    solar_convection one).
PHYSICS FRAMING for the experiment: the mechanism is convective ENTRAINMENT / penetrative
erosion of a weakly stable layer -- surface-cooled plumes sink until their entropy matches
ambient, mix, flatten the gradient, and the mixed layer grows DOWNWARD with time. That
growing mixed-layer depth is the diagnostic, not instant full-box convection. Box thermal
time ~6500s so tlim=20000 is ~3 thermal times; expect totE to DIP then recover (that dip-
then-recover IS the result). NOTE a literal "Eddington radiative everywhere" IC is NOT
usable: with kappa~rho^0.5*T^9 the Eddington relation drives grad=dlnT/dlnp up to ~4, base
T=52700K, and DENSITY INVERTS (rho_base 1.5e-7 < rho_photosphere 6.2e-7) => Rayleigh-Taylor,
not just Schwarzschild, unstable. Capped-radiative variants (grad capped at 0.44-0.80) stay
density-monotone if that route is ever wanted. Reusable: scratchpad/rad_profile.py mirrors
the C++ integration incl. the current composite opacity.

*** RESULT (2026-08-05, job 185393, run/cooling_test/): COOLING DOES DRIVE THE STABLE LAYER
CONVECTIVE. Clean POSITIVE answer, via textbook convective ENTRAINMENT. ***
IC verified stable at t=0: grad = 0.314-0.339 at every height (grad-grad_ad = -0.06..-0.09),
entropy s/cv monotonically increasing upward (36.644 base -> 37.015 at 55%).
Sequence (horizontally-averaged grad = dlnT/dlnp; photosphere at 50% = 1.0e8):
 - t=500 : instability NUCLEATES at z~35% (just BELOW the photosphere), grad jumps to 0.654
           (+0.25 superadiabatic) while everything below is still at the initial -0.08.
           Entropy there drops 36.811 -> 36.435 = the radiative-cooling entropy deficit.
 - t=1000: unstable region has eaten DOWNWARD to 15-25%; base still stable.
 - t=1500: reaches the base (5% superadiabatic). Mixed-layer base tracked: 0.51e8 (t=500) ->
           0.20e8 (t=1000) -> 0.03e8 (t=1500+).
 - t=2500: CZ WELL MIXED, entropy nearly flat 36.44-36.58, broad mild superadiabaticity
           (+0.05..+0.1). Atmosphere above stays strongly stable (grad-grad_ad ~ -0.38).
ENERGY confirms independently: totE DIPPED to min 8.50693e31 at t=1600 (-9.30%) then RECOVERED
(8.60005e31 at t=2800) -- the dip-then-recover predicted for a stable box that has to spin up
convection. Turnaround coincides with the mixed layer reaching the base, i.e. when the bottom
BC first has something to push against. CAUSE IS CLEAN: at t=500 (when destabilisation happened)
the base entropy was still 36.643 = its initial value, so the boundary had done nothing.
The mixed CZ settles slightly BELOW s_inflow (36.5 vs 36.62) => the bottom continuously pushes
entropy up = self-regulating, responding to the convection rather than driving it.
Mass bottomed -0.13% (t=1000) then rose steadily to +0.46% (t=2800) - small but a TREND, watch it.
CZ top sits at z~0.78e8, BELOW the t=0 photosphere at 1.0e8 (tau=2/3 surface will have moved as
the box cooled ~9%; recompute tau before reading anything into that gap).
NEW VIZ: run/plotcool.py (vz x2-x1 slices; writes BOTH fixed +-2km/s in plot/vz/ for
comparability with solar_convection AND auto-scaled 99.5th-pct in plot/vz_auto/ -- the auto one
is essential early, when a stable box's velocities are far below the 2 km/s granulation scale;
prints |vz|max + vrms at base/photosphere/top per snap). run/plotcool_entropy.py (s/cv and
grad-grad_ad vs z, all snapshots coloured by time, marks s_inflow; prints mixed-layer base).
NOTE plotsun.py's tau panel uses the STALE opacity (1e2, no floor, no Kramers) - don't trust it.

*** KEY OPACITY SCALING RESULT (2026-08-05, scratchpad/opscan.py) -- I WAS WRONG ONCE, THIS IS
THE CORRECTED VERSION. *** Question: can a genuinely RADIATIVE-EQUILIBRIUM envelope be stable
here, so we can test whether surface cooling drives convection WITHOUT the profile simply
relaxing to the convective state it should have had?
 - LOWERING THE OPACITY NORMALISATION C DOES NOT WORK, and I initially suggested it in error.
   Analytics: nabla_rad = kappa*p/(4g(tau+2/3)), and hydrostatic-in-tau gives p ~ g*tau/kappa,
   so the PRODUCT kappa*p is set by the tau scale and is INDEPENDENT of C. Confirmed numerically:
   C = 20, 2, 0.2, 0.02 ALL cross grad=0.4 at z=1.00e8 (the photosphere) and blow up below.
   rho_ph ~ C^(-2/3) and kappa_ph ~ C^(2/3) cancel exactly.
 - THE REAL LEVER IS THE TEMPERATURE EXPONENT, not the normalisation. The instability comes
   entirely from kappa RISING with depth (rho^0.5 T^9).
 - CONSTANT OPACITY IS THE DECISIVE TEST: gray Eddington + constant kappa gives the classic
   nabla -> 1/4 result. Numerically max grad in box = 0.234 < 0.400, NEVER unstable at any
   kappa0. T_base = 9688.8 K identical for every kappa0 (T structure depends only on tau, and
   tau(z) is pinned by putting the photosphere at 50%); rho_base ~ 1/kappa0:
     kappa0=1.0 -> rho_base 2.025e-7, K_bot 3.877e16, p_top 1.683e2
     kappa0=0.1 -> rho_base 2.025e-6, K_bot 8.354e15, p_top 1.683e3   <-- rho closest to the
                   current run's 3.01e-6, so similar dt; RECOMMENDED starting point
     kappa0=0.01 -> rho_base 2.025e-5, K_bot 1.800e15
 - PLAN for the test: set get_kapr to a constant kappa0, build the IC radiative-equilibrium
   EVERYWHERE (no adiabat switch, i.e. the radiative_only branch), match the bottom entropy as
   usual (it is now automatic - base state is computed from the IC), then see if surface cooling
   still drives convection. PREDICTION: brief adjustment transient, then NOTHING sustained,
   because in radiative equilibrium dF/dz=0 so there is no net cooling to make an entropy
   deficit. If it DOES convect that is the real result and far more interesting.
*** CONTROL RUNS COMPLETE (2026-08-05) -- DEFINITIVE NULL. Both ran the FULL t=20000, 41 dumps,
no nan, ZERO superadiabatic cells summed over every dump. ***
 - job 185888 radiative_test (kappa=const 0.1): worst max(grad-gad) = -0.043 (t=2500).
   mass -2.0%, totE -13.2%. final KE 4.53e28/1.64e28/1.53e28.
 - job 185988 powerlaw_test (kappa=0.1*rho^0.5*T^0.5, b=0.5): worst = -0.030, AT t=0 -- i.e. it
   never got closer to neutral than its own IC. mass -1.0%, totE -6.4%. KE 7.32e28/3.00e28/2.50e28.
   At t=500, the exact moment cooling_test spiked to +0.569 and nucleated, this run was at -0.086
   (MORE stable). Settles to a stationary -0.033 after t~1500.
 - cooling_test (185393) finished too: t=20000, steady statistical equilibrium from t~2000 (base
   s/cv 36.577+-0.005 for the last 18000s, mixed-layer base pinned at 0.03), sustained convection,
   BUT mass ended +8.9% (leak grew all run; NOT the entropy mismatch - suspect fixed-T top ghost).
 - powerlaw_test is the CONTROLLED pair with cooling_test (base state matched: T 12055 vs 12135,
   rho 3.13e-6 vs 3.01e-6, K 7.78e15 vs 8.03e15) => the two runs differ ONLY in whether the
   stratification can carry the flux radiatively. b=0.5 margin (9%) held; b=0.8 would likely
   have crossed.
NEW OPACITY API in cooling_convection.cpp: get_kapr(rho,T,kap0,kapa,kapb,kapr) - power law
kap0*(rho/1e-6)^kapa*(T/1e4)^kapb when kap0>0, else the composite H- law. kapa=kapb=0 =
constant opacity (backwards compatible). Params <problem>/kappa_const, kappa_a, kappa_b,
radiative_ic. Threaded through 3 call sites + get_wb_eos_arr; two_stream_RT takes LOCAL copies
(kap0/kapa/kapb) because device lambdas cannot read host globals.
NEW RUN DIRS: run/radiative_test/ (rad.athinput, athena_rad), run/powerlaw_test/ (plaw.athinput,
athena_pl). Both PIN A PRIVATE COPY of the binary (cp build/src/athena into the run dir) --
IMPORTANT because relinking build/src/athena while a job runs from it is a real hazard; it only
survived because the linker made a new inode. Always copy the binary into the run dir.
VIZ: run/plot_entropy_compare.py now takes [right_run] [tmax]; right_run = radiative_test |
powerlaw_test. Default tmax = the SHORTER run's end time so both share the colour scale.
Writes cooling_test/plot/entropy_compare_<right_run>.png.

*** BC DEFECT FOUND (2026-08-05) -- CO5BOLD inflow-entropy relaxation is GATED ON UPFLOWS. ***
cooling_convection.cpp ~L1296: `if (vx > 0.0) { ... relax entropy toward s_inflow ... }`.
Velocities start at EXACTLY zero and the seed perturbs DENSITY only, so the gate is false at the
base and the relaxation never fires. The base then cools unopposed, the entropy deficit
ACCUMULATES, and the moment any upflow appears the whole stored deficit discharges at once.
Measured in radiative_test: s_base-s_inflow went 0.006 -> -0.461 (t=1000) -> -0.674 (t=2000)
with v_rms(base) EXACTLY 0.000 the whole time, then at t=2500 snapped to -0.215 with a 1.569 km/s
kick, thereafter pinned at ~-0.163 with ~1.3 km/s sustained. powerlaw_test fired earlier/gentler
(already pinned at -0.117 by t=500).
HOW IT WAS DIAGNOSED (reusable): (1) v_rms(z) shows the motion is BOTTOM-concentrated, 1.4-1.6
km/s at z=0.03-0.10 decaying to 0.1 by z=0.75; (2) the enthalpy flux <rho*vx*cp*T'> is NEGATIVE
(-3 to -4e9, ~5% of sigma*Teff^4) at the base = COUNTER-GRADIENT = mechanically forced motion in
a STABLE layer, NOT convection (convection carries heat UP).
Does NOT invalidate any conclusion: the null runs have zero unstable cells and this flow moves
heat the wrong way; cooling_test nucleated at z=0.35-0.45 at t=500 while its base was still at
the initial entropy. But z<0.1 in the radiative runs is boundary-controlled - do not read physics
there. FIX (not yet done): apply the relaxation whenever s < s_inflow regardless of flow
direction (closest to the physical intent), or weight it smoothly by max(vx,0)/|v| instead of a
hard on/off. Check steps 4/5 of the same block at the same time - they do the mass/momentum
adjustment and may be behind cooling_test's +8.9% mass growth.

DENSITY INVERSION vs CONVECTION (physics settled 2026-08-05): convection ALWAYS sets in before a
density inversion can develop, and this is provable. rho ~ p/T so dlnrho/dlnp = 1 - nabla;
inversion needs nabla > 1, convection needs nabla > nabla_ad = (gamma-1)/gamma < 1 for ANY
gamma>1, and nabla_rad varies continuously, so 0.4 is crossed first. Equivalently:
N^2 = -g^2 rho/(Gamma1 p) - g dlnrho/dz, whose first term is always negative, so ANY hydrostatic
density inversion has N^2 < 0 automatically. Inversion is SUFFICIENT for convective instability,
never a competing condition.
BUT convection starting does NOT prevent the inversion. Efficient convection (dense, subsonic,
high heat capacity - our box, solar envelope) drives nabla_actual back to ~nabla_ad, so no
inversion. INEFFICIENT convection (low density, near/super-Eddington, e.g. the iron-bump layers
of luminous-star envelopes) carries little flux, nabla_actual stays ~nabla_rad > 1, and the layer
is density-inverted AND strongly Schwarzschild-unstable at the same time. Caveats there: Ledoux
can stabilise if there is a mu gradient (usually uniform composition in those envelopes, so
Schwarzschild=Ledoux); and MLT is invalid in the radiation-dominated regime, so the 1D inverted
hydrostatic profile is largely a 1D artifact - 3D gives porous/clumpy structure and outflows.
NOTE a framing correction I had to make: Schwarzschild instability and "Rayleigh-Taylor" are NOT
two distinct instabilities in a compressible hydrostatic atmosphere - N^2 < 0 is the single
buoyancy criterion and a density inversion is just a badly-violated case. Difference is
quantitative (growth rate/vigour), not mechanistic.

CAUSALITY (user's closing question, settled): convection is ALWAYS driven proximately by unstable
stratification (N^2<0); "cooling drives convection" is shorthand for cooling CREATING that
unstable stratification. Two routes to instability: (a) OPACITY - nabla_rad > nabla_ad as a
property of the steady state (solar envelope, H-/ionisation); (b) COOLING - a top-concentrated
heat sink generates a local entropy deficit that did not exist before (our cooling_test, whose IC
was stable at EVERY height). The control runs prove cooling ALONE is not sufficient: same
photosphere, same sigma*Teff^4 radiated, but a stratification that CAN carry the flux
radiatively never destabilises. Third, separate sense of "driven" (Stein & Nordlund): where the
KE/buoyancy work is generated - the solar CZ is unstable throughout by route (a) yet energetically
driven from the surface by tau~1 cooling. Our steady CZ shows both: nabla-nabla_ad ~ +0.05..+0.10
maintained by ongoing surface cooling against the bottom entropy input.

PHYSICS FRAMING (user asked "doesn't this mean a star/planet with a radiative surface can have
surface convection from cooling?"): the operative variable is NET COOLING vs THERMAL
EQUILIBRIUM, not "radiative surface". Our IC was NOT in radiative equilibrium (grad=0.32 cannot
carry sigma*Teff^4), hence real flux divergence -> cooling -> entropy deficit. A star in true
radiative equilibrium has dF/dz=0: the surface radiates but is exactly resupplied, nothing cools,
no convection. Genuine real-world YES cases all involve unbalanced/time-varying forcing with
cooling concentrated at the TOP of the layer (cooling at the BOTTOM stabilises - that is why
nocturnal ground cooling makes an inversion): ocean deep convection (Labrador/Greenland Sea,
stably stratified, driven km-deep by surface buoyancy loss alone - closest analog to our run),
nocturnal ocean mixed layer, stratocumulus cloud-top radiative cooling. NO cases: A-star /
hot-WD radiative envelopes. And solar surface convection is NOT an example either - there the
deep envelope is Schwarzschild-unstable on its own via H-/ionisation (nabla_rad > nabla_ad);
cooling shapes and drives it but does not create the instability.

DEFERRED EXPERIMENT (2026-08-04, superseded by the 2026-08-05 entry above): "does COOLING alone drive convection
from a RADIATIVE initial condition?" User wants to build the IC stratification RADIATIVE EVERYWHERE
(not adiabat below the photosphere) and see if the Schwarzschild instability + surface cooling spins
up convection on its own. STATUS: edit NOT YET made -- get_wb_eos_arr (solar_convection.cpp ~L1315-1341)
STILL has the Schwarzschild->adiabat switch active (interior adiabatic, = working run 179718). NOTE
solar_convection.cpp is UNTRACKED in git (no diff available). TO DO: (1) make interior use Eddington
T_rad=Teff*(0.75*(tau+2/3))^0.25 all the way down (never set convective=true) -> strongly
superadiabatic deep IC (T~tau^0.25), violent t=0 transient expected. (2) CRITICAL consistency fix:
radiative IC base (z=0) is HOTTER/steeper than adiabat, so it will NOT match the hard-coded base
state (p_base=4.38e6, T_base=13718, rho_base=2.31e-6, K_bot) used by the open-bottom BC (L733-736)
AND the CO5BOLD inflow-entropy relaxation (L1185, s_inflow=cv*ln(K_bot)). Must recompute the new
base state (via scratchpad/profile_check.py which mirrors the C++ integration) and update BOTH places,
else IC<->boundary mismatch -> the documented runaway/drain failure. (3) Judge result in EARLY
TRANSIENT: the entropy-fixing bottom re-imposes its own adiabat over a thermal time regardless of IC,
so late-time convection is partly the boundary, not the cooling.

NEXT-SESSION: discuss a CODE-modifying project (deferred). Candidates surfaced this session:
(A) realistic IONIZATION EOS (tabulated H/He) - highest physics impact, replaces fixed-mu ideal gas;
(B) generalize the hardcoded CO5BOLD open-bottom + const-energy top into a reusable athinput-config
BC module; (C) improved RT (multi-angle / opacity-binned / eventually 3D horizontal transport).
User wants to discuss next time (did not pick yet).

 - job 179708 (drained via TOP, archived): GLOBAL mean-momentum removal (single value: reduce SUM(rho*vx)
   and cell count over whole bottom plane + MPI_Allreduce -> mom_mean; ghost IM1 = interior-mom_mean)
   + fixed-entropy inflow. TOP reverted to HYDROSTATIC OPEN (user request): isothermal-hydrostatic
   density+energy extrapolation from interior (q0_i/factor_i form), zero-gradient velocity (open both
   ways) - dropped the fixed-T_top-minimum form. CHECK mass conserved? survives? convecting?
 - job 179704 (RUNNING, FAITHFUL global p_in): user insisted p_in be a SINGLE value summed over
   the WHOLE bottom plane (uniform over horizontal plane, per S&N), not per-meshblock. Implemented:
   Kokkos::parallel_reduce over all bottom-face cells -> local (sum_out=SUM rho|v| for vx<0,
   sum_in=SUM v for vx>=0), then MPI_Allreduce(2 reals, MPI_SUM) across ranks -> ONE global
   rho_in=sum_out/sum_in (clamped [0.05,20]*rho_base), p_bot=K_bot*rho_in^gamma. Incoming cells
   all use this single rho_in; outgoing zero-gradient. Added #include <mpi.h> (guarded). This is
   S&N's exact "adjust pressure so net mass flux vanishes, uniform over horizontal plane". CHECK:
   mass conserved (should be MUCH better than per-mb -41%)? survives? convecting? atmosphere quiet?

Current code state: athinput <hydro> floors dfloor=1e-11, tfloor=500K. Heating DISABLED via
if(false && x1v<z_heat_top). BC constants in HydrostaticEquilibrium: gm1, K_bot=p_base/rho_base^g
(p_base=1.5e6, rho_base=7.92e-7), T_top_fix=4860. Ghost index renamed `igc` (outer `Real ig`=1/gamma
is the pow exponent). Archived runs: prev_179689_eqIC_heat, prev_179694_eqIC_noheat,
prev_179695_SN_runaway, prev_179697_SN_massgrow (all in run/sun_test/).

(SUPERSEDED) ITERATION job 179617: "give the atmosphere ROOM" via a TALLER BOX.
Insight: photosphere settles at a fixed PRESSURE (tau=1 at p~1.2e6) => fixed HEIGHT (~1.0e8) set
by stratification from the base up, INDEPENDENT of where the top is. So extending x1max 1.3e8->
2.0e8 keeps the photosphere at ~1.0e8 but moves it from 77% of box (11 cells above) to ~50%
(~48 stable cells above). This makes the IC-consistency change UNNECESSARY (photosphere migration
stops mattering with that much atmosphere above), so IC/Tbot left UNCHANGED (clean controlled test).
Changes: (a) sun.athinput mesh+meshblock nx1 64->96, x1max 1.3e8->2.0e8 (dx1~2.08e6, ~unchanged;
64 meshblocks, 16 ranks x 4). (b) sponge top 22%->12% (zb=r1-0.12*(r1-r0)) so it no longer overlaps
the stable atmosphere (was a confound draining the whole layer). Heating layer stays f_heat=10%
(column total = Fbase automatically, balances surface cooling regardless of box height). Built OK.
Prev bottom-heat run archived to run/sun_test/prev_179610_bottomheat/. Diagnostics reusable:
scratchpad/{assess,stability,tau}.py (point at run/sun_test/bin). Check: does the atmosphere quiet
(vx_rms decay above photosphere ~1.0e8) now that there are ~48 stable cells + sponge off the layer.
NOTE dx1 grid indices shift: photosphere now ~cell 48 of 96; heating layer bottom 10 cells.
