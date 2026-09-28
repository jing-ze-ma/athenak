---
name: dhj-improvement-list-0921
description: Ranked to-do list for the cubed-sphere deep hot Jupiter (dhj) work, agreed 2026-09-21 night; status of each item
metadata:
  type: project
---

User 09-21 night: "cross out the undone improvement list sorted by importance". Ranked list + status:

1. [DONE 09-21: NO DEFECT, note retracted; open instead: bvals_cc corner slots 48-55] bvals_fc cs SEAM halo with > 1 MeshBlock per panel (every cs MHD production ran 2x2) -- agent
   running 09-21 night: tests_seam_fc/, build_cs_fc, UNCOMMITTED; the cc fix 94c7165d is already in HEAD.
2. [BUILT + committed 09-21 night as problem/ck_spherical, DEFAULT OFF; see RESULT below] r^2 SPHERICAL FORM missing in the CORRELATED-K kernel and the monolithic rt_split=false sweep
   (= what dhj runs). CONFIRMED from the code's own comment (two_stream_rt.hpp ~1790: "not converted ...
   r_out/r_in < 1.3") while the production grid is x1 9.44e9..2.0556e10 = r ratio 2.18, area 4.7.
   Agent running: problem/ck_spherical (default off), tests_ck_sph/, build_cksph_*.
3. [DONE 0eb7e2c4 for Cartesian/sp; cs gap (RaiseVel re-solves caches) agent running: tests_cs_rst/] general-EOS MHD restart caches (wder/wtemp; hydro has them since 6ae9ccbe) -- agent running:
   tests_mhd_rst/, build_mhdrst_*.
4. ck / deep-conduction HANDOVER at ck_pcut_bar=10 bar (rad_tau_lo/hi 30/300) is operator split and
   per-cell semi-implicit = the structure that gave the dt-proportional pump in the He/B boxes;
   test = a half-CFL arm (needs a GPU job: ask the user before submitting), then decide on a merged ck
   column solve (the ck mode-3 Jacobian was dropped at the merge by user choice).
5. wb_cache_every = 10 in sp_mhd_prod3 and cs_mhd_prod2 (the WB arms): the STALE WB CACHE pump
   ([[fmode-wb-cache-culprit]]) was active there -> candidate for the unexplained WB-arm death at rot
   12.3; any new WB arm must use wb_cache_every = 1.
6. dt_min = 1e-2 too tight for the semi-implicit setting (use 1e-3).
7. implicit transverse conduction on cs (remove rad_cap_ang = 0.5): docs/dev/cs_implicit_transverse.md.
8. floors at the night-side vacuum top (dfloor vs WB background, keep-temperature).
DONE/moot: rt_rad_force was never on in cs_mhd_prod3 / sp_mhd_prod3 / cs_mhd_prod2 -> the transverse
force bug 4bcdc855 never affected dhj. cs_mhd_prod3 itself ran clean to rot 283 (STOP file 09-14).
See [[m1-paused-back-to-dhj]], index-hot-jupiter.md.

RESULT ck_spherical (tests_ck_sph/README.md, 1-D CPU column on the production radial grid, NOT yet a 3-D
GPU run): the OLD form is non-conservative by O(1): column budget error +0.63, beam power deposited
0.378 of incident, per-cell identity 1.8e-2 (apply kernel uses areas, propagation does not dilute = the
"neither form alone is conservative" case -> spurious net sink ~ -2F/r). NEW: 1e-10 / 2.6e-9 / 0.991.
Effect after 3000 cycles: e/rho +13 % mean (+42 % max) at 1e-3..1 bar; above 1e-4 bar the OLD form
COLLAPSES ONTO dfloor (rho = 5e-14 exactly, e/rho down 40x), the new one stays smooth -> prime suspect
for the night-side vacuum top slab / FOFC / Alfven-dt problems of every dhj production. Cost 1.46x.
ALSO FOUND: the semi-implicit source apply swallows 36-40 % (old) / 2-5 % (new) of what the sweep asks
for (rt_desum, now reachable from the dhj pgen) -> next item. rad_flux_inner for ck = sigma T_int^4 at
the CUT face (per unit area there), dilutes 3.83x to the top: flagged, not rescaled. All dhj productions
to date (cs_mhd_prod3 to rot 283, sp_mhd_prod3, ...) carry the old form. NEXT: a GPU production-grid
A/B needs user approval. Also done 09-21: restarts bitwise for tabulated-EOS MHD (0eb7e2c4) and on the
cubed sphere (544cb36b); cc corner slots 48-55 agent running (tests_seam_cc_corner/).

UPDATE 09-21 late night -- RETRACTION + QUEUE (user-approved order):
* RETRACTED: the BEAM part of ck_spherical as committed in 3c5d4846 is WRONG (my brief): it made each
  column absorb the power intercepted at the domain TOP, A_top mu0 F* (planet absorbs pi r_top^2 F*,
  ~3x too much; the "0.991 vs 0.378" and the "+13 %/+42 % warmer, no dfloor collapse" 1-D result are
  contaminated by it). A parallel beam is not confined to a diverging column; the OLD local deposit
  kappa rho F* exp(-tau) is the right one. The THERMAL spherical form is right. The ck agent was resumed
  to restore the beam, add gates "L(r) = A F_net flat at A(x1min) sigma T_int^4 with the star off" and
  "A_top F_top,thermal = absorbed beam + internal" (user's request), and redo the e/rho comparison.
  ck_spherical stays default OFF; nothing ran with it.
* QUEUE on the same kernel, one after the other: (1) beam correction + gates [running]; (2) option B =
  PSEUDO-SPHERICAL direct beam (true chord lengths through the shells at impact parameter b = r sin
  theta0, column-local/horizontally homogeneous, both legs through the tangent point => TWILIGHT past the
  terminator, up to ~55 deg at the domain top; chord matrix precomputed once since the star is fixed in
  the grid; default off) [user: yes]; (3) REAL IMPLICIT ck column solve = restore the mode-3 Jacobian
  for ck (dropped at the merge; recover from pre-merge history), one tridiagonal in T per column summed
  over g-points, once per step, conduction folded in (removes the 10 bar split handover) [user: yes].
  Gate for (3): rt_desum sweep-to-gas gap -> round-off (today 2-5 % new / 36-40 % old geometry).
* GPU: HEAD did not compile with hipcc (DualView explicit-space templates in 9 files) -> fixed b4c6a2f2.
  A/B smoke STAGED, NOT SUBMITTED: bench/ck_sph_ab/{off,on} (restart from cs_mhd_prod3 rot 283, +1 rot,
  apudev 14.5 min, rt_desum on); binary must be REBUILT after the beam fix (README sect. 5).
* M1 for ck scoped (user asked): multigroup implicit_x1 with fixed Eddington = the two-stream at the same
  cost, no gain; transverse exchange ~100x weaker than radial on the dhj grid; beams are M1's weak spot.
  Long-term = short characteristics per band/g-point (heating = kappa (J - B) directly, beam + tensor).
* cs_mhd_prod3/rst holds 568 restart files x 82 MB = 46 GB (cleanup candidate, ask the user).
09-22 early: beam correction committed 2db7095c (thermal gates: star-off L(r) flat to 1.2e-3 vs 2.47x growth
with the old form; day column mu0 0.92 closes to 0.96, mu0 0.38 only 1.61 = NOT closed, needs longer
relaxation; T(p): photosphere night -144 K / day +20 K; p < 1e-3 bar old form collapses onto pfloor).
Option B (problem/ck_beam_sph) running with the same agent. USER: the implicit ck column solve (queue
item 3) MUST ALSO SUPPORT option B -> in the brief: the beam deposit is a T-independent source at frozen
opacity, so it enters the Newton RESIDUAL only (no Jacobian term); if opacities are refreshed inside the
Newton loop the chord-matrix tau_ray must be recomputed per iteration (state the cost); gates must include
ck_spherical x ck_beam_sph x implicit on/off combinations.
09-22: option B committed 13a7cf30 (problem/ck_beam_sph, default off, on-the-fly chords, CPU 1.12x). EXACT
PHOTON BUDGET (tests_ck_sph/photon_budget.py): for a spherically symmetric atmosphere P_new/P_exact =
1.000, the OLD tau/mu0 beam absorbs only 0.77-0.86 of the exact power (loses 14-24 %); effective absorbing
radius r_eff = 1.398e10 = r(tau_vert=1) + 4.3 H (H 3.7e8; (r_eff/r_abs)^2 = 1.27 = transit-radius effect).
Beam power reaching the ck cut is DROPPED (0.003 % here). Twilight heating down to mu0 ~ -0.5. cs corner
slots fixed ec0a9724. RUNNING 09-22: implicit ck column solve (tests_ck_implicit/), cs per-region gate +
WB/reconstruction measurement (tests_cs_regions/), L_z drift cs vs sp (tests_cs_angmom/). Vertex REAL fill
build queued after the region gate. GPU items (WB-on cs arm with wb_cache_every=1; ppmx/wenoz arm; ck A/B
smoke staged in bench/ck_sph_ab) need the user's go.
09-22 ANGULAR MOMENTUM (tests_cs_angmom/): inertial L_z from the full-3-D float dumps; cs w0 velocities are
CONTRAVARIANT gnomonic components (projection onto the zonal unit vector is exact for a linear moment). cs
cs_mhd_prod3: dL_z/L_z0 -2.4e-2 over 283 rot, plateau after rot 150, drift 100-end -8.3e-5/rot; sp arms
-2.2e-2 (230 rot) / -3.4e-2 (174 rot). => cs is NOT worse than sp; the old static 6.9e-4/rot estimate is
10-30x too pessimistic. Upper bounds only: bottom sponge (p > 50 bar, dominant), top sponge, early Rayleigh
drag, open radial BCs, outer Maxwell torque cannot be separated without an L_z history + surface torques.
Side: sp_mhd_prod3 (polar HLLE) zonal winds flip sign between dumps at mid/top shells; sp_mhd_nopole is
smooth and tracks cs. -> cs vertex list item 5 CLOSED (no action). HORIZONTAL CONDUCTION for dhj: physically
negligible (chi ~3e4 cm2/s at 10-100 bar -> 1e13 s per cell); rad_cap_ang exists for stability only and
acts as smoothing in evacuated cells; proposed apudev arm with angular conduction off (asked the user).
09-22 01:00 GPU RUNS SUBMITTED (user approved): one HEAD binary 916dc953 (md5 d34468ab..., clean-snapshot build),
bench/ck_sph_ab/{off 11931202, sph 11931203, sphbeam 11931204, prodbin 11931220 = ORIGINAL production
binary as control (mine)}, bench/cs_recon_ab/{plm_ctl,ppmx,wenoz} 11931205-7 (FROM SCRATCH, nghost 3: ppmx/
wenoz cannot restart from an nghost-2 file), bench/cs_wb_cache/{every1,every10} 11931208-9 on apu (10 h,
16 rot). bench/cs_noang/noang STAGED, NOT submitted: the agent's sbatch was DENIED by the permission layer;
do not submit it for the agent -- the user must run it or tell me explicitly.
EARLY FINDING (3 min into the runs): the `off` control (old ck form, HEAD binary) collapses dt 16.7 -> ~1 s
within 600 cycles with 109112 cells clipped by rt_de_max, while `sph` (ck_spherical = true, one-line diff)
stays at dt 12.5-12.9 like production. Either the old form is fragile to the cold EOS-cache restart or the
HEAD old path regressed vs the production binary -> prodbin arm decides. FACTS from the GPU agent: on cs
`reconstruct` governs ONLY the angular sweeps, the RADIAL sweep is hard-wired GridPiecewiseLinearX1
(mhd_fluxes.cpp:132,193-204); <mhd>/fofc was never on in production (dhj.log fofc column == 0; use
eos_dfloor/efloor/tfloor/eos_fail); cs_mhd_prod3_wb ran WITHOUT cs_wellbalanced_src (NOTES.md wrong; key
absent, default false) and died at rot 12.27 (dt 9.9e-3 < dt_min); the taper route rad_tr_tau_lo=1e30
is an exact zero for the transverse conduction.
09-22 01:30 RETRACTION (important): "all dhj productions carried an O(1) RT energy-budget error" is FALSE.
The production binaries (18b9c563 / 27ca5b13, 09-11) applied the ck source as the PLANE-PARALLEL divergence
src = -(Ft-Fb)/dx1 (two_stream_rt.hpp:1963 there), consistent with the plane-parallel sweep. The he4 commit
1159a8f3 (09-17, merged into rt-integration 09-21) made the SHARED apply kernel area-weighted while the ck
sweep stayed plane-parallel -> the non-conservative HYBRID (budget +0.63, gap 36-40 %, top collapsing onto
pfloor) exists ONLY at HEAD with ck_spherical = false; no production ever ran it. All "OFF" numbers in
tests_ck_sph (incl. the T(p) 'old vs new' comparison and the pfloor collapse) describe that merge
REGRESSION, not production. GPU evidence: every HEAD arm with ck_spherical=false is sick (off: dt 16.7 ->
1 s in 600 cycles, 109k cells clipped; from-scratch plm_ctl / ppmx died at t ~ 2.5-2.9e4 s with T ~ 1e10 K
floor cells at r 1.31e10; every1 heading the same way), while sph / sphbeam run at the production dt
(14.2 / 11.3 vs ~20 cycles/s = 1.4x / 1.75x cost). CANCELLED as invalid: wbcache every1/every10, recon
wenoz, noang (11931207-9, 11931368). Still valid: sph, sphbeam, prodbin (11931220, original binary).
TODO: make ck_spherical=false restore the production plane-parallel divergence for the ck path (asked the
implicit-ck agent, which is editing that kernel); then redo 'production plane-parallel vs spherical' T(p)
and resubmit recon / wb-cache / noang arms on a CONSISTENT form (decide with the user: plane-parallel
like production, or ck_spherical = true).
09-22 01:45: prodbin control (original binary, same restart) runs HEALTHY at dt 20.4-20.9 s, 25 cycles/s ->
the `off` collapse IS the merge hybrid. sph runs at dt 12.5-12.9 (1.6x smaller than production: reason
unknown until dumps are compared). rad_angular switch committed 213095e4. USER: resubmit with
ck_spherical = true -> RESUBMITTED: recon plm_ctl/ppmx/wenoz 11931395-7 (apudev, from scratch, nghost 3),
noang 11931398 (apudev, restart rot 283; its control = the sph arm 11931203), wb-cache every1/every10
11931399/11931400 (apu, 10 h, 16 rot). First attempts archived in each arm's hybrid_invalid/.
09-22 02:00 cs PER-REGION GATE committed (tests_cs_regions/, ~50 min, run_regions.sh). RESULTS: vertex is
real (field loop order interior 2.0 / seam 1.9 / VERTEX 1.6, ratio growing 1.75 -> 3.06; atmosphere v_t
vertex/interior 2.5-2.9). MY RANKING WAS WRONG ON #1: cs_wellbalanced_src does NOT help the vertex (8 %
worse), wb_x1 null, wb_cache_every 1 vs 10 = 0.06 % on static tests (the stale-cache pump needs a moving
background; GPU arms 11931399/400 test that). #2 CONFIRMED: angular ppmx / wenoz cut the vertex error of
the atmosphere to 0.31 L1 / 0.15 Linf, spurious max|v|/c_s 5.6e-2 -> 1.1e-2, cost 1.28x, dt unchanged;
ppm4 is BAD (2.4x worse at the vertex on rotation); nothing fixes the field-loop vertex L1 (CT/EMF at the
3-valent corner?). Remaining vertex term: div F vs geometric source cancel 1:600 interior, 1:260 at a
vertex. cs x1 reconstruction is ALWAYS GridPLM (new candidate: make it selectable). cs_wellbalanced_src
unreadable from <hydro> (gap).
09-22 04:00 GPU A/B ANALYSED (tests_gpu_ab_0922/README.md): sph's dt 12.5 s is the OHMIC limit
cfl*dx1^2/(6 max_eta) with eta pinned at max_eta 1e13 in the r/Rp 1.23-1.24 shell, because sph's NIGHT
side is 330-550 K COLDER than the production binary at 1e-6..1e-3 bar (more neutral gas -> eta cap);
prodbin is MHD-CFL limited at 19.8 s. Floors: efloor 130x, tfloor 49x prodbin; ME -28 % in 0.47 rot (but
the invalid `off` arm shows the same first-row ME step -> part is binary/restart, not the switch).
ck_beam_sph warms the terminator +721 K at 1e-6 bar and brings the night top back to within 45 K of
prodbin, recovers ~5/8 of the floors; the two together: 2.8x (sph) / 3.4x (sphbeam) throughput penalty.
noang == sph to 2e-3 and NOT cheaper (transverse conduction already negligible). Recon arms: falling dt
= normal spin-up; NO vertex-localised spurious flow in the real atmosphere (vertex bands quieter than
interiors, 0.8x), ppmx/wenoz null on floors/KE/ME and 1.44x cost -> not worth it for dhj. UNDETERMINED:
ck_spherical alone cannot be isolated (prodbin is a different binary); needs a 15-min HEAD run with
ck_spherical=false (now the repaired plane-parallel path, ed188f69) on the same restart.
09-22 05:30 offfix ISOLATION (tests_gpu_ab_0922): (i) the SWITCH: ck_spherical alone cools the night top
+533 K at 1e-6 bar, saturates the Ohmic max_eta cap (dt 19.9 -> 12.5), doubles the ME decay; per-cycle
compute only 1.22x (the 2.8x throughput loss is the dt). (ii) the BINARY LINEAGE (offfix vs prodbin, switch
off both, 364 commits): dt and T match, BUT ME decays 6x faster per rot (-35 %/rot vs -5.5 %), dE/E
-8.9e-5 vs ~0, efloor/cycle 9.5x, per-cycle cost 1.46x; first-row ME step is a restart artefact (missing
EOS cache). UNBISECTED -> NEW LIST-C ITEM: bisect the ME-decay / efloor / cost regression between
27ca5b13 and HEAD on the rot-283 restart (candidates: he4 merge, corner-slot resample ec0a9724, restart
cache commits, hipcc DualView change b4c6a2f2, floors). Must be done BEFORE prod4 launches.
09-22 07:00 BISECTION DONE (bench/bisect_me/README.md, 8 apudev jobs): the ME-decay regression (-35 %/rot
vs -5.5 %) is commit 15e9f790 (in merge 46410519, lhlld, 09-11): HLLD_BX_ZERO_TOL 1e-4 -> 1e-8 in
hlld_mhd.hpp / hlld_uct_mhd.hpp, unconditional. The test 0.5 Bx^2 < tol * pt* is a plasma-BETA cut for the
rotational (**) states; 1e-8 switches the whole beta 1e4-1e8 dhj envelope from the collapsed 3-state to the
full 5-state HLLD, whose unresolved rotational jump dissipates B_theta/B_phi. Validated on the Balsara vortex
(beta ~1-100) where 1e-4 was unstable. Confirmation arm 46410519_tolrevert (job 11932540) pending.
SEPARATE: the 9.5x efloor rate + part of the 1.46x cost enter at 6c0df36b (merge polar-average-perf),
NOT bisected inside. DECISION NEEDED: make the tolerance an input (default?) -- the two regimes want
different values; both dhj productions ran with 1e-4.
09-22 08:00 BISECTION ROUND 2 (bench/bisect_me/README.md): (1) CONFIRMED 15e9f790 (HLLD_BX_ZERO_TOL 1e-8) is the
WHOLE ME-decay regression (revert -> slope -0.167 = pre-merge); fix in progress = <mhd>/hlld_bx_zero_tol
default 1e-4 (user). (2) The 10x efloor jump = d8f00d49 "ck and picket-fence sweeps: centre-to-centre
layers" (polar-average-perf, merge 6c0df36b), behind problem/rt_layer_legacy (legacy=true -> bitwise back).
Mechanism: the old whole-cell layer offset by half a cell smeared the tenuous TOP SLAB (cells on pfloor
1e-3 / dfloor 5e-14, T -> tfloor 200 K, B_b jumping cell to cell); centre-to-centre half layers remove that
dilution, the per-cell net source exceeds e/dt, rt_de_max clips, the semi-implicit apply overcools ->
efloor 1.3e5 -> 1.45e6 per 1000 cycles. The two regressions are independent; layers cost only 2.5 %.
NOTE: the user chose to KEEP rt_layer_legacy in the clean-up -- it is the production-era layer form.
DECISION NEEDED: (a) rt_layer_legacy = true for dhj (bitwise the production layers; but it is refused by
ck_spherical / ck_beam_sph / ck_implicit -> those need the c2c layers), or (b) fix the top-slab
overcooling under the c2c layers (floors / a dilution-consistent top closure), or (c) accept the
efloor rate. Also: command-line problem/ overrides need the key present in the athinput.
09-22 10:00 TOP SLAB RESOLVED (tests_ck_topslab): the 10x efloor rate is NOT a bug in the c2c layers -- they
are the correct emission; the OLD layer averaged the top cell's B with the ghost above the domain (0.37 B)
and kept it artificially hot. The colder (correct) slab hits pfloor = 1e-3 bar, a value calibrated for the
old layers. FIX = <hydro>/pfloor 1e-5 in the dhj input (efloor 27048 -> 0 on the reproducer, T(p) unchanged
to 5 digits, dE/E -1.7e-5 = physical). -> prod4 input: pfloor = 1e-5 (check the earlier note
dhj-floors-for-1e-6-bar which chose pfloor 1e-3 with tfloor_kelvin 50 crashing). GPU confirmation arm
still to run. ck_implicit on GPU = 4.95x -> not for production; frozen operator committed as option.
09-22 11:00 GPU PROFILE AT HEAD (bench/prof_0922/README.md, rot-283 restart, 2 GPUs, 300 cycles, rocprofv3;
ck_spherical ON, ck_beam_sph off): RT = 82 % of GPU time (was 37 % on 09-12); ONE kernel, the ck sweep
(picket_fence_two_stream_RT_pass = chain + pre, two_stream_rt.hpp) = 70.6 %; conduction weights 3.2 %,
fc pack/send 2.6 %, rt_apply 2.3 %, radial implicit conduction 1.9 %; everything else (MHD, EOS, cs terms)
in the remaining ~18 %. Host/launch/MPI overhead ~21 %. ck_implicit 4.95x. => the SWEEP is the whole
game; with the spherical form the per-cycle sweep cost roughly doubled vs 09-12 (face mixing, two-pass
probe). Levers, in order: (1) occupancy/latency of the sweep kernel (VALU 5 % on 09-12: more resident
waves -- split chains further, FP32 recurrence RT_FP32, larger launches), (2) the spherical-form extra
work (profile the c2c/face-mixing arithmetic; the 2-pass probe), (3) fewer chains (physics: ck_nquad,
bands). A 2x on the sweep = 1.55x on the whole run.
09-22 12:00 SWEEP COST (tests_ck_sweep_cost): ck_spherical = 1.32x on the loop, 1.74x on the sweep kernel,
because the u/d sweeps need 2 PROBE passes (4 column passes per chain) and every pass rebuilt the layer
operator twice per cell (4 kappa lookups, 8 expm1, 24 divides per (cell,chain) where 1/1/3 suffice).
ck_sweep_cache (default 2 under spherical, 1 otherwise; bitwise) recovers only 4-7 %: the extra passes
are the cost. Occupancy: 1056 waves/GPU = 1.16 waves/SIMD, vgpr 128 in both forms (not register-limited);
only RT_NB (chain block, not bitwise) could raise it. IN FLIGHT (Opus): probe-free forms ck_sweep_form =
sd (S,D decoupled sweep) | tm (transfer matrix / adding-doubling), both exact, 2 passes, GPU-timed.
NEW UNEXPLAINED: the sweep-cost agent's base binary (HEAD 446b7b30 + uncommitted Jacobian work) runs
29.9 s where prof_0922 (a544a761) ran 23.0 s on identical cycles = 1.30x regression -> bisect between
a544a761 and 446b7b30 (candidates: 90c33ce5 frozen-op/topslab, 0aaf7f39, or the uncommitted phase-4
code) BEFORE committing phase 4.
09-22 13:00 OHMIC CAP (tests_gpu_ab_0922 "Ohmic cap in sphbeam"): sphbeam's dt 12.85 s is held by the SAME
night-side shell as sph (i=42, r/Rp 1.241, T 1406-1461 K, p 7 mbar, x_e ~1e-9), 100 % night side. eta
reaches max_eta below T_cap ~ 1450 K (alkali ionisation); the binding cells are only 6-14 K below it,
prodbin is 300 K clear. The beam warmed the TOP (1e-6..1e-3 bar) but this shell at 7 mbar is deeper
than the twilight heating reaches. LEVERS: max_eta 1e13 -> 5e12 returns both arms to the CFL limit
(sph 17.4, sphbeam 18.0 s; 9 % under prodbin) at no cost to the old form -- a dhj INPUT choice, needs
the user's physics OK (it caps the resistivity lower in the neutral night side); raising max_eta is
worse. -> DECISION FOR THE USER: max_eta = 5e12 in prod4? (memory: dhj-ideal-xe-floor-relaxation
found max_eta 1e13 beat 1e14 by 1.54x; 5e12 was not tested then.)
09-22 13:30 USER DECISIONS (prod4 draft, committed): max_eta = 5e12 (Ohmic cap), pfloor = 1e-5 (correct
top-cell emission; GPU arm bench/pfloor_ab/p1e5 confirms or not), ck_sweep_form = tm (exact probe-free
spherical sweep, 14.35 cycles/s). The CODE default of ck_sweep_form stays the 4-pass form until the
tm default flip is done (not bitwise for spherical runs; do it as its own gated commit). prod4 launch
still needs: MPI gate verdict, the 1.30x cost bisection, the pfloor GPU arm, and the user's go.
09-22 19:10 COST REGRESSION BISECTED (bench/bisect_cost/README.md, 1 apudev job): the 1.30x enters ENTIRELY at
90c33ce5 (ck_impl_frozen_op, phase 3, default OFF): its runtime bools put both branches of every
frozen-op test inside the dominant rt_chain_ck kernel (~8 sites per cell,chain) + 5 dummy Views: the
production instantiation grew 24124 -> 28391 instructions (+18 %), scratch ops +28 %; at 1.16 waves/SIMD
that serial/spill work IS the 1.30x. LESSON (standing): a default-off switch inside the sweep kernel must
be a COMPILE-TIME template tag (as ck_sweep_cache's CCH, 0.3 % when off), never a runtime bool.
FIX QUEUED: template-tag ck_impl_frozen_op (or drop it: it is a small LOSS on GPU anyway). Same audit needed
for the phase-4 bools (ck_impl_reuse_jac / ck_impl_seed) and ck_sweep_form's switch if it is runtime.
09-22 19:30 pfloor 1e-5 barye CONFIRMED on GPU (bench/pfloor_ab, tests_ck_topslab §5): efloor 7.4x down,
0 active cells on the pressure floor (was 1), tfloor safe (0 active cells at 200 K; the counter counts
ghosts), dt identical, energies identical, night T(p) within 15 K. prod4 keeps pfloor = 1e-5 (<mhd>,
barye). tm default flipped 0172f3bc. Kernel-tag fix of the 1.30x regression in flight.
