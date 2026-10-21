# gsdrag-1021: one gas source for all prod runs (conservative2 + enthalpy-flux cap on the gas charge)

Branch gsdrag-1021 off pgrst-merge-1021 077dda23 (worktree /raven/ptmp/jinma/gsdrag_1021/wt). Coordinator merges.
Builds: src copy src/s + build/{acc,he} (relocated copies of pgrst_1021 build/*_4146127b), jobs/incr.sh.

## FINDING that changed the criterion
E/(rho c^2) is NOT >> 1 anywhere in the prod setups (gam.py on outputs): accretor static t 0.002 max 4.9e-3 (floor
rho 1e-8); RSG E clip N32 max 7.0e-5; HG global clip N32 max 7.5e-7; hg145 prod (t 8.9e6) max 6.4e-7.
So a Gam = E/(rho c chat) switch would never act. Kick size matches instead dv ~ |dA|/(rho c chat) with
|dA| ~ c E (dv ~ Gam c = 5e-3 * 3e5 = 1500 km/s vs max 1175 observed): A^old = F_stored - (F0 face mean) is
NOT an enthalpy flux where the stored cell F was changed outside the face system (realisability clip of the
write-back -- which implicit_cs2_clip = rad keeps out of the write-back's own charge but which the NEXT solve's
A^old still carries --, pgen-held E/F in donor/stream light). conservative2 then hands that ~cE to the gas.

## DESIGN (opt-in key <rad_m1>/implicit_cs2_dacap = K, default 0 = off -> bitwise; K = 2 = the strict bound below)
- A physical enthalpy flux A = (v + v.D) E has |A| <= 2 |v| E (D <= 1), so a legitimate change is bounded by
  |dA| <= 2 (|v_old| + |w|) max(E^old, E') ~ 4 |w| E. Per cell, per Picard pass (ImplicitCs2DeltaA, the dA of the pass):
  s = min(1, K * 2 |w| max(E^old, E') / |dA|), |w| = max(|v|, |v + dv_vimp|), |dA| the vector norm.
- The cell's dA is used as s dA EVERYWHERE it enters: (i) stored dA (iw_cs2+3..5) -> the F0^n shift of each face
  equation (ImplicitCs2Pass, face mean of s dA: face-consistent, each face carries each neighbour's own s),
  (ii) the vimp neighbour share (Cs2Vimp nb), (iii) the cell's own gas share om_bar -> s om_bar (iw_cs2+6..8:
  write-back charges s om_bar dA; Cs2Vimp ob and cown x s; s itself in a new slot iw_cs2+9, allocated only with
  the key). s = 1 = conservative2 exactly; s = 0 = the cs2_hll = drag rule (gas pays drag only, as matched).
- Booking: ENERGY exact as before (the kinetic-energy change kew(dm) of whatever is charged leaves E: e_gas +
  (c/chat) E per cell unchanged). MOMENTUM (and L_z): conservative except sum_c (1 - s_c) om_bar_c dA_c/(c chat) V_c,
  which stays with the radiation (as the cs2_hll drag rule / matched) -- only in capped cells, i.e. where the
  stored F is not F0 + enthalpy flux.
- Under implicit_cs2_hll = inertia + key: the om_bar slot := s.

## vimp (from PGRST worker, 10-21) -- the hesdirk2 stream failures; proposal + test plan, for agreement
Data (CPU dbg_t2_admiss, bin pgrst 2de67c53, v12lzc + time2_k1_reset + _grow=1, SS from g12k2/ST rst; runs in
/raven/ptmp/jinma/pgrst_1021/runs/admk3, admk5_phi45, admk5_onf):
- After k1_reset(+grow) every remaining stage-1 failure is E_solved (S2) <= 0 in 1e-8 floor / back-sink cells right
  beside or in front of the nozzle (i 323-327, r 14.43-14.74), old vector fine (0.02..1.3), cell residual ~0. It follows
  the nozzle (stream moved to phi 45 deg: same, so not a block-edge effect); implicit_opac_newton=false: same.
- Worst cell (0,10,31,323) r 14.43, in front of the nozzle inner face, rho 1e-8, f_rr 0.48 (NOT free-streaming/clip-limited by
  that number), v_r -1.09: row TR = base 32.9 + face_F0 8.3 + opac_newton -21 + vimp_JRHS **-6.70e3** -> S2 = -16.7.
  Its r+ neighbour is the nozzle cell, E = 496 (pgen-held E = a T^4 of the 73 kK core, F = 4/3 E v_stream).
- vimp_JRHS (ImplicitVimpBuild, rad_m1_implicit.cpp ~11700-11785) = -sum_faces sg*(ph*dvk + ps - ph*pe) with
  ph = (A_f/V) dt cr E_f (1 + D_f) and dvk = face mean of the PREDICTED kick DV of the two cells. At the nozzle face
  E_f ~ the nozzle's 496 while dvk carries the floor cell's large predicted dv: the work term charges the floor cell
  E_f*dv ~ 7e3 (its own E ~ 1).
- Bisection (SS, 1 arm per key): implicit_vimp = false: NOT ADMISSIBLE 106 -> 62; gas_src = matched: 106 -> 99
  (fast cells 73k -> 9k). So both dA (conservative2) and the vimp work term feed the failures; under BE vimp is off
  (resolves to 0) and the stream is clean (0 redos).
Proposal (one criterion with your dacap, no Gam): where DV of a cell comes from the cs2 prediction (Cs2Vimp), your s
  already scales the neighbour share and cown -> check first whether implicit_cs2_dacap alone removes these failures
  (the nozzle cell's A^old = F_held - F0 mean is exactly your "pgen-held F" case, and its dA is shared to the floor
  cell's DV). If not (the matched arm still fails 99/105), the vimp work term needs the same physical bound:
  a real enthalpy-flux change at a face is |dA_f| = |E_f (1 + D_f) dv_f| <= 2 E_f |w_f| with |w_f| = max(|v|,|v+dv|)
  of the face's cells -- the floor cell's predicted dv is NOT bounded by the nozzle's velocity. Treatment: per face,
  dvk -> dvk * min(1, K * 2 max(|v_L|,|v_R|) / |dvk|)?? (too strong for physical kicks) -- better: weight the face E by
  the cell that owns the kick: ph*dvk -> nus cr (1+D_f) * 0.5*(E_L dv_L + E_R dv_R) (upwind-consistent with the
  enthalpy flux E v: each cell's kick carries its OWN E), applied only under the dacap key (same switch, bitwise off).
  I (PGRST) will TEST, not implement on your branch: as soon as your bin with implicit_cs2_dacap is built, I run the SS
  stream arms {k1r+grow, k1r+grow+dacap=2} and report NOT ADMISSIBLE + the row terms; if failures remain I implement the
  E-weighted face work on a branch off your tip (ImplicitVimpBuild only; you do not touch that function) -- reply here
  if you prefer otherwise.
- PGRST test (job 31353082, your bin g1b, SS from g12k2/ST rst, runs /raven/ptmp/jinma/pgrst_1021/runs/bis2):
  k1r+grow: NOT ADMISSIBLE 106/105 cyc, fast 72.8k | k1r+grow+dacap=2: 103, fast 13.0k (max 896) | k1r+grow+dacap=1: 103, 13.1k |
  dacap=2 without k1r: 129 + 4 dt_retry, 13.6k.  => dacap fixes most fast gas (13k vs v11 8.2k = 1.6x) but NOT the stage-1 failures.
- PGRST now implements the vimp face-work fix itself on branch pgrst-vimp-1021 (off 077dda23, ImplicitVimpBuild only):
  opt-in <rad_m1>/implicit_vimp_ewf (each cell's kick with its own E), commit b3f989f8; will test with/without your dacap.

## REVISIONS
- g1c: + implicit_cs2_dacap_vmin (default 1e-4, read only with the cap): |w| >= vmin c, so (near-)static cells keep
  their truncation-level dA (RSG with key -> bitwise to key-off).
- g1d: + write-back cap: the conservative excess over the matched drag (dt k_f F0'_f/c, half per face) is capped at
  2K |v| max(E^old,E')/(c chat) (same bound; rest stays with the radiation).  Effect on the accretor: none measurable.
- g1e: the dA cap is a separate kernel AFTER the cubed-sphere basis change (g1b-g1d judged the cs dA before it,
  spuriously capping on cs meshes = RSG/HG). sp/Cartesian arithmetic unchanged.

## RESULTS (A100, v12lzc, nx3=64, 1 GPU/arm; fast = |v| > 50 km/s; SS t 0.0034; v11 SS 8.2k)
| arm | ST fast | SS fast | SS >200 | SS max | notes |
|---|---|---|---|---|---|
| BE cs2 no key (g12be) | 26.8k | 73.3k | 1982 | 1028 | |
| BE cs2 + light off D/S/DS | 19.1k/26.8k/19.1k | 52.8k/62.0k/51.0k | | 915/893/647 | light is NOT the driver |
| BE matched | 6.6k | 8.2k | 298 | 974 | |
| BE K2 (g1b / g1d) | 18.5k | 18.9k / 18.7k | 2728 | 828 | floor r 9.6-25 cleared; residual r 9.1-9.5 |
| BE K1 / K0.1 | 18.6k / 18.7k | 18.9k / 18.9k | | | saturated |
| BE conservative (plain) | 196k | 215k | | | |
| hesdirk2 cs2 no key | 27.1k | 67.2k | | 1141 | dtr 11, notadm 192, 184 cyc |
| hesdirk2 K2 (g1b) | 18.9k | 19.3k | 2812 | 814 | dtr 1, notadm 111, 110 cyc; repeat bitwise |
Residual (K2): >200 km/s cells are equatorial (theta 80-100, peak 90), all phi, r 9.13-9.53, rho ~3e-7, |v_phi| ~300
(v_r, v_theta small): an azimuthal spin-up of the envelope-top free-streaming (|F| = cE) gas; matched: 68 such cells.
Budgets (BE K2 SS): rM 3e-10, max|Jres| 1e-6, rE/|Ebnd| -1.05e-3 (no key -2.1e-3, matched -1.0e-3).
| BE K2 (g1e) | 18.5k | 18.7k | 2744 | | g1e == g1d bitwise on sp |
| BE K0.01 vmin 0 (+/- cs2_vimp off) | 18.7k | 18.9k | 2850 | | gas charge ~ matched drag, STILL 18.9k |
| hesdirk2 K2 (g1d) | 18.8k | 19.6k | 2836 | | dtr 5, notadm 135, SS 131 cyc / 157 s |
=> NO-GO vs <= 1.5x v11 (12.3k). The key removes the floor kicks r 9.4-24 (light/clip A^old), but the residual
equatorial v_phi spin-up at r 9.1-9.5 is NOT in the gas momentum charge: with K 0.01 the charge is ~ the matched drag
and the count is unchanged. The residual comes from the radiation side of conservative2 (stored cell F / F0 handling,
lab-F face equation) vs matched's (rinr a W/c^2 fixed point, F re-forming) -- next bisection target.

## GATES (all on the final code 822d1e61 = g1e, unless noted)
- default bitwise: RSG E clip N32, HG global clip N32 (rst/hst/bins vs fxb.31348879), hg145 (md5 vs mg.ncfix), accretor
  v12lzc static P(077dda23)==N: all BITWISE.
- key on (dacap 2): RSG BITWISE (rst+hst; bins differ only by the param dump); hg145 BITWISE (rst+hst); HG global:
  <= 9e-7 rel over cycles 1-3, then one extra time2_reject retry at cycle 3 (resid 9.27e-4 both arms, borderline) ->
  path diverges (5 cycles end t 323.5 vs 378.5).
- accretor hesdirk2 + K2 (g1d == g1e on sp): 2 runs BITWISE, restart BITWISE, forced retry BITWISE (lzgb.31354588).
- IDORT (CPU g1d, PROD decks of idort1020 p arms, N vs K=2): eq1, eq2, eq3t (moving), rs_m1.2/m2/m3 n512, rw_p10t1 /
  rw_p10t100 X_64: outputs data-identical (bins differ only in the param dump); cloud r12g2 max 1.9e-5 rel (qualitative
  row) -> all rows stay PASS.

## BISECTION of the residual (coordinator 10-21) -> FOUND: realisable coupling
- matched with implicit_realisable_coupling = false: ST 18.7k (>200: 3894), SS 18.9k (2826) = the cs2+K level; matched
  (true): 6.6k / 8.2k. No other radiation-side key moved it: implicit_cs2_dbg = 1 (no face split, rad_inertia after
  the solve) K0.01 18.8k / K2 18.0k; implicit_cs2_vimp / implicit_vimp / implicit_rad_inertia off, cs2_hll = inertia:
  unchanged. => conservative2 hands the gas the full drag of the UNREALISABLE iterate flux (|F_iter| > cE); matched
  scales its deposit by min(1, cE/|F_iter|).
- FIX g1f (24810dad): under the dacap key the write-back's drag reference is the realisable matched drag (x min(1,
  cE/|F_iter|) when implicit_realisable_coupling) + the inertia excess capped at 2K|v|E/(c chat).
| arm (g1f) | ST fast (>200) | SS fast (>200) | SS max | notes |
|---|---|---|---|---|
| BE K2 | 7.4k (316) | 9.07k (562) | | 1.11x v11; rE/|Ebnd| -1.0e-3, max|Jres| 4e-7, rM 2e-10; SS 68 s |
| BE K1 | 7.0k (164) | 8.6k (398) | | |
| BE K0.01 | 6.6k (68) | 8.17k (300) | | = matched |
| hesdirk2 K2 | 8.1k (518) | 10.5k (872) | | 1.28x v11, but dtr 15, notadm 199, dt min 1.1e-6, 259 s |
- gates g1f: default bitwise (RSG/HG/hg145/accretor); key on: RSG BITWISE, hg145 BITWISE; HG global differs
  (radial momentum 1.6e-3 rel, KE 1e-3, already ~2e-6 in row 0) -- the realisable drag acts in its free-streaming
  atmosphere (no retry flip any more). accretor hesdirk2 K2: 2 runs / restart / forced retry BITWISE.
