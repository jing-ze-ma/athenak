# LATEST_KEYS (INTEG, 2026-11-21: v20.3): production keys for HG / RSG / R136 runs on clean-cc-stag2-latest

**Code:** LATEST v20.3 = **e688cb7bc** (INTEG7, awaiting the USER's push decision; push target origin/clean-stagger2-latest, fast-forward from v20.2 a0c19e459). v20.3 = v20.2 + four commits + three key changes (HG/RSG):
- **L_z fix, key ON for HG/RSG (USER 2026-11-20 ~14:55 PST): `rad_m1/clean_cs_lz_tel=true`** with the LZTEL3 form 74c8a8bb1 (here e688cb7bc). Cubed sphere + staggered only (FATAL otherwise); default false = bitwise to v20.2. The radiation's L_z is that of its face momenta, and the gas takes the face-form absorption torque. The net lateral-stress torque defect X = sum xs (a truncation defect of the cs lateral face stress) goes back to the gas each solve as the least-kinetic-energy solid-body spin dOmega = X/I (I = sum rho V R^2). L_z of gas + radiation + ledger is then exact; the ledger holds only the radial boundary torque. **Caveat:** this is a global conservation fixer, not a local source fix. X moves between radii: on RSG X ~ -2e46 per solve, always negative, which is a deep rigid counter-spin of ~2.7 m/s at R_ph per 1000 cycles (~300 m/s per 1e5 cycles) balanced by the outer layers. The well-balanced lateral-stress source fix (LZ4) will replace it. Compare Lz_tot only within one setting: the radiation L_z changes definition at the switch. Evidence: CPU Lz closure 1e-11 per cycle, restart bitwise; Caltech 600 cycles HG dt 107.7 vs 96.4, CW 0 vs 10, Lz 3e-12; RSG CW 10 vs 13, s/cyc -0.5 %, L_top and Mdot equal to off; G_rsg 1000 cycles Lz -1.5e-11. Cost: two global sums per write-back, no measurable s/cyc change.
- **EMA named for HG/RSG (USER 2026-11-20): `rad_m1/vet_rotate_avg_ema=128` + `rad_m1/vet_rotate_avg_ema_change_form=step`.** DOC CORRECTION: the K1 args named ema `0`. Since v20 the keys files did not name it, because UNIFORM wrongly listed it as inert (it is read under clean, rad_m1.cpp ~3370). So a restart took its file's value: hegiant.00184 / rsg.01510 carry 128 + step (EMA ring active), while the CALT2 gate restarts rsg.01295 / hegiant.00139 ran the K window (log: "mean over the last 4/8 rotated sets"). It is now named so every restart runs one form. R136 is unchanged (r136.00097 carries no EMA: K4 window). `vet_rotate_avg_ema_wb` stays OFF (USER: a Caltech wb on/off pair decides).
- **Fresh decks:** the three keys are FATAL on a fresh run's command line (as floor_pay is). Name them in `<rad_m1>`. Restarts add them with a NOTE. On sp fixed-ray decks, `vet_rotate_avg_ema=0` + `clean_cs_lz_tel=false` parse and have no effect (BSG smoke hst bitwise to the deck without them).
- **Code pieces:** EMAWB 9129916f2 (`vet_rotate_avg_ema_wb`, read only when named, default false = bitwise; FATAL without ema >= 1) + a6845503c (EMA fill0: outside the steady start the EMA starts as the current set alone; production restarts read the ring, bitwise); FTIEFIX f3239f458 (= f43326d06; CleanStagFaceTie wrote stage-slope slot CLF_KD past the face arrays without clean_stag_fkin: cs uniform blowup without consistency; production with vet_gd has fkin on = no-op; a bounds-checked sweep finds no other OOB).
- **Order ladder note (v20.2+):** the theta-wall fix honours a named vet_gd_theta_wall, so the LAG/ROTWB deck-D order ladder needs `vet_gd_theta_wall=mirror` (its mk.py names vacuum, which v20.1 overrode to mirror). With vacuum walls the EMA path keeps an O(dt) term (EMAWB side finding; sp only).
- **New bench rows:** rhse_rotwb (rhse wedge at rest, rotating rays + EMA 128 + ema_wb, 400 cycles, max v/c_s <= 1e-7) and rhse_rotwb_off (info, the 1/N pattern drift 5e-4).
- **Not in v20.3 (pending):** RSGDIE clean_retry_budget, the THINF / BEAMHG beam fix.
- **Gate (e688cb7bc):** CPU k20 bench i7v2030: 420/420 bins + hst bitwise to v20.2 (i6v2020); GATES FAIL 8 (known), CW 2 (= reference, forced decks); new rows rhse_rotwb 6.42e-9 PASS, rhse_rotwb_off 5.2e-4 (info). Stars plain == poison == k20ref. pzgate keys off 10/10 PASS. pzgate with the v20.3 keys (EMA 128 + step + clean_cs_lz_tel, live ring, xfade, rotation every ~2 cycles; decks integ/v203/decks) 10/10 PASS (poison 1-3, ranks 4 = 12, restart 6 -> 12; 0 CW). IDORT suite 40 runs (I7D head defaults, I7K k20 keys) data-bitwise to SUITE2 SD/SK, 0 CW, rc 0. Fresh sp BSG smoke: ema=0 + lz_tel=false named == unnamed bitwise, 0 FATAL/CW. Raven A100 gpudev R136 rt bitwise, hst == V202. Caltech H200 (1 node, 4 cycles, rsg.01295 / hegiant.00139): V203_rsg/hg (v20.2 keys) rt 5/5, hst == V202F bitwise; V203L_rsg/hg (v20.3 keys) rt 5/5, 0 CW, EMA active from a K-window restart; V203LP_rsg poison == V203L bitwise. s/cycle H200: RSG 2.68 (v20.2 keys 2.70), HG 4.35 vs 3.92 with v20.2 keys (2 steady cycles, the window -> EMA switch; GPU mem 102 vs 117 GB) -- watch on the first production segment. Binaries: Orion speedhead/athena_e688cb7bc md5 a21a2c46, speedhead2_he c5916c42; Raven integ3x_cu 28b21f49; Caltech B_v203 7df2b648.

### v20.2 block (superseded by v20.3 above)
**Code:** LATEST v20.2 = **f4cefb181** (INTEG6, awaiting the USER's push decision; push target origin/clean-stagger2-latest, fast-forward from v20.1 25ea76366). v20.2 = v20.1 + two code fixes + one key change:
- **Keys change vs v20.1 (USER 2026-11-20 ~13:20 PST): `hydro/floor_pay=rad_column` for ALL stars (HG, RSG, R136).** The hydro C2P EOS floor energy is paid from the radiation E (ELEAK; closes the rule-1 energy residual: RSG AB9 per-cycle residual ~1e-8 L_in dt vs up to 5e-7 off; long HG LF: CW 3 vs 20, wall -10 %). Result-changing where EOS floors fire. Added to uniform/keys/common_scheme.txt; applied by apply_v202.sh. On a restart the key is added from the command line; a FRESH start must name it in the deck's <hydro> block (command-line keys absent from the input file are FATAL on fresh starts).
- **Theta wall honoured** ca0d8fd7b (GDCHECK): CleanInit no longer overrides a named rad_m1/vet_gd_theta_wall on spherical polar; the mesh-following default (reflecting wall -> mirror, periodic -> wrap) now applies only when the key is unnamed (rule 12, no silent override). Bitwise for unnamed keys and the production decks.
- **rhse history maxima** f4cefb181: vcs_max (and vr/vlat max/deep) are MPI_MAX over ranks; the history writer summed them (4-rank value = 4 x). Diagnostic only, layout invariant.
- **Not in v20.2 (pending the USER):** BEAMFIX/BEAMFIX2 beam fix (Caltech gate running); LZTEL3 L_z fix 74c8a8bb1 (G_on_hg pending; counter-spin caveat); CONV3 time2_k1_xlim stays off.
- **Gate (f4cefb181):** CPU k20 bench i6v2020: 419/420 bins + hst bitwise to v20.1 (Si5v201d0); the one change is rhse_cs24 (4 ranks) user.hst vcs_max = old / 4.0 (the fix); theta-wall fix bitwise over the whole bench; GATES FAIL 8 (known), CW = reference. Stars plain == poison == k20ref; floor_pay on (in deck) == off bitwise on the fresh star decks (no EOS floors in 2 cycles), 0 CW. pzgate 10/10 (poison 1-3, ranks 4 = 12, restart 6 -> 12). Raven A100 gpudev R136 rt: keys v20.1 bitwise, hst == V201D; floor_pay on rt bitwise, 0 CW, hst identical to off (R136 EOS floors create ~1e16 erg in 50 cycles: inert). Caltech H200 (CALT2, B_v202, 4 cycles, 0 CW): keys v20.1: RSG and HG rt 5/5, hst == V201D bitwise, poison == plain; floor_pay on: RSG and HG rt 5/5, poison == plain; energy residual per cycle (dEtot - (L_in - L_top + P_esrc - P_sponge) dt)/(L_in dt): RSG off -8e-9, 9.4e-8, 3.9e-7 -> on -8e-9, 1.1e-8, -1.7e-8; HG -2.6e-8 both (= the boundary flux E_hflx; no floors in 4 HG cycles); L_top/Mdot equal to 6 digits over 4 cycles; long HG (LF, v20.1e): CW 3 vs 20, wall -10 %. H200 s/cycle: RSG 2.73 (2.72 on), HG 3.89 (3.92 on).

### v20.1 block (superseded by v20.2 above)
**Code:** LATEST v20.1 = **4ecbae182** (code: merge of the v20 cut 924fb24c6 into integ-v201 3efbbbd17, src tree == 3efbbbd17; INTEG5, awaiting the USER's push decision; push target origin/clean-stagger2-latest, fast-forward from v20 924fb24c6 and v19 13ea625fe). v20.1 = v20 code (6d3da47b2) + five pieces, every key default off, so the keys files are unchanged from v20:
- **KERN kern-fw** 3ee23a639 b9c0d1739 (= 12730bdd2 7ee3e80dd): the clean face weights computed once per face (CleanFaceSource takes the flux's phi), and the clean_halo_direct pack skips entries deeper than the exchange depth. Bitwise; Caltech H200 RSG -3.4 %, HG -5.4 % per cycle (2 swapped rounds).
- **hst layout invariance** 3f4369df7 (= a3b72fdcc): with mesh/layout_invariant the hydro / user HistSums per-MeshBlock partials are summed exactly (ReproSum) and merged over the ranks. The history no longer depends on the rank count (CPU 4 vs 12 ranks and H200 LI4 vs LI8: hst byte-identical); state unchanged; hst vs v20 differs at round-off (<= 1.4e-14).
- **GHOSTGEO** d70782673 (= b596d5398): CleanStagMemUp takes the donor's outflow over the owner's volume (CLG_F + 3), not a seam ghost's own volume, so both copies of a cubed-sphere seam face take one decision (symmetry / layout). Bitwise off seams; the k20 bench is bitwise; HG cs evolution changes where a seam decision flips.
- **CONV3 time2_k1_xlim** 774262b56 (= 0e9ca76bb): rad_m1/time2_k1_xlim = th (default 0 = off, read only when named) scales the exchange pair of the hesdirk2 stage increment per cell so the stage start keeps gas eint >= th eint and E >= th E. Exact per cell (E + gas), no new state (restart bitwise).
- **VGD_CS_LB default 8 on CUDA** 3efbbbd17 (was 6; HIP keeps 2): the vet_gd cs shell kernel's launch bound. Bitwise; H200 RSG -1.1 %, HG -0.8 % (2 swapped rounds); env VGD_CS_LB still overrides.
- **Key-on / default recommendations pending the USER (all OFF / unchanged as cut):** (1) hydro/floor_pay=rad_column (ELEAK; closes the RSG energy budget to noise, 0 CW, no slowdown); (2) rad_m1/clean_cs_lz_tel=true (LZTEL2 fix; RSG 300 cycles 0 CW, L_z drift 8e-12, HG dt unchanged); (3) code default rad_m1/time2_face_fct=true for scheme = clean (production names it already); (4) rad_m1/time2_k1_stiff on (~free on H200, -16 % HG passes); (5) rad_m1/time2_k1_xlim=0.1 (HG prod kappa 1, 20 cycles: eint<0 -> 0, stage-1 E<0 1.7k-14k -> ~400, 0 CW, HG -5.7 % wall / -18 % passes, RSG -1.6 % / -7 %; conv/CONV3.md).
- **Gate (v20.1 774262b56 + 3efbbbd17):** CPU k20 bench: bins + hst 420/420 identical to the hst-fix head a3b72fdcc (kern-fw, GHOSTGEO and xlim-off bitwise on the bench; the hst fix vs Sk20r2: bins identical, rhse_cs6/24 hydro hst round-off), GATES FAIL 8 (known), CW 2 (= reference); stars plain == poison == k20ref; pzgate 10/10 (poison 1-3, ranks 4 = 12, restart 6 -> 12). Caltech H200: V201D RSG rt 5/5 + poison == plain, HG rt (romio) 5/5, 0 CW; LI4d vs LI8d hst byte-identical + bins identical (layout invariant incl. hst); kern-fw alone (V20NK) hst == v20n bitwise; GHOSTGEO alone (V201B) == V201D: production cs changes RSG from the first row (hydro 2.2e-12, user 2.4e-8) and HG from row 1 (1.1e-11 / 2.2e-6) = seam decisions that now agree on both panels; xlim off inert; V201E (LB 8) hst == V201D. Speed vs v20 (H200, per cycle): RSG 2.85 -> 2.73 s (-4.2 %), HG 4.18 -> 3.93 s (-6.0 %). Raven A100 gpudev R136 rt bitwise 0 CW, hst == v20.1 hst-fix head (V201); LB 8 (V201E) hst == V201D, wall unchanged (47.0 vs 46.9 s for 23 cycles).

### v20 block (superseded by v20.1 above)
**Code:** LATEST v20 = **2efe38c34** (code: merge of v19 13ea625fe into integ-trial-v20n 6d3da47b2, src tree == 6d3da47b2; INTEG4 cut prep 2026-11-21, awaiting the USER's push decision); the docs commit on top updates docs/clean_stag2/ (this file, the keys files, apply_v20.sh), on branch clean-cc-stag2-latest-v20 (push target origin/clean-stagger2-latest, fast-forward from v19 13ea625fe). Candidate = integ-trial-v20n 6d3da47b2 (INTEG4) = v20g 6b48f352c + restart I/O pieces <= 1 GiB + read-ahead buffer fix (240d9eb32 35a5b37e6) + CONSERV lzbook (sp L_z ledger, hst M_hflx E_hflx Lz_hflx Lz_tot; evolution bitwise) + H200 build knobs + P4 (xdonor diag only with clean_nl_print) + P5 code with ATHENA_CUDA_BODY_GLOBAL off by default (H200 +2-4 % with it) + LZCS/LZTEL2 clean_cs_lz_tel (default off) + force_reference=clean torque booked + ELEAK hydro/floor_pay (default off); v20g = v19 code fe6e94f03 + clean_conv_kappa default 1.0 (9c92823f5) + preflight "no effect" list for the 5 inert keys (49200b4c1) + GEOM momentum unification deb308830 619c3973d 694a3f9aa + FOFC count print ef89bf2b7 + THINF clean_stag_thinf (key, default off) e844a8f1f cb61f83bd + GEOM cs curved-face radial area d7f90a658 + GEOM clean_ray_flux_in removed 08fbea5d9 + LAG ring fill0 the only form 23b99d15a + TIE cs face tie 6b48f352c.
- **One scheme key set for HG, RSG and R136 (UNIFORM, USER 2026-11-20; stag2/uniform/UNIFORM.md):** every keys file = uniform/keys/common_scheme.txt + uniform/keys/<star>_problem.txt (applied by integ/keys/apply_v20.sh). Every scheme value is named explicitly, because restarts keep their old defaults (recipe::NewDefaults is false on a restart). Changes vs v19: `rad_m1/clean_conv_kappa=1.0` (also the code default now; RSG -21 %, HG -6 %, R136 insensitive), `rad_m1/time2_k1_stiff=0` (all stars, USER; HG had 1.0, ~7 % cost), `time/restart_refill_ghosts=false`, `rad_m1/clean_rot_fuse=true` for all stars (bitwise to false; memory/speed only).
- **Retired keys (no effect under scheme = clean; ImplicitSolve is not registered):** time2_adm_poslim, time2_dtr_poslim, zalesak_cap_rounds, vet_col_lat_posfix, time2_vet_x2_fields (now in the preflight "no effect" list), plus time2_vet_col, time2_accept_fallback, time2_posfix, vet_rotate_avg_ema*, vet_gd_twin*, vet_source_noesrc, implicit_*. Inert pairs on A100 (HG, R136): hst bitwise with and without them.
- **Rules:** tfloor = max(1000 K, the opacity table's bottom). Hardware keys (uniform/keys/hardware_rules.txt) are not scheme keys: `time/dt_retry_host=m1` only when the device peak exceeds ~90 % (HG 12 ranks / 3 A100 nodes); name it in the job's hardware deck.
- **Momentum unification (GEOM, USER):** u0 M1_F = the face momenta's cell vector on every mesh (one write-back); the faces' Phi cell flux is written as m1_phi1..3 after m1_f1..3 in the "m1" output (staggered). clean_ray_flux_in is removed (momentum is the only form; naming it is FATAL). Re-baseline: sp/Cartesian state at the 1e-8 (hydro) to 1e-5 (L_top) level, cs unchanged; bench output-only rows (m1 dumps, rhse_sp user.hst).
- **Deferred (USER):** K8 / rotation clock, geo3 rays, noinflow outer BC.
- **In v20, key off:** THINF clean_stag_thinf (beam fix; not in the keys files).
- **Code-only re-baselines in v20:** cs curved-face radial area (cs rows), ring fill0 the only form (fresh ring/rot decks; clean_rot_ring_fill0=false FATAL), cs face tie (rule-1 energy; multi-block cs round-off).
- **Not in v20:** REALIZ fix 3 (parked), clamp, wall a_ref pgen fix.
- **Resolved (zero-CW rule):** the CPU fresh-start R136 deck r136cL_c showed 2 CLEAN WARNINGs (stage fallback, cycles 1 and 3) because the test deck did not name rad_m1/time2_face_fct (off for scheme=clean unless named; production sets true). With it: 0 CW, 0 fallbacks, stage-start E<0 = 0 (CONV2; rsgs_tf 24 cycles also 0 CW at kappa 0.01/0.1/1). The v20 bench keys k20 and the k20 star decks name it.
- **Default-off keys with a USER recommendation pending:** hydro/floor_pay=rad_column (ELEAK: closes the RSG energy budget to noise, 0 CW, no slowdown); rad_m1/clean_cs_lz_tel=true (LZTEL2 fix: RSG 300 cycles 0 CW, L_z drift 8e-12; HG dt unchanged).
- **Bench keys k20** = k18 + hydro/fofc=true + mesh/nghost=3 + rad_m1/time2_face_fct=true (reference Sk20r2).
- **Key-on / default recommendations pending the USER (all OFF / unchanged in v20 as cut):** (1) hydro/floor_pay=rad_column (ELEAK; TQS localised the rule-1 energy remainder to the hydro EOS floor's created energy, which this closes to noise: AB9/AB9rt RSG, AB10 HG 0 CW, no slowdown); (2) rad_m1/clean_cs_lz_tel=true (LZTEL2 fix: RSG 300 cycles key on 0 CW, L_z drift 8e-12; HG 50 cycles dt identical to key off; rt bitwise); (3) code default rad_m1/time2_face_fct=true (production names it; tests that omit it ran without FCT: CPU stalls + R136 fresh 2 CW, CONV/CONV2); (4) code default rad_m1/time2_k1_stiff=0 (all keys files name 0).
- **Rule-1 blockers cleared:** TIE cs face tie exact on production (TQK); energy remainder = hydro EOS floor energy (TQS) -> floor_pay recommendation above. CONSERV lzbook books sp L_z at open boundaries + hst M/E/L_z boundary fluxes. LAG fill0 R136 fresh pair bitwise.
- **Gate (v20n 6d3da47b2):** CPU k20 bench (k18 + fofc + nghost=3 + time2_face_fct) 66 hst == Sk20r2 (the new k20 reference, reproducible k20ref == k20refb), bins data-identical, GATES FAIL 8 (known prod-value/ref-missing rows), CW 1 = forced flz qFB deck; stars plain == poison 0 CW; pzgate 10/10; restart I/O pieces: rst md5 == v20g incl. a 2^10-piece build, k = 1 read-ahead restart bitwise. Caltech H200 (4 ranks/node, romio): RSG rt 5/5, poison == plain == repro, HG rt 5/5 (the >= 2 GiB piece case that failed on v20g), 0 CW, hst chain V20G == V20J (cols 1-26) == V20L == V20M == V20N; speed RSG 2.84 HG 4.18 s/cycle (v20l vs v20g SPL 2.71 vs 2.73). Raven A100 gpudev R136 rt bitwise 0 CW (hst == V20K).
- **hst:** he_star_m1 user hst has 4 new columns at the end (M_hflx, E_hflx, Lz_hflx, Lz_tot; 30 columns); columns 1-26 and hydro hst are bitwise to v20g.


### v19 block (superseded by v20 above)
**Code:** LATEST v19 = **3225da996** (code; full 3225da996b7141312be92b0d9ad491219feadeb6; USER go 2026-11-20, INTEG3), on branch clean-cc-stag2-latest-v19; the docs commit on top adds docs/clean_stag2/ (this file + the keys files) (push target origin/clean-stagger2-latest). = merge of v18 6fbe1532c (tree unchanged) into fe6e94f03 = SPEEDHEAD 8f3c095ba (8f023a658 + REALIZ d96bbc478 + FRESH dffcb17a9 3351ee5ed + HALOX chain) + GEOM a80cc1336 (B1-B3, bitwise refactor) + 458df0a0a (R1) + HistSums 8c46615f7 (= 79a074719) + cs L areas 868771d51 (= 458de9f68) + donor-form removal fe6e94f03 (= 0eb0c13a9).
- **Robustness bundle is the cubed-sphere DEFAULT (R1, USER):** clean_cs_lz_abs=2 (gas phi momentum from the faces' absorption, its L_z booked in the CLEAN LZ ledger) and clean_stag_csvec_radial (cell vector exact for a radial field) are the only cs form; the old torque-row path is gone and naming the old values is FATAL. No keys needed. Evidence: A100 RSG 300 cycles 0 CW / 0 DT_RETRY, -32 % wall vs control (lz + csvec); HG 100 cycles 0 CW.
- **Keys change vs v18:** HG `rad_m1/clean_rot_fuse=true` (HG only). Column sweep NOT used (time2_face_fct_sweep stays false). clean_conv_kappa at the code default 0.1 (no key; the kappa=1.0 rule below still holds). Applied by integ/keys/apply_v19.sh (SWEEP=0).
- **Removed:** time2_face_fct_form=donor (FATAL if named); REALIZ fix 3 clean_stag_freal_kin not included (parked).
- **hst re-baseline (one, at v19):** (1) HALOX c126b6af5: ranks' history sums in rank order (MPI_Gather) -> multi-rank hst last bits differ from v18 / 8f023a658 (A100 RSG/HG: max rel 3.3e-15; state bitwise); (2) HistSums: the sums in a fixed cell order, no longer OMP_NUM_THREADS dependent; (3) 458de9f68: cubed-sphere L_in / L_top use the clean face areas (hst only); (4) cs runs: the R1 forms change the cs state (intended; == v18-era code with both keys on, bitwise). Compare v19 hst to v19 hst only.
- **Bench row:** cs feed1 L_z total includes the CLEAN LZ ledger (int2/gate/reg/cs_tests/ct.py).
- **Gate (v19 = fe6e94f03 tree):** CPU trimmed bench: all bins bitwise to G 458df0a0a (354/354), which is 66/66 hst bitwise to v18-era code + both cs keys, sp/Cartesian bitwise to keys-off; 15 known FAIL rows unchanged; stars plain+poison; pzgate 10/10 (poison 1/2/3, ranks 4=12, restart 6->12); A100 R136 (sp) bitwise; A100 3-node RSG rsg.01510 + HG hegiant.00184 (v19 HG keys): restart round trip bitwise, poison == plain, 3 repeated RSG runs hst identical, 0 CW; RSG 2.62 s/cycle (8f023a658 keys-off: 4.52).

(previous header) **Code:** LATEST v18 = **6fbe1532c** (full 6fbe1532c456887848c9314e751f254f5d5cae21), on branch clean-cc-stag2-latest-v18 (push target origin/clean-stagger2-latest). It supersedes v17b 6bb6d0844 and earlier. **Keys change in v18 (USER: the old prod limiter form): rad_m1/time2_k1_redist=0 for HG/RSG (was 0.1; R136 already 0), backups *.bak_v17b.** rad_m1/time2_face_fct_sweep is in the code and stays OFF. The v17b key decisions (xdonor=lag, xsrc=true, gas_keep=1, gas_cons=true, mg_color=mirror, zalesak_active, layout_invariant) stand. The keys are the same for all three, except vet_gd_band_ring=false for R136.
**HG 3-node variant (USER decision 2026-11-20, HGMEM):** HG on 12 ranks / 3 A100 nodes = the HG LATEST keys + `time/dt_retry_host=m1` (dt_retry clones of the M1 state in pinned host memory; A100 HG 12r fit + rt BITWISE with m1, HGMEM). Default (4+ nodes) keeps dt_retry_host at its default. The key exists since SPEEDHEAD b863534e8 (HGMEM 7128fb2f8); the 3-node memory fit also needs HGMEM 621450c63 + 7a640c113 (SPEEDHEAD after 7b3003931); in LATEST from v19.
**clean_conv_kappa RULE (USER 2026-11-20): `rad_m1/clean_conv_kappa=1.0` must never be used without `rad_m1/time2_face_fct_sweep=true`** (A100: kappa without the sweep cascades, 7b3003931 RSG dt 1000 -> 62.5, 14 CW; with the sweep clean over 75 cycles). kappa is OUT of v19 (sweep only); revisit after the DeltaAI kappa+sweep long run.
**v19 keys plan (USER 2026-11-20, applied at the v19 cut by integ/keys/apply_v19.sh):** column sweep DROPPED (USER, later on 11-20: time2_face_fct_sweep stays false for HG and RSG; MI300A +3-5 %, A100 RSG +1.4 %, no robustness gain); `rad_m1/clean_rot_fuse=true` (HG only; RAYS 3f8e5e37f, default off on all backends; A100 HG 12r 39.48 GB, rt bitwise, hst = off, HG 24r -2.5 % mid); the HG 3-node variant (+ time/dt_retry_host=m1) inherits it. Not RSG / R136.
**History sums (INTEG 2026-11-20):** up to v18 the hst columns depend on OMP_NUM_THREADS at round-off (the state does not: dumps bitwise 8 vs 12 threads). Fix 79a074719 (branch clean-cc-stag2-histsum, repro::HistSums, fixed order) enters at the v19 cut (one hst re-baseline). KNOWN LIMITATION after it: hst sums still depend on the rank count / MeshBlock-to-rank layout (fixed order per rank + MPI_Reduce); compare hst only at equal rank counts.
**Radiation output (v19, USER spec 2026-11-20):** see stag2/OUTPUT.md. Under clean, m1_vet = the full Eddington tensor m1_d11 d22 d33 d12 d13 d23 of the closure the step's solves used, in the frame of m1_f1..3 (radial first; cubed sphere rotated from Cartesian, unlike old prod). P = m1_e * m1_d. Launch recipe for the 12-rank/3-node HG host map: stag2/halox/PLACEMENT.md (HALOX; -1.5 %, bitwise; 24-rank/6-node default is already optimal).
(previous header) **Code:** LATEST v17b = **6bb6d0844**, on branch clean-cc-stag2-latest-v17b (it supersedes v16 and earlier; v11/v12 do not compile for CUDA). **Keys change in v17b (USER decisions 2026-11-19): clean_stag_xdonor=lag, clean_stag_xsrc=true, implicit_pos_floor_gas_keep=1, implicit_pos_gas_cons=true, clean_stag_mg_color=mirror** (backups *.bak_v16). The R136 A100 rows of v9-v11 are PENDING until Raven is back; rule 2 is closed on MI300A. The keys are the same for all three, except vet_gd_band_ring=false for R136.

**Machine-readable files** (one key=value per line, `key=value` tokens; append run-control keys yourself):
- stag2/hg_latest.txt
- stag2/rsg_latest.txt
- stag2/r136_latest.txt

R136: r136_latest.txt is for restarts from **r136.00097** (Raven /raven/ptmp/jinma/stag2/r136rst). That file's own header has nx1 = 738, the stretched grid and IC ic_r136a1_B_dfl17_ext165, so the mesh/IC keys repeat the file's values. The '831' deck (r136_half_rt_sc.athinput, 831x128x64 fresh start) is the A100 RT gate deck only, not a production key set. Table paths are the viper paths (HG/RSG) and the Raven paths (R136); DeltaAI has its own copy, delta/rsg_latest_delta.args.

## How the keys were made
- **Source:** the K1 production args.
  - HG / RSG: viper vgpu/k1_{hg,rsg}.args, the base args that VGPU gives every gpuq run.
  - R136: Raven gh200lead/r136_keys.txt + the kline rarms.sh staggered set. Where a key appeared twice, the last value wins, so time2_k1_redist=0 and he_ic_file=..._ext165.
- **Decisions applied:**
  - **clean_xt_exact=cons** (S3B).
  - **implicit_fb_local / _rounds dropped:** the code is deleted from the clean solve (USER 2026-11-19).
  - **clean_stag_kpos not named** (default off). Both S3B's kpos and S3P's P4 7fe820fbc are parked: P4 is bitwise with kpos off, and kpos=on is not usable yet.
  - **implicit_pos_floor_cons=true**, the K1 value. It is **PENDING S3B_hg767342_pg**: flip HG+RSG to false only if that HG smoke holds dt 97.5 with 0 CW. RSG with false is already proven (S3B_rsg684506_head: dt 888.9, 0 CW).
  - **implicit_pos_gas=true**, named explicitly. This is the K1 default and also the value in every prod restart file.
  - **MPI-IO independent** = the default of job/mpiio_collective=false (RSTIO). No key is needed. Name job/mpiio_collective=true only to read or write collectively.
  - **vet_gd_band_ring:** false for R136/A100 only (SPD2: bitwise, -47 ms/cycle). HG/RSG keep true (the default; false costs +25 GB/rank on RSG).
- **MEMTRIM defaults:** no key. clean_ring_fields stays at its default 'all'.
- **Dropped keys:**
  - The rad_m1/implicit_* keys that CleanPreflight lists as no effect: it removes them from a fresh run's input, and ABLATE P16 showed that naming them is bitwise. A restart takes its file's value anyway.
  - time2_inc_kecap: time2.cpp:429 reads it only without clean (ABLATE k2).
- **KEPT although CleanPreflight calls them "no effect":** time2_vet_col, time2_accept_fallback, time2_posfix, time2_k1_stiff, time2_k1_redist, time2_face_fct and time/dt_retry_t2_be.
  - All of them are still read under clean (rad_m1_time2.cpp 282/542/577/629/666/677, dt_retry.cpp 488).
  - A fresh hgs run with them dropped DIFFERS from K1 at cycle 2 (keq, INTEG). The preflight label is wrong for at least one of them; this was reported to the coordinator as a rule-12 item.
  - Verified (stag2/integ/keq, binary 3007b797b): hg_latest on the hgs deck and rsg_latest on rsgs_tf, fresh 3 cycles, give hst BITWISE vs the full K1 args + the same decisions. hgs twice = identical. Both RSG arms show the same 3 CW (dt_retry NON-CONVERGED at cycle 1 of a fresh small-deck start), so K1 and latest agree on it.
- **GPU node counts (MI300A, v3):** K-line RSG fits on 2 nodes ONLY WITHOUT restart output (MEMTRIM: peak 98.9 GB/rank). RSG runs that write a restart OOMed on 2 nodes up to v4 (needed 3); from v6/v7 on (RSTIO zero-copy + windows, MEMTRIM O) restart-writing RSG FITS on 2 nodes (INTEG_rsg684506_v7); dump staging is the next memory target (MEMTRIM/RSTIO). HG: 4 nodes, incl. rtcheck.
- **GPU facts (coordinator):**
  - The RSTIO I/O line on RSG (lib 684506, 6 MI300A): start-up 24-32 s, final 124 GB dump 16 s.
  - K1 on rsg.01658 held dt 1000 (0 CW, 0 floor clips).

### HG (hg_latest.txt)

| key | value | reason / decision source |
|---|---|---|
| rad_m1/zalesak_active | `true` | ADDED v15: bitwise to za off (CPU bench 66/66 hst + MI300A HG/RSG hst identical to v14), adopted for speed (HG -20 %) |
| rad_m1/time2_vet_col | `extrap2` | KEPT: read under clean (time2.cpp:282) though preflight names clean_closure; unverified dead |
| rad_m1/time2_accept_fallback | `true` | KEPT: read under clean (time2.cpp:677); unverified dead |
| rad_m1/time2_posfix | `4` | KEPT: read under clean (time2.cpp:666); unverified dead |
| ~~rad_m1/implicit_work_fix~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/cs_vertex_mono | `true` | K1 production value (hg) |
| rad_m1/time2_auto_be | `16` | K1 production value (hg) |
| rad_m1/time2_k1_stiff | `1.0` | KEPT: read under clean (time2.cpp:542); unverified dead |
| rad_m1/vet_gd_twin | `false` | K1 production value (hg) |
| rad_m1/vet_gd_symref | `none` | K1 production value (hg) |
| rad_m1/vet_gd_twin_fuse | `false` | K1 production value (hg) |
| rad_m1/vet_gd_twin_state | `mean` | K1 production value (hg) |
| rad_m1/vet_source_noesrc | `false` | K1 production value (hg) |
| time/dt_retry_host | `false` | K1 production value (hg) |
| rad_m1/vet_gd_gl_nmu | `0` | K1 production value (hg) |
| rad_m1/vet_gd_geo | `2` | K1 production value (hg) |
| rad_m1/vet_rotate_avg | `8` | K1 production value (hg) |
| rad_m1/vet_rotate_advance | `true` | K1 production value (hg) |
| mesh/cs_seam_blockend | `extend` | K1 production value (hg) |
| mesh/cs_vertex_3cell | `true` | K1 production value (hg) |
| time/cfl_number | `0.3` | K1 production value (hg) |
| problem/he_ic_file | `/viper/ptmp2/jinma/prod_1107/tables/hgglob_1021/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt` | K1 production value (hg) |
| problem/he_opac_table | `/viper/ptmp2/jinma/prod_1107/tables/hgglob_1021/files/hg145-prod-1018/rosseland_mesakap_afe2dd8da2bd.txt` | K1 production value (hg) |
| problem/he_planck_table | `/viper/ptmp2/jinma/prod_1107/tables/hgglob_1021/files/hg145-prod-1018/planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt` | K1 production value (hg) |
| rad_m1/cs_cross_mono | `true` | K1 production value (hg) |
| ~~rad_m1/implicit_opac_newton~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| mesh/cs_seam_rho_guard | `4` | K1 production value (hg) |
| rad_m1/implicit_predictor_tcap | `1.1` | K1 production value (hg) |
| rad_m1/time2_k1_redist | `0.1` | KEPT: measured LIVE under clean (keq: dropping it changes hgs from cycle 2) although CleanPreflight lists time2_k1_* as no effect (reported) |
| ~~rad_m1/implicit_fb_local~~ | ~~`true`~~ | DROPPED: DELETED from the clean solve (USER 2026-11-19, S3B a8b7a44fc); no effect |
| ~~rad_m1/implicit_fb_local_rounds~~ | ~~`4`~~ | DROPPED: with implicit_fb_local (deleted) |
| ~~rad_m1/implicit_opac_lag_passes~~ | ~~`3`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/vet_gd_tangential | `true` | K1 production value (hg) |
| mesh/layout_invariant | `true` | K1 production value (hg) |
| hydro/wall_wb_onesided | `true` | K1 production value (hg) |
| hydro/split_c2p_once | `true` | K1 production value (hg) |
| rad_m1/vet_gd_rotate_dt | `1170` | K1 production value (hg) |
| rad_m1/vet_gd_rotate_align | `true` | K1 production value (hg) |
| rad_m1/vet_gd_rotate_axis | `3d` | K1 production value (hg) |
| rad_m1/vet_rotate_avg_ema_change_form | `step` | K1 production value (hg) |
| ~~rad_m1/implicit_tol~~ | ~~`1e-10`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_am_cs_kecap~~ | ~~`0.5`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_clip_kecap~~ | ~~`0.5`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| time/dt_retry | `true` | K1: dt_retry on (checkpoint + retry) |
| rad_m1/time2_vet_x2_fields | `15` | K1 production value (hg) |
| ~~rad_m1/implicit_face_mmat~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/vet_col_lat_posfix | `lag` | K1 production value (hg) |
| ~~rad_m1/implicit_cs_e0_metric~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_cs2_wbsolve~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| hydro/eos_logt_max | `8.0` | K1 production value (hg) |
| ~~rad_m1/implicit_hr_cmem~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_physflux~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_blend_dtc~~ | ~~`1e-6`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/zalesak_cap_rounds | `4` | K1 production value (hg) |
| ~~rad_m1/implicit_lzx~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_vimp_lag_local~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| time/tlim | `1e30` | K1 production value (hg) |
| time/ndiag | `1` | K1 production value (hg) |
| output2/dt | `1e30` | K1 production value (hg) |
| output3/dt | `1e30` | K1 production value (hg) |
| output4/dt | `1e30` | K1 production value (hg) |
| output5/dt | `1e30` | K1 production value (hg) |
| output6/dt | `1e30` | K1 production value (hg) |
| rad_m1/scheme | `clean` | K line: clean solve |
| rad_m1/clean_layout | `staggered` | K line: staggered |
| rad_m1/clean_rays | `sc` | K line: sc rays |
| rad_m1/clean_init | `given` | K1: given on restarts (HG), steady (RSG) |
| rad_m1/vet_gd_rotate_xfade | `true` | K line: ring cross-fade |
| rad_m1/clean_rot_ring | `live` | K line: live rotating ring |
| rad_m1/time2_face_fct | `true` | KEPT: read under clean (time2.cpp:629, Time2FaceFct); preflight label no effect is wrong/unverified |
| hydro/tfloor_kelvin | `1000` | K1 floor (HG/RSG 1000, R136 5000) |
| hydro/eos_logt_min | `2.9` | K1 production value (hg) |
| rad_m1/vet_rotate_avg_ema | `128` | v20.3 USER 2026-11-20: named (was `0` in the K1 args; unnamed and called inert since v20, so a restart took its file value: hegiant.00184 / rsg.01510 carry 128 + step, rsg.01295 / hegiant.00139 ran the K window). EMA ring over the rotated sets; vet_rotate_avg is then not used |
| rad_m1/clean_cs_lz_tel | `true` | v20.3 USER 2026-11-20 ~14:55 PST: LZTEL3 fix 74c8a8bb1 (cubed sphere, staggered; FATAL otherwise); net lateral-stress torque defect -> least-KE solid-body spin; counter-spin caveat in the v20.3 block |
| ~~rad_m1/time2_inc_kecap~~ | ~~`0.5`~~ | DROPPED: DELETED (ABLATE k2 ba7703940/0c7bb3041, USER); time2.cpp:429 reads it only without clean |
| rad_m1/time2_adm_poslim | `false` | K1 production value (hg) |
| rad_m1/time2_dtr_poslim | `false` | K1 production value (hg) |
| rad_m1/clean_xt_exact | `cons` | cons: S3B decision (exact linearised exchange + work, closure prints); RSG GPU S3B_rsg684506_head dt 888.9 0 CW |
| rad_m1/implicit_pos_gas | `true` | true = K1 default, named explicitly (coordinator 11-19); also in every prod restart file; 061849266 CW fires only when off |
| rad_m1/implicit_pos_floor_cons | `true` | true = K1 (held HG dt 97.5); PENDING S3B_hg767342_pg: flip HG+RSG to false only if HG holds dt 97.5 with 0 CW (RSG false proven: S3B_rsg684506_head) |

### RSG (rsg_latest.txt)

| key | value | reason / decision source |
|---|---|---|
| rad_m1/zalesak_active | `true` | ADDED v15: bitwise to za off (CPU bench 66/66 hst + MI300A HG/RSG hst identical to v14), adopted for speed (HG -20 %) |
| rad_m1/time2_vet_col | `extrap2` | KEPT: read under clean (time2.cpp:282) though preflight names clean_closure; unverified dead |
| rad_m1/time2_accept_fallback | `true` | KEPT: read under clean (time2.cpp:677); unverified dead |
| rad_m1/time2_posfix | `4` | KEPT: read under clean (time2.cpp:666); unverified dead |
| time/dt_retry_t2_be | `0` | KEPT: read by dt_retry.cpp:488; unverified dead |
| ~~rad_m1/implicit_work_fix~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/cs_vertex_mono | `true` | K1 production value (rsg) |
| rad_m1/time2_auto_be | `16` | K1 production value (rsg) |
| rad_m1/vet_gd_twin | `false` | K1 production value (rsg) |
| rad_m1/vet_gd_symref | `none` | K1 production value (rsg) |
| rad_m1/vet_gd_twin_fuse | `false` | K1 production value (rsg) |
| rad_m1/vet_gd_twin_state | `mean` | K1 production value (rsg) |
| rad_m1/vet_source_noesrc | `false` | K1 production value (rsg) |
| time/dt_retry_host | `false` | K1 production value (rsg) |
| rad_m1/vet_rotate_advance | `true` | K1 production value (rsg) |
| mesh/cs_seam_blockend | `extend` | K1 production value (rsg) |
| hydro/merge_lat_blockend | `exchange` | K1 production value (rsg) |
| mesh/cs_vertex_3cell | `true` | K1 production value (rsg) |
| problem/he_ic_file | `/viper/ptmp2/jinma/prod_1107/tables/rsg_E_1020/files/rsg-E-1020/ic_rsg_g25_deep.txt` | K1 production value (rsg) |
| problem/he_opac_table | `/viper/ptmp2/jinma/prod_1107/tables/rsg_E_1020/files/rsg-g25-prod/rosseland_mesakap_4c0e8b64d26a.txt` | K1 production value (rsg) |
| problem/he_planck_table | `/viper/ptmp2/jinma/prod_1107/tables/rsg_E_1020/files/rsg-g25-prod/planck_ext2_x0.653848_z0.019757_d38e5526afef.txt` | K1 production value (rsg) |
| rad_m1/cs_cross_mono | `true` | K1 production value (rsg) |
| ~~rad_m1/implicit_opac_newton~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| mesh/cs_seam_rho_guard | `4` | K1 production value (rsg) |
| rad_m1/implicit_predictor_tcap | `1.1` | K1 production value (rsg) |
| rad_m1/time2_k1_redist | `0.1` | KEPT: measured LIVE under clean (keq: dropping it changes hgs from cycle 2) although CleanPreflight lists time2_k1_* as no effect (reported) |
| ~~rad_m1/implicit_fb_local~~ | ~~`true`~~ | DROPPED: DELETED from the clean solve (USER 2026-11-19, S3B a8b7a44fc); no effect |
| ~~rad_m1/implicit_fb_local_rounds~~ | ~~`4`~~ | DROPPED: with implicit_fb_local (deleted) |
| ~~rad_m1/implicit_opac_lag_passes~~ | ~~`3`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/vet_gd_tangential | `true` | K1 production value (rsg) |
| mesh/layout_invariant | `true` | K1 production value (rsg) |
| hydro/wall_wb_onesided | `true` | K1 production value (rsg) |
| hydro/split_c2p_once | `true` | K1 production value (rsg) |
| time/ndiag | `1` | K1 production value (rsg) |
| time/tlim | `1e30` | K1 production value (rsg) |
| output2/dt | `1e30` | K1 production value (rsg) |
| output3/dt | `1e30` | K1 production value (rsg) |
| output5/dt | `1e30` | K1 production value (rsg) |
| output6/dt | `1e30` | K1 production value (rsg) |
| problem/he_ic_equil | `false` | K1 production value (rsg) |
| rad_m1/scheme | `clean` | K line: clean solve |
| rad_m1/clean_layout | `staggered` | K line: staggered |
| rad_m1/clean_rays | `sc` | K line: sc rays |
| rad_m1/clean_init | `steady` | K1: given on restarts (HG), steady (RSG) |
| rad_m1/vet_gd_rotate_xfade | `true` | K line: ring cross-fade |
| rad_m1/clean_rot_ring | `live` | K line: live rotating ring |
| rad_m1/time2_face_fct | `true` | KEPT: read under clean (time2.cpp:629, Time2FaceFct); preflight label no effect is wrong/unverified |
| hydro/tfloor_kelvin | `1000` | K1 floor (HG/RSG 1000, R136 5000) |
| hydro/eos_logt_min | `2.9` | K1 production value (rsg) |
| rad_m1/vet_rotate_avg_ema | `128` | v20.3 USER 2026-11-20: named (was `0` in the K1 args; unnamed and called inert since v20, so a restart took its file value: hegiant.00184 / rsg.01510 carry 128 + step, rsg.01295 / hegiant.00139 ran the K window). EMA ring over the rotated sets; vet_rotate_avg is then not used |
| rad_m1/vet_rotate_avg_ema_change_form | `step` | v20.3 USER 2026-11-20: named with ema 128 (= the restart files' value) |
| rad_m1/clean_cs_lz_tel | `true` | v20.3 USER 2026-11-20 ~14:55 PST: LZTEL3 fix 74c8a8bb1 (cubed sphere, staggered; FATAL otherwise); net lateral-stress torque defect -> least-KE solid-body spin; counter-spin caveat in the v20.3 block |
| ~~rad_m1/time2_inc_kecap~~ | ~~`0`~~ | DROPPED: DELETED (ABLATE k2 ba7703940/0c7bb3041, USER); time2.cpp:429 reads it only without clean |
| rad_m1/time2_adm_poslim | `false` | K1 production value (rsg) |
| rad_m1/time2_dtr_poslim | `false` | K1 production value (rsg) |
| rad_m1/time2_stage_local_be | `false` | false: USER rule (stage-local BE leaks energy) |
| rad_m1/clean_xt_exact | `cons` | cons: S3B decision (exact linearised exchange + work, closure prints); RSG GPU S3B_rsg684506_head dt 888.9 0 CW |
| rad_m1/implicit_pos_gas | `true` | true = K1 default, named explicitly (coordinator 11-19); also in every prod restart file; 061849266 CW fires only when off |
| rad_m1/implicit_pos_floor_cons | `true` | true = K1 (held HG dt 97.5); PENDING S3B_hg767342_pg: flip HG+RSG to false only if HG holds dt 97.5 with 0 CW (RSG false proven: S3B_rsg684506_head) |

### R136 (r136_latest.txt)

| key | value | reason / decision source |
|---|---|---|
| rad_m1/zalesak_active | `true` | ADDED v15: bitwise to za off (CPU bench 66/66 hst + MI300A HG/RSG hst identical to v14), adopted for speed (HG -20 %) |
| rad_m1/vet_gd_band_ring | `false` | ADDED for A100/R136 only (SPD2 2026-11-19): key-only, R136 hst md5 BITWISE, -47 ms/cycle (C1 0.953 -> 0.903 s), memory 24.6 of 40 GB; HG/RSG keep the default vet_gd_band_ring=true: on RSG false costs +25 GB/rank (93.2 vs ~69 GB, MEMTRIM) |
| problem/he_opac_table | `/raven/ptmp/jinma/stag2/gh200lead/pkg/opac/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt` | K1 production value (r136) |
| problem/he_planck_table | `/raven/ptmp/jinma/stag2/gh200lead/pkg/opac/planck_tops_x0.58_z0.0064_a09_T7.txt` | K1 production value (r136) |
| rad_m1/time2_adm_poslim | `true` | K1 production value (r136) |
| rad_m1/cs_cross_mono | `true` | K1 production value (r136) |
| ~~rad_m1/implicit_opac_newton~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| mesh/cs_seam_rho_guard | `4` | K1 production value (r136) |
| rad_m1/implicit_predictor_tcap | `1.1` | K1 production value (r136) |
| ~~rad_m1/implicit_fb_local~~ | ~~`true`~~ | DROPPED: DELETED from the clean solve (USER 2026-11-19, S3B a8b7a44fc); no effect |
| ~~rad_m1/implicit_fb_local_rounds~~ | ~~`4`~~ | DROPPED: with implicit_fb_local (deleted) |
| ~~rad_m1/implicit_opac_lag_passes~~ | ~~`3`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/vet_gd_tangential | `true` | K1 production value (r136) |
| mesh/layout_invariant | `true` | K1 production value (r136) |
| problem/he_wb_grav_window | `true` | K1 production value (r136) |
| output2/bin_double | `true` | K1 production value (r136) |
| output3/bin_double | `true` | K1 production value (r136) |
| output6/bin_double | `true` | K1 production value (r136) |
| ~~rad_m1/implicit_wb_vkick_reject~~ | ~~`1e9`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| mesh/layout_invariant_mode | `block` | K1 production value (r136) |
| problem/he_sponge_rate | `0` | K1 production value (r136) |
| problem/he_bc_outer | `noinflow` | K1 production value (r136) |
| output2/dt | `2000` | K1 production value (r136) |
| output3/dt | `2000` | K1 production value (r136) |
| rad_m1/floorcap_rho | `1e-14` | K1 production value (r136) |
| rad_m1/floorcap_v | `3e7` | K1 production value (r136) |
| rad_m1/floorcap_report | `25` | K1 production value (r136) |
| time/dt_retry | `true` | K1: dt_retry on (checkpoint + retry) |
| problem/he_ic_file | `/raven/ptmp/jinma/stag2/gh200lead/pkg/ic/ic_r136a1_B_dfl17_ext165.txt` | K1 production value (r136) |
| mesh/nx1 | `738` | K1 production value (r136) |
| mesh/x1min | `1.3531365000e+12` | K1 production value (r136) |
| mesh/x1max | `1.0435500000e+13` | K1 production value (r136) |
| mesh/use_grid_stretch_r_poly | `true` | K1 production value (r136) |
| mesh/f_stretch_r_c1 | `7.554612447083e+00` | K1 production value (r136) |
| mesh/f_stretch_r_c2 | `7.750389587028e+00` | K1 production value (r136) |
| mesh/f_stretch_r_c3 | `-2.971290112463e+00` | K1 production value (r136) |
| mesh/f_stretch_r_c4 | `4.271691922033e+01` | K1 production value (r136) |
| mesh/f_stretch_r_c5 | `-6.710490864599e+00` | K1 production value (r136) |
| mesh/f_stretch_r_c6 | `-8.763839994759e+01` | K1 production value (r136) |
| mesh/f_stretch_r_c7 | `-8.853219696626e+01` | K1 production value (r136) |
| mesh/f_stretch_r_c8 | `1.211561964156e+02` | K1 production value (r136) |
| mesh/f_stretch_r_b1_amp | `1.757003053390e+01` | K1 production value (r136) |
| mesh/f_stretch_r_b1_x | `7.000677767464e-01` | K1 production value (r136) |
| mesh/f_stretch_r_b1_w | `2.003554846351e-01` | K1 production value (r136) |
| mesh/f_stretch_r_b2_amp | `2.234756275431e+00` | K1 production value (r136) |
| mesh/f_stretch_r_b2_x | `8.297907168131e-01` | K1 production value (r136) |
| mesh/f_stretch_r_b2_w | `6.280572958365e-02` | K1 production value (r136) |
| mesh/f_stretch_r_p_amp | `-5.501037773041e+00` | K1 production value (r136) |
| mesh/f_stretch_r_p_xa | `-3.301071048718e+01` | K1 production value (r136) |
| mesh/f_stretch_r_p_xb | `7.096507145028e-01` | K1 production value (r136) |
| mesh/f_stretch_r_p_w | `6.546917834856e-02` | K1 production value (r136) |
| mesh/f_stretch_r_p2_amp | `-3.663380392738e-01` | K1 production value (r136) |
| mesh/f_stretch_r_p2_xa | `1.892800773211e-01` | K1 production value (r136) |
| mesh/f_stretch_r_p2_xb | `9.350578461191e-01` | K1 production value (r136) |
| mesh/f_stretch_r_p2_w | `3.480611062510e-02` | K1 production value (r136) |
| rad_m1/time2_vet_x2_fields | `15` | K1 production value (r136) |
| rad_m1/time2_k1_redist | `0` | KEPT: measured LIVE under clean (keq: dropping it changes hgs from cycle 2) although CleanPreflight lists time2_k1_* as no effect (reported) |
| ~~rad_m1/implicit_face_mmat~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/vet_col_lat_posfix | `lag` | K1 production value (r136) |
| hydro/eos_logt_max | `7.0` | K1 production value (r136) |
| ~~rad_m1/implicit_hr_cmem~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_physflux~~ | ~~`false`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| ~~rad_m1/implicit_blend_dtc~~ | ~~`1e-6`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/zalesak_cap_rounds | `4` | K1 production value (r136) |
| ~~rad_m1/implicit_lzx~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/time2_stage_local_be | `false` | false: USER rule (stage-local BE leaks energy) |
| ~~rad_m1/implicit_cs2_wbsolve~~ | ~~`true`~~ | DROPPED: no effect under clean: CleanPreflight drops the implicit_* keys the clean solve does not read (ABLATE P16: naming them is bitwise); a restart reads its file value anyway |
| rad_m1/scheme | `clean` | K line: clean solve |
| rad_m1/clean_layout | `staggered` | K line: staggered |
| rad_m1/clean_rays | `sc` | K line: sc rays |
| hydro/tfloor_kelvin | `5000` | K1 floor (HG/RSG 1000, R136 5000) |
| rad_m1/time2_face_fct | `true` | KEPT: read under clean (time2.cpp:629, Time2FaceFct); preflight label no effect is wrong/unverified |
| rad_m1/vet_gd_rotate_xfade | `true` | K line: ring cross-fade |
| rad_m1/clean_rot_ring | `live` | K line: live rotating ring |
| ~~rad_m1/time2_inc_kecap~~ | ~~`0.5`~~ | DROPPED: DELETED (ABLATE k2 ba7703940/0c7bb3041, USER); time2.cpp:429 reads it only without clean |
| rad_m1/clean_xt_exact | `cons` | cons: S3B decision (exact linearised exchange + work, closure prints); RSG GPU S3B_rsg684506_head dt 888.9 0 CW |
| rad_m1/implicit_pos_gas | `true` | true = K1 default, named explicitly (coordinator 11-19); also in every prod restart file; 061849266 CW fires only when off |
| rad_m1/implicit_pos_floor_cons | `true` | true = K1 (held HG dt 97.5); PENDING S3B_hg767342_pg: flip HG+RSG to false only if HG holds dt 97.5 with 0 CW (RSG false proven: S3B_rsg684506_head) |
## RULE CHECKLIST per LATEST (USER 2026-11-19); v1 = 3007b797b, v2 = fc2c2c8d4
| # | rule | evidence | v1 | v2 |
|---|---|---|---|---|
| 1 | energy / momentum / L_z / mass conservation at round-off | flz q t 32.88 dE/E 6.7e-16 dLz 3.9e-16 (K line had 7.3e-5: fixed by S3B 3df072231); flz Q 2.9e-16 / 9.4e-16; flz q reject=retry 6.7e-16; cs feed1 E 1.1e-16 Lz 0; cs closed sum EV 2.2e-16; pulse0 closed-box E 8.9e-16 / 4.4e-16; drag rows x7 PASS | green | green |
| 1 | (open items) | OPEN: **flz q face FCT** dE/E 1.3e-3 FAIL. It is a prod-value K-line row, but time2_face_fct=true IS a production key (read under clean, 96b6a3da9) -> owner S3B. OPEN: flz q BE only dE/E 8.5e-9, 50 CW, rc 1 (K-line prod-value row, = S3B gate) -> owner S3B. OPEN: no dedicated mass row in the trimmed gate (floors booked by the ledger rows only) -> owner INTEG to add | OPEN | OPEN |
| 2 | restart bitwise + dt_retry / forced retry bitwise | 11 RT rows (stg hcyl/sp3/flz, ringC/S, csCnr/csSnr, rot1 C/S, FCT, reject=retry) + 8 dt_retry rows fac 1 bitwise; hgs 4 ranks x 3 I/O modes (indep / coll / indep file read coll) dump + hst bitwise; ring hg_s a=b=c3/c4/c6; sp ring a=b=c3..c6; csSnr forced retry bitwise (v2, + bounds-check build rt + dtr bitwise, 0 hits); R136 A100 half: restart 5->10 and forced retry bitwise | green | green |
| 2 | (open items) | OPEN: GPU RSG restart ~1e-13 diff (owner RSTIO, investigating) | OPEN | OPEN |
| 3 | reproducibility | CPU: the hgs K1-key run twice -> identical hst (keq hlt2 = hlt_rd); every bench row = the previous gate. GPU: R136 A100 rerun pending (gpuq rtcheck straight arms) | CPU green, GPU OPEN (INTEG) | CPU green, GPU OPEN (INTEG) |
| 4 | positivity reported, not silent | stars print pos/mark/neg/guard counts (0); CLEAN OLDVEC / KPOS / floor counts in the GPU .res; S3B 061849266 counts eint-below-floor write-backs as CLEAN WARNING; pos_un_drop reported (and not wiped by the visitor) | green | green |
| 5 | well-balancing | rhse_sp 0.18, rhse_sp2 1.4e-9, rhse_cs6 5.6e-9, rhse_cs24 (4 ranks) 2.2e-8, wb.sh fv1/fv/lld/cart, wbpoly ideal 1.9e-16 / general 2.7e-6: all '=' | green | green |
| 6 | accuracy / order | radshock M2/M5/M7 N256 + M7 N512, thin_r/dif/lat N32, Marshak a-d, nltek1, hcyl, fs: all '='. The deck D order test (ordD) is not in the trimmed gate; WBSIDE ran it on 1645edc41 (rate 1.96), and no changed row touches it | green (ordD not rerun) | green (ordD not rerun) |
| 7 | symmetry | decomposition rows '=' (sp3 1 vs 2 blocks, cs thin1 vs thin1mb, mbmr 2 ranks, cellf 1 vs 2 blocks/ranks 0 FAIL). OPEN: no dedicated mirror-symmetry row -> owner INTEG to add | OPEN | OPEN |
| 8 | CPU/GPU parity | R136 A100 runs, but there is no CPU-vs-A100 compare row this version. MI300A: gpuq INTEG_*_latest2 (v1) pending; DeltaAI GH200 RSG smoke after the v2 build | OPEN (INTEG) | OPEN (INTEG) |
| 9 | no silent fallbacks | 0 CW on all bench rows except flz q BE only (known, row 1); stars 0 CW 0 retry; cwin 0. Preflight labels WRONG in v1/v2 (time2_k1_* / face_fct called no effect; INTEG keq found k1_redist live) -> fixed in S3B 96b6a3da9 = v3 | OPEN (labels) | OPEN (labels; fixed in v3) |
| 10 | build provenance | v1: Orion integ/integ_he md5 452f6ac0 / 12ec406f, src 70f4a3d5, rev 3007b797b; Raven integ_cu d2b31171. v2: Orion integ2 d478db66, integ2_he ea58bd29, integ2_bc 4ed808dc, src e8323b4c, rev fc2c2c8d4; Raven integ2_cu f1eb4b46 (.prov rev fc2c2c8d4). All binaries are fresh fastbuild slots with .prov; none were copied | green | green |

### RULE CHECKLIST v3 = 7497211f3 (stag2/integ/gate *_lat3*; Orion integ3 beeda6f4 / integ3_he 54d7a3ae / integ3_bc a8ef304c, src 99dedcef; Raven integ3_cu d49ca0dc)
| # | v3 | evidence |
|---|---|---|
| 1 conservation | green; **OPEN 1 (S3B, PRIORITY): flz q face FCT dE/E 1.3e-3. time2_face_fct=true is a PRODUCTION key, so this is a rule-1 violation on production keys.** OPEN: flz q BE only 8.5e-9 / 50 CW (S3B); no mass row (INTEG) | bench 69/69 '=' vs v2 (all conservation rows unchanged) |
| 2 restart / retry bitwise | green; OPEN: GPU RSG restart ~1e-13 (RSTIO; v1 and v2 MI300A rtchecks DIFF) | 11 RT + 8 dt_retry rows '='; csSnr forced retry + bounds-check rt/dtr bitwise (0 hits); hgs 3 I/O modes bitwise; ring/sp bitwise; R136 A100 straight = forced = restart bitwise |
| 3 reproducibility | CPU green; GPU: R136 A100 v3 straight = v2 straight BITWISE (all bins + hst; kpos off + SPD2 speed bitwise by construction) | |
| 4 positivity reported | green | stars counts 0 (identical to v2) |
| 5 WB | green | WB rows '=' |
| 6 order | green (ordD not rerun) | order rows '=' |
| 7 symmetry | OPEN (INTEG: no mirror row) | decomposition rows '=' |
| 8 CPU/GPU | OPEN (INTEG): MI300A HG/RSG on v3 queued; DeltaAI GH200 RSG smoke 3631287 (v2) pending | |
| 9 no silent fallbacks | preflight labels FIXED (S3B 96b6a3da9: time2_k1_*, face_fct*, dt_retry_t2_be listed as read); OPEN: flz q BE only 50 CW; rsgs_tf fresh start with the full K1 keys gives 3 CW (dt_retry NON-CONVERGED at cycle 1, identical with K1 and latest keys) | |
| 10 provenance | green | fresh fastbuild slots, .prov src 99dedcef = 7497211f3 |

### OPEN item 4 (INTEG) closed on v3 7497211f3, 2026-11-19
- **Mass row (new gate row 'mass (all hydro.hst)', int2/gate/sb/stab.py massrow):** max |dM/M| over all 35 hydro histories of the bench is 1.8e-15 (v2 and v3). The threshold is 1e-13. PASS.
- **Mirror-symmetry row (new row 'mirror x2 rhse_sp / sp2', tools/mirror.py):** this measures the x2 equator mirror asymmetry of the final m1 dumps of the symmetric rhse_sp / rhse_sp2 decks.
  - Staggered v3: **6.3e-8 / 1.8e-10**.
  - Cell layout (CEL base): **1.4e-14 / 6.2e-14** (round-off).
  - **So the staggered layout breaks mirror symmetry beyond round-off. Rule 10 is OPEN**; owner to be assigned (stag / S3x). The size is at the closure/solver-tolerance level; the likely suspects are the sc/fv1 direction set or the face-solve ordering.
  - The row passes at its current tolerance (1e-6 / 1e-8) so that it catches regressions.
- **GPU reproducibility:** the R136 A100 half-deck straight run (v3, pool 32490421, integ3_cu d49ca0dc) was repeated: all bins (hydro_w/m1/m1_vet at every output) and every hst row are BITWISE. PASS.
- **CPU vs A100 parity (rule 11):** R136 half deck reduced to 831x16x8 (1 block), 10 cycles. CPU is Orion integ3_he 54d7a3ae (1 rank x 16 threads); GPU is 1 A100 integ3_cu. The data files have identical md5.
  - dt is identical every cycle; 0 CW on both.
  - Per-column max|cpu-gpu|/max|col|: mass 2.6e-15, tot-E 4.4e-13, E_rad 3.3e-13, e_gas 1.7e-12, L_top 2.2e-14, L_bot 4.8e-8, 1-mom 5.4e-7, KE 1.3e-9.
  - Dumps after 10 cycles: dens/eint identical, vel <= 2.3e-8, m1_e 1.1e-11, m1_f 1.8e-5 (relative to the small transverse flux components).
  - The only large relative numbers are near-zero columns (2-mom 3.6e-3 of a ~0 column, P_esrc).
  - Far below the dt-halving error. PASS. Files: stag2/integ/parity/parity.txt.

### v3 MI300A results (viper gpuq, 2026-11-19)
- **HG** (INTEG_hg767342_v3, from hegiant.00184 to 767346, 4 nodes):
  - dt 97.5 every cycle, 0 CW, 0 retry, 8.95 s/cycle.
  - **rtcheck: final dump + hydro/user hst BITWISE.** Restart writing fits on 4 nodes. The HG GPU restart row is green.
  - Also closes the S3B_hg767342_pg question as far as K1 floor_cons=true goes: HG holds dt 97.5 with 0 CW.
- **RSG** (INTEG_rsg684506_v3, 684506 -> 684510; 3 nodes, after the 2-node OOM when writing a restart):
  - dt 888.9 every cycle, 0 CW, 0 retry, 7.6 s/cycle.
  - **rtcheck: final dump + hst DIFF.** This is the known OPEN rule-2 item, owner RSTIO.

### Owners update (coordinator)
- Mirror symmetry (staggered 6.3e-8 vs cell 1e-14) -> **SYMM**.
- flz q face FCT dE/E 1.3e-3 (production key time2_face_fct) and flz q BE only 50 CW -> **S3C** (S3B has left).

### Rule-3 item CLOSED (REPRO, 2026-11-19): MI300A HG v3 IS run-to-run reproducible when the keys match
- The two runs differed only in keys: INTEG passed clean_xt_exact=cons and SPD3 passed none. A rerun with INTEG's keys is BITWISE (REPRO_hg_v3_xtc).
- Root cause fixed: viper gpuq base args k1_{hg,rsg}.args were the OLD K1 keys. They are now stag2/{hg,rsg}_latest.txt (backups .bak_k1), so keyless gpuq requests run the LATEST keys.

### RULE CHECKLIST v7 = c5be823d2
- Every bench value is '=' vs v6. Rows 2/3/5/6/8/10 are as for v6.
- **Rule 9 OPEN (S3C):** 13 rows now FAIL on counted CW from the HLLC-NaN -> BE path (S3C a830f4765 makes it visible, 135 CW in the bench). The data are unchanged; the HLLC guard is a USER decision.
- **Rule-9 note (MEMTRIM3 question):** the sp ring row shows CLEAN WARNING 4 in arms a, b/a2, c3-c6 and h6 from v7 on; v6 had 0.
  - The warnings are S3C a830f4765's counted "stage start not finite in N cells (gas or E; the solve fails and the step is redone)", at cycle 7 t2s 3.
  - It is the same HLLC-NaN -> BE path as the 13 bench rows.
  - The sp hst/bins stay bitwise across all arms; the row passes on data.
  - Part of the S3C rule-9 OPEN item. It is not a MEMTRIM3 effect.

### MI300A v7 (viper gpuq)
- **RSG on 2 nodes, WITH restart writing (rtcheck):** it now fits (no OOM, 4 ranks) after RSTIO zero-copy + windows and MEMTRIM O. dt 888.9, 0 CW, 9.6 s/cycle. rtcheck final dump + hst DIFF (OPEN rule 2, RSTIO).
- **HG on 4 nodes:** dt 97.5, 0 CW, 8.1 s/cycle. rtcheck BITWISE.
- **Note:** Raven ssh master expired (2026-11-19 afternoon). v9's R136 A100 rows are PENDING until the USER reopens it.
- **Stage fallbacks on production keys (rule 9, S3P2 item):** the MI300A end summaries of INTEG_hg767342_v7 / v8 and INTEG_rsg684506_v7, including the rtcheck arms rs/r1/r2, show **stage fallbacks=0 and accepted-with-fallback stages=0**. So production HG/RSG on LATEST take no backward-Euler stage fallback.
  - A counted CLEAN WARNING for every such fallback (Time2Restore) is on branch clean-cc-stag2-latest-v10 (4a9fc8c54, report only) and will be gated in v10.
- **Gate-harness fix (INTEG):** int2/gate/stag_r.sh (backup stag_r.sh.bak_integ_rstrace) deletes the .rst dumps only after ALL parts have finished. Each part leaves a marker and the last one deletes.
  - Before the fix, a part that finished first deleted another part's restart dumps between its dump and its restart. That made v9's RT rot1 stag row fail with "Unable to open restart file".
  - The code was not at fault; the row was rerun on its own.
- **(Superseded in v15, see below)** Rule 12, mislabeled key: rad_m1/zalesak_active was DO NOT USE before v15. It is documented as bitwise but is NOT bitwise.
  - SPD4: the hgs hydro momentum differs by ~1e-2 relative with time2_face_fct=true, and the MI300A RSG hst differs.
  - It is not in any latest keys file. Owner: the zlocal-1113 code.

### Rule 2 (GPU RSG restart) CLOSED on v11 64b048253 (RSTIO 3882c2ce8 + 3fc90fae1 + 7a54f5162)
- **INTEG_rsg684506_v11:** 2 MI300A nodes, restart writing. dt 888.9, 0 CW, 0 retry, 9.55 s/cycle. **rtcheck: final dump + rsg.hydro.hst + rsg.user.hst BITWISE.**
- **INTEG_hg767342_v11:** 4 nodes. dt 97.5, 0 CW, 7.9 s/cycle. rtcheck BITWISE.
- All restart rows are now bitwise: CPU, A100 R136 (up to v8; v9-v11 pending Raven), MI300A HG and MI300A RSG.
- **OPEN (speed, not a rule):** v9's MEMTRIM3 CleanGeo made the RSG face kernels (rhsaf / fs0 / fwc / drag) about 13 % slower on MI300A. SPD4 is fixing this without giving up the memory saving.

### CUDA build broken in v11 and v12 (INTEG, nvcc compile check on Orion)
- Since RSTIO 3882c2ce8 (in v11 as 56964103f), the generic lambda write_group in src/outputs/restart.cpp has `const PackFn &pk = {}`. nvcc 12.1 rejects the brace-enclosed default argument, so restart.cpp does not compile for CUDA.
- **v11 64b048253 and v12 5330903d1 cannot be built for A100. Use v13 e4b2d568c or later** (fix: `= PackFn()`).
- HIP (MI300A) and CPU builds were unaffected.
- **CLOSED (speed):** the MEMTRIM3 CleanGeo face-kernel slowdown is fixed by SPD5 in v13. MI300A RSG 9.26 -> 7.79 s/cycle and HG 7.74 -> 7.16 s/cycle; memory peak unchanged (83873.5 MiB, SPD5).
- **CUDA compile check (nvcc 12.1, sm_80, Orion, no MPI):** v13 e4b2d568c: full build links, rc 0.
- **CUDA compile check:** v14 510d00b0a full build links, rc 0.

### rad_m1/zalesak_active=true ADOPTED in v15 (bb4735136): bitwise to za off, adopted for speed
- SPD6 fixed the race and made the key bitwise.
- Evidence:
  - The CPU bench run twice: all 66 history files identical to v14 with za off and with za on (incl. the face FCT rows).
  - The za-on arm's restart, retry, stars, ring, sp and cellf rows also equal v14.
  - MI300A HG and RSG with za on: rtcheck bitwise and hst byte-identical to the v14 runs.
- Speed: HG 7.07 -> 5.63 s/cycle, RSG 7.45 -> 7.17 s/cycle.
- The key is in hg/rsg/r136_latest.txt and in viper vgpu/k1_{hg,rsg}.args (backups k1_*.args.bak_preza). The do-not-use mark is removed.
- **CUDA compile check:** v15 bb4735136 full build links, rc 0.
- **CUDA compile check:** v16 65f4321f3 full build links, rc 0.

### Evidence arm v16 + clean_stag_xdonor=lag + clean_stag_xsrc=true (lat16xd; for the USER's decision, not a key change in v16)
- **Better:**
  - flz q face FCT: dE/E 1.3e-3 -> 6.7e-16.
  - flz q BE only: 8.5e-9 / 50 CW / rc 1 -> 7.8e-16 / 0 CW / rc 0.
  - Bench CW 173 -> 49; GATES FAIL 20 -> 18.
- **Neutral:** WB, drag, mirror, sp3 and cs uni_gd rows stay within tolerance.
- **New:**
  - cs closed: sum E V rel 6.3e-10 (v16 2.2e-16), a rule-1 concern.
  - RT ringS: 4 CW (data same). RT flz q face FCT: 4 CW (data same).
- Table: stag2/integ/gate/tab_lat16xd_vs_Slat16.md.

### R136 A100 on v15 bb4735136 (Raven, 2026-11-19 evening, after the socket was reopened)
- **Half-deck RT** (r136_half_rt_sc + zalesak_active=true; pool 32490421 on GPUs 2,3; integ15_cu md5 b3068257): straight = forced dt_retry = restart 5->10, BITWISE (bins + hst), 0 CW. The straight run equals v13's straight bitwise; this deck has no face FCT, so za is inert here.
- **Production speed** (r136.00097 restart, 4 ranks x 1 A100, pool 32490420, the r136_latest keys incl. band_ring=false + za=true, 15 cycles):
  - **0.894 / 0.897 s/cycle steady (v8: 0.90)**, 0 CW, 0 FATAL.
  - The two runs give identical hst (md5 62df1e0d / c868c829), so A100 run-to-run reproducibility holds.
- Note: pool 421's GPU0 is used by a prod step; steps there must pin CUDA_VISIBLE_DEVICES (with --overlap a step does not get exclusive GPUs, so OOM otherwise).

### v17b 6bb6d0844: RE-BASELINE with the USER's decisions (2026-11-19 morning)
- **New keys** (in hg/rsg/r136_latest.txt, integ/keys/k1_*.args and viper vgpu/k1_{hg,rsg}.args):
  - rad_m1/clean_stag_xdonor=lag + rad_m1/clean_stag_xsrc=true: the P4 route; kpos stays off.
  - rad_m1/implicit_pos_floor_gas_keep=1 (A): fixes the flz rule-1 leaks.
  - rad_m1/implicit_pos_gas_cons=true (C).
  - rad_m1/clean_stag_mg_color=mirror: rule 10, mirror at round-off.
  - Code: HLLC guard 31051e6cb (B); GHOST layout_invariant CleanSums; XDFIX cs-closed seam fix.
- **Not adopted:** time2_face_fct_rounds=-4 (D, not bitwise, withholds energy) and the cs two-set ring fuse.
- **Rule rows:** see the v17b block in STATUS.md.
- **Re-baselined:** drag b0.1 R1/R3 t10 err +12 % from the HLLC guard (reference values 0.0025 -> 0.0028, 0.00052 -> 0.00059).
- **OPEN (owner XDFIX):** 35 counted CW (list in STATUS.md).
- **Checklist change (2026-11-19 allocation swap):** MI300A HG/RSG rows are not available after v17b. A100 multi-node HG/RSG rtchecks on Raven pHG262 32307710 (A100HR) replace them.
- **A100HR correction (2026-11-19):** every A100 HG/RSG rtcheck so far (v18) restarted at the dump cycle (leg b = 0 cycles + rst + 4). So it shows **rst round trip bitwise, not a mid-run restart**; the mid-run rows T4_hg_v18b / R12_rsg_v18b (fixed arm.sh: 1 cycle + rst + 4) both PASS BITWISE (v18 code d8ea9f4d4, LATEST keys). Absolute HG/RSG s/cycle from those arms were 1.25x too low. The R136 A100 RT and all CPU RT rows were real mid-run restarts.
- **CUDA compile check:** v17b 6bb6d0844 full build links, rc 0.
