## ===== RESUME HERE (10-21 ~05:00, coordinator session 7ad91e13; CLAUDE ACCOUNT SWITCH next) =====
Read memory first (MEMORY.md index; esp. handover-2026-10-21, globals-stay-4de8168d, fast-tests-bundle, build-incremental-cpu-node,
no-rerun-known-arms, gate-current-prod-decks, idort-page-paper-only, tau-free-blend-study, be-spurious-growth-0926, fmode-lag-wb-cache).
Workers of this session cannot be messaged by the next one: read their dirs/NOTE files.
### Code
- fin51-r-fix = 077dda23 (PUSHED, = prodmerge-1019b): 4de8168d + pgrst-merge-1021 (pgen registers restored on the hesdirk2->BE redo,
  Finalize frees (rc 134 shutdown abort fixed), accretor booking A, opt-in time2_k1_reset(+_grow)). 4de8168d = 1cb24f82 (memcut + opt-in
  implicit_cs2_clip) + cbg/cbh BE-fallback fix. TRAIN.md (amdef_1019) has the gate history.
- Every train tip also builds a plain A100 bin (built-in pgens, IDORT on GPU). Default gate set = CURRENT prod decks (RSG E clip N32,
  HG global clip N32, hg145, accretor v12lzc), not bsg/e1/hg.
### Production (USER: globals stay on 4de8168d + hesdirk2 + conservative2; only the accretor is on 077dda23)
| chain | bin | deck | state |
| RSG E R2 Raven 31316848 (+4) | he_a100 4de8168d 9af6efbb | rsg_E_start1020r_clip_rst5d 6d28a2d0 | pending, est 10-22 ~04:47 |
| HG global R Raven 31337226 (+4) | he_a100 4de8168d | hgglob_fresh_start1020_raven_clip 59ff97b9 | pending, est 10-22 ~05:30 |
| HG wedge hgP145 31263626/27 | 82275582 (old defaults, USER) | hg145_prod | pending, est 10-21 ~21:30; not held |
| Accretor run12 31291344 (+3) | acc nofma 077dda23 7e7ea3d8 | acc3d_prod_v12lzc_be_gm c6af572f (v12lzc + time_scheme be + implicit_gas_src matched) | released, FRESH t=0, pending |
| Viper V2 12291385 / HG V 12278622 | he_apu 4de8168d_r71 0a0d8746 | clip decks (V2 rst 5 d) | pending, no estimate |
| DeltaAI RSG E | its 4de8168d build c56693cb | clip, rst 5 d (44b2ef21), output on /projects/bivj | RUNNING (started 10-20 18:56 CDT, 2.25-3 s/cycle, dt 1e3, ~15 sim d / wall h, 52.9 GB/GH200) |
| Caltech RSG E 5126654 | its 4de8168d build 9a56b4a1 | clip_rst5d a7e734c4 | pending |
| Caltech accretor | 077dda23 (Caltech to build) | template docs/handover/acc3d-prod-1018/acc3d_prod_v12lzc_be_gm.athinput.in 5de0635e | V12_READY NOTE cdf71b98; FRESH; races run12 |
- Notes on bsg-files-1009 (wt-bsg): RSG-E rst5d template, DeltaAI disk plan (rst 5 d, keep 2, /projects), T_gas tool (rsg-page/tgas.py +
  eos_dump_rsgE.athinput, NOTE 224bb271), rsg-page pipeline (c542d6b0), CORRECTION 14db112d / 6052161e: RSG stays on 4de8168d.
- OPEN USER questions: (1) DeltaAI dumps every 5 d vs daily + thinning (disk fills in ~1 d at ~35 GB/h); (2) Viper chains 6 vs 12 nodes
  (memcut_1020/PLAN.md: 2 blocks/APU on 6 nodes, same cost; meshblock changes only at a fresh start); (3) time2_k1_reset / new gas-source
  key default (later).
- Prod log: /raven/ptmp/jinma/prod_monitor_1020.md. Swap tooling: /raven/ptmp/jinma/prod_swap_1021 (swap_all.sh has a keep-deck mode).
### Accretor method (USER: eventually hesdirk2 + conservative2; for now BE + matched)
- GO screen BE+matched = v11 speed (1.00x/1.02x), 0 redos, budgets closed, fast floor cells = v11 (6.6k/8.2k).
- conservative2 kicks floor gas (73k fast cells, max 1028-1175 km/s). NOT held light and NOT E/(rho c^2) (max 4.9e-3): 99.5 % of the fast
  cells are clip-limited free-streaming floor cells (|F| = cE) at the gainer surface / theta walls (r 9-9.6); conservative2 hands the change
  of the stored F (limiter clip, holds) to the gas in the next solve. Fix in progress: opt-in <rad_m1>/implicit_cs2_dacap = K (cap dA at
  2K|w|max(E); K=2 strict) — worker gsdrag, /raven/ptmp/jinma/gsdrag_1021 (NOTE.md, branch gsdrag-1021 off 077dda23; IDORT job gsdidort
  31353150 running). Then screen v12lzc hesdirk2 + conservative2 + dacap (+ vimp hook): 0 SS redos, fast <= 1.5x v11, budgets, cost <= 1.6x.
- hesdirk2 stream failures: stage-1 inadmissible every stream cycle; k1_reset fixed the static phase (68 -> 1); remaining cells beside the
  nozzle dominated by vimp_JRHS (vimp off: 106 -> 62). Worker pgrst, /raven/ptmp/jinma/pgrst_1021 (+ acc3d_1017/acctrc/NOTE.md); hook the
  vimp work term to the same dacap / |F|=cE criterion, coordinate via gsdrag_1021/NOTE.md 'vimp' section.
- Accretor energy history: booking A merged (boundary radiation per solve); B / K1 booking NOT merged (over-corrects; pgrst-1021 only).
### IDORT page
- https://claude.ai/artifact/DvXFF6iQyjBXMru4NPi6WV = PAPER TESTS ONLY (builder idort1020/page/mkcontent.py filter + build19.py).
- (P, tau) linear-wave grid RUNNING at 05:05 (job rwgrid 31352365, 2 h limit, idort1020/rw cases_grid1021.txt, bin 4de8168d CPU): NOT yet added to the page — run
  page/rwconv.py then page/rwgrid.py (writes fig/conv_rw_grid.png + conv_rw_ordermap.png), add both to section s432, rebuild, republish.
  CHECK FIRST: P 33.3 tau 0.1 N256 has l1 = 1.5e4 (measured omega 36.1 + 33.8i vs ref 36.0 - 0.52i): a blow-up or wrong root.
### TODO 10-22 (USER)
- tau-free central/half-range blend criterion study (memory tau-free-blend-study): idort_f vs idort_a vs f vs new H/J + Eddington-tensor
  anisotropy variant on IDORT paper tests + the order grid.
### Artifacts (account switch)
- Backup of the 6 pages updated 10-19..10-21: /raven/ptmp/jinma/artifacts_backup_1021 (INDEX.md, sha256.txt); older backups _1019, _1014.
  New account: republish from these (strip the service skeleton first line + trailing </body></html>; DvXF is clean), record new URLs in
  memory page-urls.md. Claude Docs page 'Shock-Supported RSG Atmospheres' (4QPzw) needs a manual export from its page menu.
- 10-21 ~06:30 radwave P_rad/P_gas 33 tau 0.1 blow-up (my grid extension, c_light 10): ROOT CAUSE = conservative2 write-back 3-pass loop under
  time2_vstage (rad_m1_impl_wb_kernel.hpp ~525-576), gain om(1+f_xx)E/(rho c chat) = ~1.3 at E/(rho c^2)=1. USER: not a Newtonian regime
  (prod E/(rho c^2) <= 4.9e-3, chat = c) -> fix implicit_cs2_wbsolve PARKED on pgrst-vimp-1021 @ 455a5265 (opt-in, unmerged, untested);
  caveat: would matter with RSLA chat << c. On the IDORT page mark P_rad/P_gas 33 as outside the test's Newtonian validity (E/rho c^2 = 1).
- 10-21 06:15 IDORT page v7: (P,tau) linear-wave grid added (P 33 omitted, E/(rho c^2)=1). FINDING: thin waves (tau<=0.1) E order 0.66-1.6 -> first item for the 10-22 blend study.
- 10-21 ~06:45 USER went to bed ("leave you to it"). Night plan: workers continue accretor hesdirk2+cs2 diagnosis (pgrst: combined branch
  pgrst-comb-1021 d9b4ef14 = k1_reset=pos + vimp_ewf + dacap, stream arms sA-sD job 31355869 + gpudev; gsdrag: bisect radiation side of
  conservative2 vs matched for the r 9.1-9.5 equatorial v_phi spin-up). NO production change without the USER (accretor stays BE+matched,
  globals stay 4de8168d). Coordinator watches prod starts (HG wedge est ~14:00, R2 ~19:00, HG R / run12 ~22:00) and site NOTEs.
- 10-21 ~07:00 gsdrag FOUND residual spin-up driver: drag on an unrealisable cell flux (|F_iter| > cE). matched is clean only via implicit_realisable_coupling=true (drag x min(1, cE/|F_iter|)); matched with it false = cs2+K2 level (18.9k). Fix g1f (dacap on: cs2 write-back gas charge capped at the realisable matched drag + rad-inertia excess <= 2K|v|E/(c chat)) in gpudev job 31356947.
- 10-21 ~07:30 gsdrag g1f (24810dad): BE + conservative2 + dacap K2 = GO level for the floor gas: SS 9.07k fast (1.11x v11), budgets = matched, cost = matched. hesdirk2 K2: 10.5k fast but still 199 NOT ADMISSIBLE / 15 dt_retry (pgrst domain, open). Gates on g1f running (31357315/16/17, IDORT 31357323). USER decision later: accretor BE+matched -> BE+conservative2+dacap? (one gas-source method) — not applied.
- 10-21 ~08:00 pgrst: k1_reset=pos acts on ~6 cells/cycle but changes nothing on the accretor stream (103/103 cycles BE); remaining failures likely eint_writeback <= 0 (solve's kick work > floor-gas eint). Next (coordinator GO, opt-in only): confirm with dbg_t2_admiss; check why existing implicit_pos_gas limiter doesn't prevent it (ordering vs admissibility test?); else opt-in energy-conserving kick-work cap. Accretor prod stays BE+matched.
- 10-21 ~08:30 pgrst CORRECTION: accretor stream stage-1 failures are the vimp operator's own positivity fallback (E_solved <= 0, min -33, back-sink floor cells beside the nozzle; ~2/cycle; admissibility rejects via vimp_now != impl_vimp), NOT gas eint. Admissibility (~15300) checks the raw solve before write-back limiters. Next test: implicit_vimp_ewf extended to the vimp Jacobian (pgrst-comb-1021 bbd6b430; arms sAe/sDe/sEe, interactive 31357631 + gpudev pgg14).
- 10-21 ~09:00 pgrst final report (NOTE.md in pgrst_1021): stream hesdirk2 still fails every cycle after dacap, vimp_ewf (+Jacobian bbd6b430), k1_reset=pos, vstage off. Remaining: vimp positivity fallback beside nozzle (~2/cycle) and stale gas-energy K1 at the stream front. Coordinator: worker continues with a-posteriori local re-solve (T2POS time2_posfix ace7cd39 ported onto pgrst-comb-1021), opt-in.
- 10-21 ~09:30 gsdrag FINAL (NOTE.md in gsdrag_1021): fix 24810dad (local only) = opt-in implicit_cs2_dacap=K incl. realisable matched drag reference. BE+cs2+K2: 9.07k fast (1.11x v11) GO, budgets/cost = matched. Gates: key-off bitwise RSG/HG/hg145; key-on RSG, hg145 bitwise; HG global changes 1e-3 by design (realisable drag in its free-streaming atmosphere); accretor det/rst/retry bitwise; IDORT PASS (cloud fields 1.2e-2, row measures to 3 decimals). USER DECISIONS PENDING: (1) accretor deck matched -> conservative2 + dacap=2 (BE)? (2) dacap in HG global? hesdirk2 still needs the pgrst redo fix.
- 10-21 ~10:30 pgrst: implicit_vimp_kill (local vimp drop, a1c5d379) partly helps (sAk 15 stage / 130 BE; 132 remaining rejections converged with S2<=0). time2_posfix ported (pgrst-posfix-1021 f542399b, pgen registers saved/restored around re-solve); build+cbg/cbh on interactive 31360181, stream arms on gpudev pgg16 (sAkp, sAp, sAekp, sGkp).
- 10-21 ~11:30 pgrst: posfix (e1612750) acts but does not cure cbg/cbh (cbg worse: 28 NOT ADMISSIBLE, Picard divergence). Insight: hesdirk2 REJECTS stages with an operator positivity fallback (vimp/offdiag dropped) while the BE redo ACCEPTS the same fallback -> rejection only discards the 2nd-order stage. New opt-in time2_accept_fallback (5b525abb): accept converged stage with fallback if E,T,eint > 0 (counted). Tests: cbg/cbh interactive 31366882; accretor gpudev pgg18 (sAf, sAfp, sHf, sIf). USER should judge: accepted stages run with the operator term dropped.
- 10-21 ~15:30 NEW ACCOUNT. Pages republished (memory page-urls). Caltech: RSG E link1 IB fault (requeued), rdiff H200 NOT bitwise (R!=C, R2!=R) -> TASK 1ba37054 (UCX_RNDV_THRESH=inf / UCX_IB_GPU_DIRECT_RDMA=no). Accretor hesdirk2: best patch arm 92 stage/82 BE; pgrst resumed on new direction: pgen holds as Dirichlet rows in the implicit stage solve (opt-in), + gsdrag dacap 24810dad.
- 10-21 ~15:45 USER: queue HG global on Caltech -> TASK 1d352abb: 4de8168d H200 build, template docs/handover/hgglob-1020/hgglob_fresh_clip.athinput.in (= Viper HG V deck + @FILES@), 2x4 H200 (3 blocks/GPU), hold link 1 until the UCX test (TASK 1ba37054) picks the env.
- 10-21 ~16:00 pgrst: implicit_hold (pgen holds as Dirichlet rows at the 12 EFIX sites; ry_per_accretor marks nozzle/neck/donor lobe) implemented on pgrst-posfix-1021 2c53b9a2, opt-in; booking of held energy via face fluxes TODO; build interactive 31372635 -> gpudev pgg19 arms h1-h4. cbg + accept_fallback: 23 NOT ADMISSIBLE, 0 accepted (genuine S2<=0).
