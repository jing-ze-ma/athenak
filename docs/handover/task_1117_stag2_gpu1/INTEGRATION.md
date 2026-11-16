# clean-cc-stag2-gpu1 integration log

Base: 5e33f056a (origin/clean-cc-push16 head, viper). Worktree /orion/ptmp/jinma/clean_1115/int/wt.

## 1. Cherry-picks and conflict resolutions

### Rounds 1-17 (d6b6ed4e2..788e3391c)
- 8e47cdc55 (Newton change relative to |E|): **dropped**, viper's 0a3df762d does the same for both layouts (and
  ff74a049e's per-cell norm sits on top). Conflict resolved to HEAD; the commit became empty.
- 7ffb67be7 (extrap history in staggered restarts): Time2RstNchWant conflicted with viper's clean_rot_ring channel
  (b7964ecaa). Resolution: the extrap channel count without `cln_cell` (ours) and viper's ring line after it, same
  order as the pack/set code.
- 5c2ebac6f (fkin memory d in restarts): Time2RstNchWant/Pack/Set conflicted with viper's clean_rot_ring channels.
  Resolution: both, ring channels first (viper's order), then the 3 fkin d channels; the channel count follows the
  same order. Cell-layout restart files keep viper's layout (the d channels need clean_stag_fsolve + fkin).

### clean-cc-stag2-fix (f884c4fa5), applied after the rounds
- dropped as duplicates of viper's commits (identical `git patch-id`): 795db8814 = 57a08c09f (CleanGasEint),
  70374d536 = 4cacdae21 and 2af745af7 = 8db6f6141 (init audit).
- dropped as a pair: 76ca3548f and its revert 52d1025ce (net zero; d0986a72d is the kept area form).
- all other fix commits applied without conflict.

### clean-cc-stag2-wb (6dc54b2e3)
- 8b6979103 dropped (= 5e33f056a, identical patch-id); 6dc54b2e3 applied cleanly. The other wb commits are the
  fix branch's (already in).

### clean-cc-stag2 positivity (4934af164, c78474e98, 14c3f6a6c)
- applied without textual conflict.

### clean-cc-stag2-fix2 (2b73c9eae .. 83f251687)
- applied without textual conflict, BUT a latent slot clash: fix2's CleanDMem (196159442) kept the stage-start d
  in face slot CLF_N + 4, which gflux (c78474e98) uses as CLF_SLO (the low-order face mark). Each solve resets
  CLF_SLO = 0, so the saved d would have been wiped and a BE redo would have restored marks into d. Fixed at the
  dcdee83b1 resolution below (distinct slots).

### clean-cc-stag2-dp (dcdee83b1, 2ceddaef7, bb6bfb7ac, 8e5c4aaca)
- dcdee83b1 (vimp): conflicts in rad_m1.hpp (declarations, both kept), CleanMakeArgs (x.gf and x.drag both kept),
  CleanInit's face-slot count and rad_m1_clean.hpp's slot constants: CLF_SLO (gflux) and CLF_DB (vimp) were both
  CLF_N + 4, and CleanDMem used CLF_N + 4 as well. Resolution: named staggered slots
  CLF_DS = CLF_N + 4 (d save), CLF_SLO = CLF_N + 5, CLF_DB = CLF_N + 6, CLF_DQ = CLF_N + 7; nfc = CLF_DQ + 1 with
  vimp, CLF_SLO + 1 with gflux, CLF_DS + 1 with fsolve + fkin, CLF_N + 4 with fsolve (only allocation sizes change,
  no arithmetic). CleanFaceWeights: theta_d (vimp) and the gflux low-order flag combined; the low-order beta-term and
  advection forms take the drag-modified theta. CleanStagMemUp: the memory it takes upwind is now the one
  CleanFaceSource drops, (1 - w) theta [(1 - b) G + d] (vimp scales G by 1 - b); before, gflux + vimp (never run
  together on the sibling branches) would have taken the unscaled G upwind.
- bb6bfb7ac (upwind): CleanFaceWeights advection conflicted with gflux's low-order donor form. Resolution: with
  clean_stag_upwind the matrix takes upwind's donor cell 2 s a_up E_up on every face, and CleanStagAdvectCorr
  (the lagged rest) is skipped on a gflux low-order face (that face must stay pure donor-cell, an M-matrix with an
  E-independent part that cannot drain an empty donor); without upwind gflux's low-order form is unchanged.
- 8e5c4aaca (stag mg): CleanInit's multigrid levels conflicted with viper's clean_pc = prodmg (3543b578b), the same
  "x1 never halved, two coarse levels" hierarchy. Resolution: levels for (cell and clean_pc != x1line) or
  clean_stag_pc = mg, the prodmg shape for cell prodmg or stag mg, CLS_N components for the staggered mg. rad_m1.hpp
  members of both kept.

### clean-cc-stag2-perf (coordinator addition: the bitwise commits S2, S3, M2, M1, S6a)
The branch was rebased since the message (same patch-ids): 019d77dd4 = 4f627ad41, 7d2b29b56 = 51f652535,
5447412ad = c100efc38, f9eb668f0 = 3bbd41406, e22acf314 = b11fe322d. Picked by the listed hashes.
- 019d77dd4 (S2, F-solve line factors once): clean.
- 7d2b29b56 (S3, CLF_QL face weights once per pass): conflicts with the gflux/vimp slots and members.
  Resolution: CLF_QL = CLF_N + 4 .. + 6 (perf's), and the other staggered slots move up: CLF_DS = CLF_N + 7,
  CLF_SLO = + 8, CLF_DB = + 9, CLF_DQ = + 10. Staggered nfc = CLF_DQ + 1 (vimp), CLF_SLO + 1 (gflux),
  CLF_DS + 1 (fsolve + fkin), else CLF_QL + 3 (perf's CLF_N + 7); cell unchanged. CleanFaceCross with fwc reads
  the stored cross factor, and that factor (PC setup) now takes vimp's theta_d = theta + b (1 - theta), as the
  unstored form does; the gflux low-order return-0 stays first. Checked the order: CleanStagDrag(false) and the
  gflux marking run outside RHS -> PC setup -> BiCGStab, so the stored weights are those of the pass.
- 5447412ad (M2, one ghost scratch storage): empty here, viper's 98c731e30 already has it (skipped).
- f9eb668f0 (M1, no vet_gd band array with FV rays): skipped, identical patch-id to viper's 43eab024a (in base).
- e22acf314 (S6a, fused ghost write-back): conflicted with fix2's D_ab ghost parity (579faf103), which perf did
  not have. Resolution: the fused kernel multiplies the clamped copy by the product of the parities of every
  clamped side. Same values as the sequential x1, x2, x3 copies (each later copy reads the earlier one's ghost,
  so a corner gets s1 s2 a(edge)); multiplying by +-1 is exact. Also the reflect flags became an int bitmask
  (rf[6] was a runtime-indexed local array captured in a device lambda, an A100 NVVM hazard).

### 746e5dc98 cleanup (see section 2)

### clean-cc-stag2 7b74bf6db (beam gap fix, coordinator addition)
- applied cleanly. New closure slot CLN_JR = 49, CLN_N 49 -> 50. Checked: CLN_FH (39..42) and CLN_G (47, 48) end
  below it, the cell layout's CLC/CLF slots are separate arrays; cln_s, cln_mn/mc/mp, the ring arrays and the
  dt_retry visit of cln_s are all sized by CLN_N. Consequence: the restart channels of the extrap history
  (2 CLN_N + 1) and of clean_rot_ring (vra_k CLN_N + 1) grow by 2 and vra_k: restart files with those channels
  written by push 16 (5e33f056a) cannot be read by this branch (files written and read by it are fine; the reg
  RT tests cover that). Listed in the TASK.

### 2a85e23b8 cell-layout restart compatibility (coordinator requirement)
- CleanRstNq(): the cell layout writes/reads the first CLN_JR = 49 closure components per record of the extrap
  history and clean_rot_ring (push 16's format); staggered all CLN_N. Cross-version test: see the table (XV rows).

### Bug found by the first test round: every gas-coupled staggered run failed (fixed, 5da9b48fa -> 6a88a47ab)
- Symptom (round 1, binary 3cff6ccbb): rhse_sp/rhse_sp2/sp3/flz/rw/drag/shock all went E < 0 within a few cycles,
  dt collapse; Cartesian no-gas decks fine. Bisect by build (bis/): rounds, fix, positivity, fix2, dp heads good;
  dp + perf S2 (019d77dd4) bad, also with clean_stag_pc = line (the phys worker's "pc = mg" lead was a coincidence of
  its cases).
- Cause: an interaction of two sibling branches. Positivity 4934af164 stores the write-back's h div Phi in CLW_R0
  before the F solve; perf S2 stores the F solve's per-line right-hand side in CLW_R0 (and other E-solve slots).
  The write-back then took the F-solve data as E's residual. Fix: the write-back uses CLW_Q (unused by the F solve).
  After the fix: drag, sp3, flz q, rhse_sp2 run with 0 warnings (dbg/b2/c_*).

### Additional picks (coordinator)
- ba8612464 (perf S1, F-solve stop rule, CLEAN FSOLVE line): one conflict at the CleanSolve counters (cln_nclip,
  removed by the cleanup) -> only cln_fsn/cln_fsit/cln_fsres resets kept.
- a47fe0673 (M7, <rad_m1>/clean_stag_opac, opt-in): conflict in rad_m1.cpp's staggered key block -> the opac read
  after the clean_stag_pc read. Kept OPT-IN (USER/coordinator: worse with gflux on Marshak c/d); experimental.
- clean_pos_limit default off for staggered (USER decision; cell default true unchanged): a087a26c5.
- 25fb5b741 (gflux Marshak boundary momentum share): conflict with the cubed-sphere package's per-side cf/wf
  bookkeeping in CleanMomentumChange -> the new boundary share db goes to dm and to cf[side].

### Rebase onto 0678599eb (viper's final push-16 candidate) -> 72 commits replayed
- c91e5d859 (staggered fref on acceptance) vs viper 4ac0c0e49 (the same rule for the cell layout): one rule for both
  layouts (viper's, without `!cell`); our staggered t2_fail terms kept (cell part identical to viper's).
- 8a6ce1f2c (fix2 dt_retry checkpoint) vs viper d320347b3 (cell subset: cln_fvi, cln_fref_set): ours kept (the full
  cln_* block contains viper's two members; kDtrMinimal keeps cln_fc1..3 + cln_fvi).
- viper's clean_fv_joint opt-in taken as is.

### cubed sphere (coordinator: bf0cac59a, squashed cs package) -> ced70776d, clean pick.

### dt_retry visitor (crash found by the positivity worker, fixed 3e7f105dc)
- The cln_* block was regenerated with tools/dtretry_gen.py for the integrated members (841d51a76: + fsolve stop
  counters, clean_stag_opac, cs). Regenerating made it visit cln_hbuf (viper's shared exchange scratch), which the
  cln_halo exchange objects alias through unmanaged views: restored to the cycle-0 (pre-CleanInit) handle, the next
  exchange wrote freed memory ("double free or corruption" right after the first DT_RETRY, M7 Marshak decks).
  cln_hbuf is now skipped (and in the generator's SKIP_NAMES, with its reason). Reproduced and fixed on
  d_a_G.athinput (rc 134 -> the same rc/warnings as the positivity worker's own binary).

## 2. Cleanup (746e5dc98 and later)
- Removed study keys/diagnostics: clean_stag_bclip (limiter variants of the bracket, the cln_bfct kernels),
  clean_stag_tau_pos, clean_stag_pos_gas (cln_posgas/cln_poslim0), clean_diag_force (CleanStagForceDiag, CLEAN FORCE
  line), the CLEAN SPLIT line with its momentum-state counts (cln_nr3, cln_fh, cln_nclip, cln_eneg, cln_emax). The
  staggered "CLEAN cycle" line (r1/r2/first/solves) and the new CLEAN FSOLVE line (S1) are kept.
- Physics keys kept, all read only when named. Defaults for clean_layout = staggered = the recommended set:
  clean_stag_split true; clean_consistency true where the rays are not opt-in (not on the cubed sphere; on spherical
  polar only with vet_gd); clean_stag_fsolve true; clean_stag_fkin + clean_stag_fwork true (with fsolve and
  consistency); clean_stag_gflux true (with consistency); clean_stag_vimp true (with split; off on cs, notice);
  clean_stag_upwind true; clean_stag_pc = mg; clean_pos_limit = false (USER); clean_rays = fv1 on sp with vet_gd;
  clean_stag_opac off (opt-in, experimental). Cell-layout defaults and output headers unchanged (the staggered keys
  are read with DoesParameterExist; clean_pos_limit's GetOrAdd default is cln_cell, i.e. true for cell).
- Coupling notice (one line, only when a named key is dropped): "clean_stag_fkin and clean_stag_fwork need each
  other, clean_stag_fsolve and clean_consistency; clean_stag_gflux needs clean_consistency, clean_stag_vimp
  clean_stag_split: ignored". clean_stag_pc other than mg | line is FATAL.
- Lint: cpplint (repo CPPLINT.cfg) on every changed file, compared message-by-message with 0678599eb: 0 new
  errors except two long lines in the generated rad_m1_dtretry_visit.hpp (generator output, as the base's).
- GPU form: CleanGhostFill's reflect flags as an int bitmask (no runtime-indexed captured bool[6]).

## 3. Tests (final round: branch head 9689f4216, CPU binary bench/bin/athena_g1d md5 bd94f61e; base 0678599eb
build_base md5 8ef32cef; Orion gcc/13 openmpi/4.1 Release -ffp-contract=off NMAX 4096; arms concurrent per node)

Arms: CEL = cell layout at 0678599eb (bench, wb, rst) or 5e33f056a (phys: rw/drag/shock/fs, bitwise the same
for those decks per the reg harness); G1 = the recommended staggered set written out (incl. clean_pos_limit =
false); G1D = only `clean_layout = staggered` (the defaults); G1O = G1 + clean_stag_opac = 0.
s/cycle = AthenaK cpu time / cycles with all arms of a batch sharing one node (ratios indicative).

### 3.1 Keys off: reg harness (0678599eb-era copy, int/reg, arms N2 = G1 binary, O2 = 0678599eb)
| check | result |
|---|---|
| cmp.sh N2 O2 (cc ladders, cs incl. 2-rank thin1mbmr, flz q/Q, rhse, gtg, rw, rst runs) | 41/41 dirs `differ 0` |
| x1line variants (clean_pc = x1line: hcyl, cb, flz q, cs thin1mbmr 2-rank, sp3mbmr 2-rank RT) | 5/5 `differ 0` |
| RT (hcyl, sp3mb, feed, thin1, feedx1, sp3mbmr, feedmr) | 7/7 `hydro_u same m1 same` both arms |
| wb.sh fv1 / fv / lld / cart | identical: 1.21e-9 / 2.96e-9 / 4.75e-9 / 7.38e-9 |
| CLEAN WARNING, rc | 0, all rc 0 |

### 3.2 Staggered recommended set vs cell (bench, int/bench: tab_cart.md, tab_sp.md)
| deck | G1 score | CEL score | G1 s/cycle | CEL s/cycle | warn G1/CEL | min E G1 |
|---|---|---|---|---|---|---|
| cb x0.5 E_L1 (FWHM) | **0.46** (1.0 x4) | 0.46 | 0.178 | 1.26 | 0/0 | 0 |
| cbg x0.5 | **1.10** | 1.10 | 0.383 | 1.21 | 0/0 | 0 |
| ws66 x0.5 | **0.9209** | 0.9209 | 0.406 | 1.89 | 0/0 | 3e-9 |
| shd x0.5 (umbra 0.25) | 0.007134 (19.1) | 0.007139 (18.9) | 0.613 | 7.49 | 0/0 | 1e-4 (CEL -3.5e-2) |
| shdg x0.5 (umbra) | 0.02911 (0.595) | 0.02909 (0.565) | 2.05 | 8.93 | 0/0 | 1e-4 (CEL -2.5e-2) |
| hcyl x0.5 E_L1_out (F_L1) | 0.01265 (0.0211) | 0.01206 (0.0174) | 0.724 | 1.68 | 0/0 | 0.133 |
| nltek1 x1 L1 | 9.106e-5 | 9.120e-5 | 0.303 | 0.705 | 0/0 | |
| thin_r N32 L1(E)/L1(F) | 0.1176 / 0.1090 | 0.1179 / 0.1418 | 0.114 | 0.132 | 0/0 | |
| dif N32 | 0.01554 / 0.02471 | 0.01536 / 0.02406 | 0.124 | 0.158 | 0/0 | |
| lat N32 | 0.2914 / 0.2907 | 0.2915 / 0.2904 | 0.092 | 0.107 | 0/0 | |
| pulse0 N64 L1(E); closed-box energy | 5.487; **4.4e-16** | 6.016; 8.9e-16 | 0.974 | 0.956 | 0/0 | |
| rhse_sp max v/c_s | 0.1794 | 0.1851 | 0.033 | 0.056 | 0/0 | |
| rhse_sp2 WB max v/c_s | 1.2e-8 (pc = line: 2.0e-9) | 1.5e-9 | 0.0148 | 0.0172 | 0/0 | |
| sp3 1 vs 2 blocks: Er spread t=2 / m1 maxrel | 1.5e-9 / 6.5e-8 | 2.9e-10 / 7.8e-9 | 0.124 / 0.142 | 0.159 / 0.163 | 0/0 | |
| flz q t 32.88: cycles, dE/E, dLz/Lz, wall v_r | 27, 6.7e-16, 5.9e-16, +0.003/+0.010/+0.009 | 40, 1.0e-15, 5.9e-16, -0.003/-0.007/-0.007 | 0.854 | 1.33 | 0/0 | |
| flz Q t 80.39: dE/E, dLz/Lz | 5.9e-16, 1.1e-15 | 5.9e-16, 1.25e-15 | 0.368 | 0.439 | 0/0 | |
| wb.sh rhse at rest fv1/fv/lld/cart | 5.7e-9 / 1.8e-10 / 1.7e-9 / 4.2e-9 | 1.2e-9 / 3.0e-9 / 4.8e-9 / 7.4e-9 | | | 0/0 | |

G1D (defaults) vs G1 (keys written out): data-bitwise on cb, rhse_sp2, flz q, rw +-amp, drag x2, shock M2, fs A/B.

### 3.3 Gas-radiation physics (int/phys/score)
| family | G1 | CEL | verdict |
|---|---|---|---|
| rw 24 runs (P 1/3, 3.3; tau 0.1/10/1000; N64/N128; +-amp) | L1(E) = R16 full set to 4 digits (ratio 1.00 every case); orders 1.79-2.45; 0 warn | 1.76-2.45 | green |
| stiff drag 18 cases (R 1, 3; tau 1-100; icfl 0.4/4) | 0 warn / 0 retry, 1.1-1.9 solves/stage, dP, dE <= 4e-11 | 6 runs with 1 warn + retry | green |
| radshock M2/M5/M7 N256+N512 | L1 at or below CEL (M7 N512: Tg 1.9e-4 vs 1.4e-3); 0 warn; 5-7x cheaper/cycle | | green |
| radial free streaming (fs) | dm r^2 spread 4.1e-2, Phi_r r^2 spread 1.3e-2 | 6.1e-2, 7.2e-2 | green |
| M7 Marshak decks a-d (positivity worker's mw), chain off | warn/retry 4/5, 2/3, 0/0, 0/0, rc 0, no crash | 10/10, 5/6, 4/5, 3/4 | green (no crash after 3e7f105dc) |
| T-dependent relaxation eq (audit m7, kappa ~ T^-3, icfl 4) | 9 cycle-0 warnings + retries, rc 0 (64^2: 14) | 0 | known limit (c), see TASK |

### 3.4 Restarts, dt_retry, cross-version (int/rst)
| test | result |
|---|---|
| RT staggered hcyl 2x4, sp3 2 blocks, flz q (16 vs 8 + restart + 8) | `same same`, 0 warn |
| RT cell hcyl, sp3mb, feed, thin1 (G1 binary) | `same same` |
| RT with clean_rot_ring (live) | DIFF 3e-3 for the cell layout, identical with the base 0678599eb (pre-existing, viper); staggered DIFF 1e-1 -> known limit |
| cross-version: 0678599eb writes cycle 8, both binaries continue (sp3mb, feed, feedx1, thin1, ring) | 5/5 `same same` (cell restart format kept) |
| forced dt_retry at cycle 4, fac 1 (staggered q, sp3, hcyl; cell q, sp3) | data-bitwise the straight run (5/5) |
| forced dt_retry fac 0.5 | 1 retry, runs on, 0 warn (5/5) |

### 3.5 clean_stag_opac (M7) cost (G1O vs G1)
Constant-opacity decks: bitwise no-op (sp3, flz q, rhse_sp2, cb, hcyl, nltek1), s/cycle within the node noise
(+-8 %). T-dependent eq deck: fails at cycle 0 (dt_retry exhausted) where G1 recovers. -> stays opt-in.

### 3.6 GPU readiness
- Review of every new device kernel (A100 NVVM classes, HIP team scratch, implicit this): no must-fix. One
  should-fix: cls_galerkin's local am[3]/ap[3] indexed by the trip-3 loop d (rad_m1_clean_solver.cpp ~481-507);
  notes: cls_restrict pos[3] (same form as the accepted cell code), CleanMomentumChange's Real *mag / cf[side] /
  wf[side] pointer to locals (cs package, runtime side), cln_fwc writes al/ar through View references, shadowed
  names in cln_fs2f. No team scratch, no TeamPolicy, no atomics, no member use in device lambdas.
- CUDA compile (cuda/12.1, sm_80, nvcc_wrapper, --fmad=false): see section 4 (status at hand-off).

## 4. CUDA compile and pushes
- CUDA (cuda/12.1 + nvcc_wrapper, sm_80, Release, MPI, --fmad=false, NMAX 4096; cuda_tests/build_cuda.sh;
  build dir int/build_cuda): rad_m1_clean_solver/_cell/_closure/time2 and the rest of rad_m1 compiled, 0 errors
  (tree 841d51a76). rad_m1_clean.cpp still in cicc after 45 min at hand-off (reported in the TASK).
- Verdict of section 3: all GREEN (known limits listed in the TASK; none is a regression of the cell layout).
- Pushed: branch clean-cc-stag2-gpu1 = 9689f4216 to origin (fork). TASK on bsg-files-1009 (hash in the report).
