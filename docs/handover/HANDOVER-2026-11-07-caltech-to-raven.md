# HANDOVER 2026-11-07 Caltech -> Viper: continue RSG E, HG, R136a1, WASP-121b and BSG tests on Raven

**USER 11-07 ~09:35 PST:** "prepare a complete handover to viper to run on raven. we should try to continue rsg, hg, r136, wasp121b, and tests for bsg there." Also: "once those pass we should make them default" (= the r7 candidate keys and the full Eddington tensor on cs, after their open gates; section 4).

Why: every Caltech GPU allocation is pending on priority, estimated start ~11-17. Nothing runs on Caltech GPUs now. All five runs are stopped at a clean rst.

**Nothing may be submitted on Raven until the USER says go.** The rst copy is user-driven (MPCDF password + OTP, section 6). Order per run: files arrive, md5 check, build, smoke, `compare_hst`, user go, chain.

J = `/resnick/groups/carnegie_poc/jingze` (Caltech). Support files (exact command lines, cfgs, file list, R136a1 patch) are in `docs/handover/caltech-to-raven-1107/` on this branch.

Earlier notes this builds on: NOTE-2026-11-05c (your Raven run dirs `/raven/ptmp/jinma/raven_gap_1105`, binary `athena_he_a100_nofma_23a55658_fmadoff` md5 620c8552, scripts), NOTE-2026-11-07-viper-caltech-raven-readiness (queue ~42 h, run.cfg out of date, no rst, ~65 % of MPA usage for 2 x 6 nodes x 4 days). This note answers the "please push" list in that note: runG.cfg, runH.cfg, last commands, binary hash and recipe, rst md5s.

---

## 0. Summary table

| Run | Newest rst on Caltech | Size (bytes) | md5 | Sim time | Production code (binary -> commit) | Caltech layout | Raven plan |
|---|---|---|---|---|---|---|---|
| RSG E | `$J/rsg_E_1020/prodH/rst/rsg.01295.rst` | 45,360,586,583 | 6a47bdb7f91fe8875932f5877a8443f0 | 5.4648e8 s = 6325 d, cycle 589788 | `athena_gpu_he_star_m1_prod_23a55658_fmadoff` (md5 38231ce9f0e39a4672cc74400f18b6e4) = prod-1102 23a55658 | 24 blocks, 8 ranks on 8 H200 (3 blocks/GPU), ~57 GB/GPU | 6 nodes, 24 ranks x 1 block |
| HG (He giant) | `$J/hgglob_1021/prodG/rst/hegiant.00139.rst` | 97,543,425,853 | 6109201a71e5cad1d123116733218db5 | 5.0112e7 s = 580 d, cycle 624961 | same binary as RSG (23a55658) | 24 blocks, 12 ranks on 12 H200 (2 blocks/GPU), ~49 GB/GPU | 6 nodes, 24 ranks x 1 block |
| R136a1 deep B | `$J/r136a1_1104/runB/prod/rst/r136.00008.rst` | 12,078,201,051 | 3a729cf95e2a67ab5f46739ebb10b2f3 | 648608.44 s = 7.507 d, cycle 147106 | `athena_gpu_he_star_m1_prod_dc2dc20c_fmadoff` (md5 d483f8af3bb04c08a8779e0a0cd412e2) = prod-1102 dc2dc20c | 4 blocks (sp, 1x2x2), 2 ranks on 2 H200 | 1 node, 4 ranks x 1 block (memory to be measured) |
| WASP-121b r1m30 (MHD, 30 G) | `$J/w121r1_1101/runs/r1m30/rst/dhj.00118.rst` | 339,657,930 | 9bfda8dcf008846c679b7b4f71ed275b | 1.2998172e7 s = 150.4 d = 118 rot, cycle 3381468 | `athena_gpu_deep_hot_jupiter_rt_wallwb_0f1bbb708_fmadoff` (md5 9da0418b115765acb2a5979c532e0ee2) = prod-1102 0f1bbb708 | 6 blocks (cs C32, 1 per panel), 2 ranks on 2 H200, ~18 GB/GPU | 1 node, 2 ranks x 3 blocks (as Caltech) |
| WASP-121b p1 (hydro, test state) | `$J/w121r1_1101/runs/p1/rst/dhj.00302.rst` | 294,022,648 | c0294b764505aae982e10a715e57bc7d | 3.3266514e7 s = 302 rot (finished, 300-rot goal) | test state for the multi-D WB A/B only | 6 blocks, 1-2 GPUs | 1 node |
| BSG prod_lhxo (tests only) | `$J/doface_1105/prod_lhxo/rst/bsg3d.00058.rst` | 6,249,790,621 | f91401ca699c76b100086dc16975b33d | 4.95936e6 s = 57.4 d (finished at tlim) | `athena_gpu_he_star_m1_v3_fmadoff` (md5 20efaed08d5ab3789073112597af3ed9) = vstarfix-1106 3b1d3dc19 (UNPUSHED) | 16 blocks (sp 256^3, 1x4x4), 4 ranks on 4 H200 | 1 node, 4 ranks x 4 blocks |

Small test states (for 1-GPU gates; the full HG/RSG states need 6 Raven nodes even for a test):

| State | Path | Size | md5 |
|---|---|---|---|
| HG small ("hgb", C32, 704x32x32 x 6, used by every HG gate on Caltech) | `$J/prod_1103b/he/runs/hgb/rst/hegiant.00001.rst` | 6,670,021,763 | 596fe2d4525248df08f75f8bf84036aa |
| RSG small ("rsb") | `$J/prod_1103b/he/runs/rsb/rst/rsg.00001.rst` | 3,968,947,892 | 4c1c93d040472e3b26e9f0e03e1f52e1 |
| R136a1 rst 5 (start of all r136open/r136stage arms, 5.000 d) | `$J/r136a1_1104/runB/prod/rst/r136.00005.rst` | 12,734,409,174 | 86295c8bd33b8b354e90d8a7f739b011 |

Dependencies bundle (decks, IC/opacity tables, hst files, launch scripts, cfgs; relative paths from J): `$J/raven_copy_1107/raven_deps_1107.tar`, 1,167,595,520 bytes, md5 426a2d8e106c022f3f97fbf36cb0c276. File list: `caltech-to-raven-1107/deps.lst`. Unpack it under one root on Raven (e.g. `/raven/ptmp/jinma/caltech_1107/`) so the tree matches J.

**Absolute paths stored in the rst parameter blocks.** The he_star rst files carry `problem/he_ic_file`, `he_opac_table`, `he_planck_table` (BSG also `he_ic_linhao_file`) as absolute Caltech paths. Name each on the Raven command line with the unpacked copy (same file, same md5) unless you know the pgen does not re-read them on restart. dhj (r1m30, p1) uses relative `ck/` and `ic/` paths: make the run-dir symlinks `ck -> pkg/ck`, `ic -> pkg/ic` as on Caltech, then check `pkg/MD5SUMS`.

| rst | keys to repoint (Caltech value) |
|---|---|
| HG (also hgb) | he_ic_file `hgglob_1021/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt`; he_opac_table `.../rosseland_mesakap_afe2dd8da2bd.txt`; he_planck_table `.../planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt` |
| RSG (also rsb) | he_ic_file `rsg_E_1020/files/rsg-E-1020/ic_rsg_g25_deep.txt`; he_opac_table `rsg_E_1020/files/rsg-g25-prod/rosseland_mesakap_4c0e8b64d26a.txt`; he_planck_table `.../planck_ext2_x0.653848_z0.019757_d38e5526afef.txt` |
| R136a1 | he_ic_file `r136a1_1104/ic/ic_r136a1_B.txt`; he_opac_table `r136a1_1104/opac/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt`; he_planck_table `r136a1_1104/opac/planck_tops_x0.58_z0.0064_a09_T7.txt` |
| BSG | he_ic_linhao_file `bsgli_1014/files/linhao_profile_full.txt`; he_ic_file `bsgli_1014/files/ic_ma2026_f0.txt`; he_opac_table `bsgli_1014/files/rosseland_ma2026_x0.7_z0.008_padded.txt`; he_planck_table `bsgli_1014/files/planck_ma2026_x0.7_z0.008_padded.txt` |

The EOS is the in-code table (`general_eos = table`, built at start-up); no file.

---

## 1. Code

### 1.1 Production commits (all in `origin/prod-1102`, pushed, tip 6b744f88b)
- RSG and HG: **23a55658**. Your Raven binary `athena_he_a100_nofma_23a55658_fmadoff` (md5 620c8552) matches this source. Use it for the exact continuation.
- R136a1: **dc2dc20c** (prod-1102 ancestor). It has `rad_m1/time2_adm_poslim` (production ran with it since 11-06 16:47). No Raven build yet.
- r1m30: **0f1bbb708** ("wallwb v2", prod-1102 round 24). No dhj A100 build on Raven yet.
- BSG prod_lhxo: **vstarfix-1106 3b1d3dc19**, NOT pushed (latke-1106 points at the same commit).

Caltech build recipe (all `_fmadoff` binaries): CUDA Release, MPI on, `-ffp-contract=off --fmad=false`, HOPPER90, Kokkos d8e9af03, gcc 13.2 + cuda 12.9 + hpcx 2.17.1. For A100 use your 11-05 recipe: AMPERE80, gcc/13 cuda/12.6 openmpi_gpu/5.0, `Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF`. PROBLEM = he_star_m1 (RSG, HG, R136a1, BSG), deep_hot_jupiter_rt (W121), ry_per_accretor (accretor gates).

### 1.2 Merge branch prod-1107 (being prepared on Caltech now)
`prod-1107` = prod-1102 6b744f88b + merges of layoutfix-1106, dhjwball-1106, rsgdiag-1107, seamE-1106, cs-tan-1106 (via cs-tan-seam-1106), rstbit-1107, + 4 default flips (USER 11-07 "make some defaults"):
- `<rad_m1>/rst_state_ext` default true (hesdirk2 bitwise restarts; false = old record);
- `<mesh>/cs_seam_pos_global` default true (layout-independent seam weights);
- `problem/ck_impl_global_ctl` default true (dhj ck pass control reduced over all ranks);
- `problem/wb_lobe_cut` default true in deep_hot_jupiter_rt and ry_per_accretor (no WB of any kind outside the L1 equipotential).

**Status at 09:50 PST 11-07: local tip c41b749ca (60 commits on prod-1102), NOT pushed to origin.** The coordinator will push it and send a short NOTE. Until then none of the branches below are on origin except prod-1102 and meshfb-1106:

| branch | local tip | on origin? | content |
|---|---|---|---|
| prod-1107 | c41b749ca | no (coming) | merge + 4 default flips |
| rstbit-1107 | 4bdd4e82f | no | rad_m1/rst_state_ext (+ dbg_rst_fp) |
| layoutfix-1106 | 383b4fd56 | no | cs_seam_pos_global, ic_rad_precond=line, cs_precond_line, ck_impl_global_ctl |
| dhjwball-1106 | be0475baf | no | hydro/wb_dirs=all on dhj, problem/wb_lobe_cut, probwb face-potential fix c26c8e7f4 |
| rsgdiag-1107 | 57aa6254f | no | implicit_predictor_ecap (246e5e400), implicit_tstart_old, CPU diag aliasing fix 89d801876 |
| seamE-1106 | 7de51b303 | no | cs_cross_mono, guard bit 8, implicit_predictor_tcap, time2_k1_redist, time2_k1_fct, seam dump (based on meshfb-1106) |
| meshfb-1106 | 8f7aa1d17 | yes | implicit_fb_local (+ rounds, fb_diag) |
| cs-tan-1106 / cs-tan-seam-1106 | 4fd91eb7d / f5d002cad | no | rad_m1/vet_gd_tangential on cs (full Eddington tensor) |
| faceextrap-1106 | 53b63aba3 | no | face opacity / extrap2 test keys (not merged) |
| vstarfix-1106 (= latke-1106) | 3b1d3dc19 | no | BSG binary v3 |
| athpp-doface-1105 | bef89c7bf | no | rad_m1/doppler_face (BSG, not merged) |

After prod-1107 is merged, every merged key except the 4 flips is opt-in (read only when named; keys-absent bitwise vs prod-1102 on GPU per each branch's gate). Not on any branch: the R136a1 diagnostic `problem/he_wb_src_linear_rmin`, saved as `caltech-to-raven-1107/r136_srclin_he_star_m1.patch` (diff against dc2dc20c).

### 1.3 Which binary for the continuations (recommendation)
- **First Raven link of RSG, HG, R136a1: the exact Caltech commit** (23a55658 / dc2dc20c), Caltech keys verbatim. Then `compare_hst` against the Caltech hst can show ~1e-12..1e-10 agreement. The prod-1107 flips change results in the last bits (cs_seam_pos_global), so the smoke compare would no longer be clean.
- Move to prod-1107 (+ r7 keys / full tensor) only at a later restart, after section 4 gates pass and the USER approves. Name every new key on that restart command line.
- r1m30: see section 2.4. It should not simply resume on 0f1bbb708.

---

## 2. Per-run sections

Common to all runs:
- **Stop rules** (every launcher; Caltech pattern in `r136a1_1104/runB/resume_poslim.sh` `bad()`): STOP on any `FATAL`; on the run-summary `NON-CONVERGED=[1-9]` only (never the per-pass "NON-CONVERGED after" lines that a DT_RETRY recovers; user rule after the 10-16 false stop); on `nan` in the last 60 log lines; on rc != 0. rc 139 (segfault): auto-restart from the newest rst, at most 2 times.
- **Watcher**: one background poll every 10 min (user rule). Report: link ended, STOP, rc != 0, and the HG-type stall: dt < 1 s on 2 consecutive checks. No rolling watches beyond that.
- **Restart files**: keep the newest 3 per run. Never prune hardlinked `rst_keep` files (HG every 10 d, RSG every 25 d of sim time).
- **Kill rule**: stop a run only by exact PID, after reading the full command line (run dir, rst, GPUs). Never by pattern: on 11-06 a pattern kill took R136a1 down with a test. A plain scancel loses work back to the last rst. There is no signal handler.
- **GPU env** (Caltech): `OMP_NUM_THREADS=1 OMPI_MCA_io=romio321 UCX_TLS=^xpmem` (+ `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1`, harmless on CUDA). The he_star 1-node launchers also set `OMPI_MCA_coll_hcoll_enable=0` (11-04 hcoll SEGV in ImplicitGFBuild's MPI_Allreduce). Use your validated Raven OMPI/UCX set.

### 2.1 RSG E (`rsg_E_1020/prodH`)
**Physics.** Global 3-D RSG envelope ("RSG E", g25 model, R_ph ~1018 Rsun), full sphere, cubed sphere C128 (6 x 2x2 blocks of 548 x 64 x 64; radial geometric head below old face 110). Implicit M1 + VET (vet_gd, short-characteristics closure) radiation with the general EOS table and MESA opacities, hesdirk2 time scheme, cfl 0.3, dt ~920 s (cap 1000 s). Started fresh 10-25; prodH restarted from 1675 d with the vertex fix on 10-29. Purpose: RSG convection, pulsation and mass loss. Science found so far: F mode (~455 d) dominates since ~2570 d over a 1O mode (~215-245 d) (rsgmode_1106); R_ph p2p ~0.1 R0.
**Now:** 6325 d (rst 1295, written 11-06 18:24 PST at STOP for tests).
**Last restart command** (11-03 22:22 PST, from rst 958, ran to rst 1295): `caltech-to-raven-1107/rsg_last_cmd.txt`. **Next command from rst 1295** (output last_time values from `bsg_1009/files/rst_info.py`): `caltech-to-raven-1107/rsg_next_cmd.txt`. Keys = runH.cfg XKEYS: `rad_m1/time2_vet_col=extrap2 rad_m1/time2_accept_fallback=true rad_m1/time2_posfix=4 time/dt_retry_t2_be=0 rad_m1/implicit_work_fix=true rad_m1/cs_vertex_mono=true rad_m1/implicit_pos_floor_cons=true rad_m1/time2_auto_be=16 rad_m1/vet_gd_twin=false rad_m1/vet_gd_symref=none rad_m1/vet_gd_twin_fuse=false rad_m1/vet_gd_twin_state=mean rad_m1/vet_source_noesrc=false time/dt_retry_host=false rad_m1/vet_rotate_advance=true mesh/cs_seam_blockend=extend hydro/merge_lat_blockend=exchange mesh/cs_vertex_3cell=true mesh/cs_seam_rho_guard=0`. Full cfg with its history: `caltech-to-raven-1107/rsg_runH.cfg`. Deck (fresh only): `rsg_E_1020/rsg_E_1025_fresh.athinput`. Your raven_gap_1105 run.cfg is from 11-04: diff it against rsg_next_cmd.txt.
**Ranks/GPUs:** Caltech 8 ranks on 8 H200 (2 nodes, `best2n_wrapF.sh`, 1 GPU per rank, PIN=1, `UCX_IB_GPU_DIRECT_RDMA=yes`, srun `-c 8 --mpi=pmix --kill-on-bad-exit=1`). Memory: nvidia-smi ~57 GB/GPU at 3 blocks, so ~19 GB/block including Kokkos pools (not a high-water mark).
**Speed (Caltech):** 1.35 s/cycle at dt ~925 s on 8 H200, i.e. ~28 star-days per wall hour.
**Raven:** 24 blocks, so 24 ranks x 1 block = 6 nodes x 4 A100 (40 GB). 2 blocks/A100 (~38 GB) is too tight. The cs decomposition 6 x 1 x 2 x 2 allows 6, 12 or 24 ranks; 12 ranks would need 2 blocks/GPU. Your estimate: ~30-50 star-days/wall h, unmeasured. rst every 432000 s (5 d) plus the wall-limit rst; hst dt 250 s.
**Known issues:**
- L_top 6-13 % high from the coarse photosphere (dr 5.06 Rsun, 1.1 cells/H_p, harmonic face kappa). It converges at dr 1.2-1.6 Rsun; Kirchhoff faces recover ~3/4 (rsgphot_1106/REPORT.md). Grid change = USER decision, not for the gap.
- Mesh-wide fallback 1.6-1.9 per cycle: E<=0 in cold H-recombined downdrafts (r/R_ph 0.90-1.00, tau_top 100-1e4, not seams). `implicit_predictor_ecap=2` halves it; energy not measured yet (rsgdiag_1107/REPORT.md). Do NOT use `time2_k1_stiff` or `implicit_tstart_old` on RSG.
- hesdirk2 restarts are not bitwise in their first step (the restart lacks kappa_P/kappa_E and the tau lateral slots). Fix: `rad_m1/rst_state_ext=true` (rstbit-1107; default true in prod-1107). It helps only after a run has written an rst with it on; the 23a55658 binary does not have it.
- `implicit_fb_local` on RSG: floors x2.6 (m1_pos_flr 3044/3031 vs 1173 global path) in both tensor arms, unexplained (cstan_1106).
- Old STOP file: `prodH/STOP` exists on Caltech. It must stay there while Raven runs RSG (section 6.4).

**RSG r7 test 11-07 (Caltech, 600 cycles)**
_(placeholder: the coordinator adds the result. Arm `$J/rsgr7_1107/r7`, rst 1295, 4 H200, keys `rad_m1/cs_cross_mono=true rad_m1/implicit_opac_newton=false mesh/cs_seam_rho_guard=4 rad_m1/implicit_predictor_tcap=1.1 rad_m1/time2_k1_redist=0.1 rad_m1/implicit_fb_local=true rad_m1/implicit_fb_local_rounds=4 rad_m1/implicit_fb_diag=300`, binary md5 85587996, ended rc 0 at 09:41 PST.)_

### 2.2 HG = He giant (`hgglob_1021/prodG`)
**Physics.** Global 3-D He giant (14.5 Msun stripped-star model, MESA r15140 14p5M wide binary, X 0.0009 Y 0.978), R_ph ~103 Rsun, T_eff ~7800 K, L_top ~1.38e38 erg/s. Full sphere C128 (6 x 2x2 blocks of 704 x 64 x 64), r 3-200 Rsun, same radiation recipe as RSG plus geodesic 42 directions with vet_rotate_avg 8, cfl 0.3 (explicit), dt 80-110 s. Fresh 10-22 (N^2 surface seed). Purpose: He-giant envelope structure, convection 17.8-94 Rsun, a 53 d radial mode, settling (hg_relax_1106: r > 15 Rsun relaxed; deep part needs ~500-900 more star days).
**Now:** 580 d (rst 139, 11-06 18:59 PST, STOP for tests).
**Last restart command** (11-06 03:23 PST from rst 130): `caltech-to-raven-1107/hg_last_cmd.txt`. **Next from rst 139:** `caltech-to-raven-1107/hg_next_cmd.txt`. Keys = runG.cfg XKEYS: `time/nlim=-1 output2/dt=28800 output3/dt=28800 rad_m1/time2_vet_col=extrap2 rad_m1/time2_accept_fallback=true rad_m1/time2_posfix=4 rad_m1/implicit_work_fix=true rad_m1/cs_vertex_mono=true rad_m1/implicit_pos_floor_cons=true rad_m1/time2_auto_be=16 rad_m1/time2_k1_stiff=1.0 rad_m1/dbg_t2_admiss=20 rad_m1/vet_gd_twin=false rad_m1/vet_gd_symref=none rad_m1/vet_gd_twin_fuse=false rad_m1/vet_gd_twin_state=mean rad_m1/vet_source_noesrc=false time/dt_retry_host=false rad_m1/vet_gd_gl_nmu=0 rad_m1/vet_gd_geo=2 rad_m1/vet_rotate_avg=8 rad_m1/vet_rotate_advance=true mesh/cs_seam_blockend=extend mesh/cs_vertex_3cell=true mesh/cs_seam_rho_guard=0 time/cfl_number=0.3`. **Since your 11-04 run.cfg: `time/cfl_number=0.3` is new (USER 11-05).** HG has run 0.3 explicitly since the 11-05 20:56 restart; rst 139's parameter block already says 0.3, and the key is still named at every restart. No other change. cfg: `caltech-to-raven-1107/hg_runG.cfg`. Deck (fresh only): `hgglob_1021/hgglob_1022_n2seed.athinput`.
**Ranks/GPUs:** Caltech 12 ranks on 12 H200 (3 nodes, 2 blocks/GPU), ~49 GB/GPU, so ~25 GB/block. **Speed:** 0.97 s/cycle at dt ~80 s on 12 H200; 60-64 star-days per wall day in practice.
**Raven:** 6 nodes, 24 ranks x 1 block (25 GB fits 40 GB). rst every 432000 s; hst dt 25 s (large hst: 167 + 434 MB).
**Known issues:**
- Seam E<=0. Causes: a one-sided cross term at seams, the opacity-Newton RHS, a hot T predictor, vertex holes after cs_seam_rho_guard=0 (11-03), and the stage-1 explicit increment. The r7 candidate fixes most of these: redos 1 vs 153, wall -32 % (section 4.1; seamE_1106/REPORT.md + AUDIT.md).
- `implicit_opac_lag=pred` with 2 passes (deck): the result depends on the path, and it drives the layout excess. passes=3 cuts the lag error ~9x at +3 % wall.
- Wall layer: L_bot 4.7 % above L_in, drifting +0.23 %/100 d; odd-even v_r in the first cells. Open.
- Same hesdirk2 restart issue as RSG. `time2_adm_poslim` is harmful on HG (dropped). Corner fix vtxpos-1105 and taucut: USER "leave both off".
- Old `prodG/STOP` exists on Caltech: keep it (section 6.4).

### 2.3 R136a1 deep B (`r136a1_1104/runB/prod`)
**Physics.** Very massive star R136a1 (MESA, X 0.58, Z 0.0064), deep model "B": spherical-polar wedge 831 x 128 x 128 (4 blocks 831 x 64 x 64), tau = 1 at ~34.5 Rsun, wb_rmax 27.5 Rsun, top sponge, he_bc_outer=outflow. Implicit M1/VET, hesdirk2, cfl 0.3, dt ~5-8 s. Purpose: near-Eddington envelope, inflated layers and outflow of the most massive known star. tlim 1.728e6 s (20 d).
**Now:** 7.507 d (rst 8, wall-limit rst, 11-07 08:18 PST, rc 0). The Caltech renewal `r136B2` (job 6264757, pending) would continue from rst 8 automatically (section 6.4).
**Command:** `caltech-to-raven-1107/r136a1_cmd.txt` (last start + next from rst 8). **`rad_m1/time2_adm_poslim=true` MUST be named at every restart (USER 11-06).** It is also carried in the rst parameters. Evidence: rejections 7.8 % -> 0.6 %, 1.24 -> 1.16 s/cycle. Launcher pattern: `r136a1_1104/runB/resume_poslim.sh` (in the deps tar; STOP rule + rc 139 auto-restart). Deck (fresh only): `r136a1_1104/decks/r136a1_128_B.athinput`.
**Ranks/GPUs:** 2 ranks x 2 blocks on 2 H200, srun `-c 8 --mem=220G`. 1.15 s/cycle at dt ~5-6 s: ~0.18 star-days per wall hour (rst 5 to rst 8, 2.5 d, took 13.8 h). GPU memory was not recorded. Estimate (inferred from rst/block vs HG and RSG): 18-30 GB/block.
**Raven:** 1 node, 4 ranks x 1 block (sp 1x2x2 allows 4). Measure GPU memory in the smoke. If a block does not fit 40 GB, there is no finer split without changing the mesh.
**Known issues (r136open_1106/REPORT.md):**
- `HeStarGravity` applies the WB walked-pressure gravity source in every cell and ignores `hydro/wb_rmax`. Above 27.5 Rsun (H << dr) it over-pulls the gas by ~1000 cm/s2; the sponge hides this. Diagnostic fix (opt-in `problem/he_wb_src_linear_rmin`) in the patch file: bitwise without the key; with it -267 vs -260 cm/s2 expected at 99 Rsun, while 50-70 Rsun is still 1.4-1.9x too strong. The same code path runs in HG and RSG tops. A proper fix (honour wb_rmax) needs GPU gates and the USER.
- Real E_tot drops at 6.641 d (6.8e42 erg in E_rad = 6e-5) and 6.997 d. Sequence: a dt collapse to 0.134 s on a floor cell (r 35.4 Rsun, rho = dfloor, printed T 5e14 K), then a widening time2_posfix (33k cells), then NOT ADMISSIBLE, then a BE redo. Where the energy goes is not established (check the counters_rollback line at the end of a run). 7 dt collapses so far, all on hot low-density cells at 35-40 Rsun.
- `problem/he_sponge_mode=all`: same dt, removes the hot v_t / table-top-T cells. Candidate for a later restart (USER decision). No sponge setting raises dt; dt is set at 26-31 Rsun.
- Pockets of near-vacuum cells grow; the limiter fires on them every step.

### 2.4 WASP-121b (`w121r1_1101`)
**r1m30 (production MHD).** Ultra-hot Jupiter WASP-121b, 1x metallicity, deep_hot_jupiter_rt with correlated-k two-stream RT (ck tables in `pkg/ck`, Exo-FMS band fluxes), cs C32 (6 blocks of 112 x 32 x 32, 1 per panel), MHD with a 30 G dipole (re-initialized 11-04 from r1m8 rst 89; inert since: `problem/b_restart_done=true` in the rst), v_A cap 1e7 cm/s, cfl 0.5, dt ~2.1 s. Goal: 300 rotations. Purpose: magnetic drag on the jets and day-night flow.
**Now:** 118 rot (rst 118, 11-06 18:27 PST). Stopped for tests (USER).
**Command:** `caltech-to-raven-1107/r1m30_cmd.txt`. Launcher `alfven100_1104/launch_r1m30_resume.sh` -> `w121r1_1101/prod_rk.sub` (both in the deps tar). `mhd/wall_wb_onesided=true` must be named (USER 11-06). Env `BVALS_CC_SPLIT=1` (bitwise, 3.6 % faster). 2 ranks x 3 blocks on 2 H200, ~18 GB/GPU, srun `--mem=120G`. hst dt 1101.54 s, rst every 110154 s (1 rot).
**USER decisions pending before any resume:**
- (a) **Don't resume on the analytic tangential tide without the user.** Production has `roche_wb_tan=false` and no wb_dirs, i.e. x1-only WB + analytic -grad TidePot. USER 11-06: "we don't want the analytic tangential tide. we just want the real multi-d wb which gives the correct tide".
- (b) Resume with `hydro/wb_dirs=all` + `problem/wb_lobe_cut=true` (dhjwball-1106, merged into prod-1107; lobe cut default there), after the section 4.4 gates pass and the USER approves. USER rule: no WB of any kind outside the Roche lobe in any run; W121 x1max/r_L1 = 0.921 and 99.5 % of the top shell is cut (3e-9 of the mass, p < 2e-7 bar).
- Evidence (dhjwball_1106/REPORT.md, hydro p1 rst 302, 10 rot): the analytic-tide control spins up horizontal KE +15 % (still rising); wba stays flat; E conservation 3.9e-6 vs 5.0e-5; dt +7 %; cost +8 % per sim time. MHD r1m30 0.4 rot: stable, rms dT 42 K at 1 bar. A loosened Mach ramp is not recommended.
**Raven:** 1 node. Keep 2 ranks x 3 blocks for the first link. The old binary's ck pass control is rank-local (layout-dependent, fixed by ck_impl_global_ctl, default in prod-1107), so a different rank count changes results. Needs a dhj A100 build (CUDA, deep_hot_jupiter_rt).
**p1 (hydro):** finished at 300 rot (11-02); rst 302 is the A/B start state. A/B script pattern: `dhjwball_1106/hab2.sh` (binary dhjwball 0a029f4b3/be0475baf, keys `problem/roche_wb_tan=false hydro/wb_dirs=all problem/wb_lobe_cut=true time/tlim=... output3/bin_double=true`).

### 2.5 BSG (`doface_1105/prod_lhxo`) - tests only
**Physics.** Blue supergiant envelope (Ma+2026 / Linhao Ma's Athena++ setup) in AthenaK: sp wedge 256^3 (1x4x4 blocks of 256 x 64 x 64), his grid, floors, cfl 0.3, GM 19.78, no FOFC, `doppler_vel = vstar_nbe`, `he_bc_linhao_fis = iter` (vstarfix-1106). Fresh 11-06 14:51, finished at tlim 57.4 d 11-07 06:14 (page v24).
**Status:** USER 11-06 **dropped exact reproduction of Linhao's run** ("hopeless"). No sink/leak mimics. latke-1106 and vstarfix-1106 stay unmerged unless the USER asks. Known: base break-up at 3-5 d; KE 7.6e42 at 10 d vs his 1-5e41; light-curve amplitude ~1/4 of Ma+2026 (USER).
**Command:** `caltech-to-raven-1107/bsg_cmd.txt` (fresh start + continuation template from rst 58). 4 ranks x 4 blocks on 4 H200, `--mem=440G`, hcoll off. Deck `doface_1105/prep3d/bsgLI_lhx_exact_o.athinput`.
**Raven:** 1 node, 4 ranks x 4 blocks (rst 6.2 GB; GPU memory not recorded, small). Needs vstarfix-1106 pushed + an A100 build.
**Candidate BSG tests (none approved to run now):**

| # | test | state | USER status |
|---|---|---|---|
| B1 | Light-curve arm T: `rad_m1/vet_gd_tangential=true` from prod_lhxo rst 30 (`doface_1105/lc_1106/run_arm.sh`) | prepared 11-06 | asked for 11-06, then **disarmed**; needs a new go |
| B2 | Light-curve arm D: `time/cfl_number=0.1` from rst 30 (same script) | prepared 11-06 | **disarmed**; needs a new go |
| B3 | No-vstar arm: prod_lhxo deck without `doppler_vel = vstar_nbe` and `he_bc_linhao_fis = iter` (`prep3d/bsgLI_lhx_exact_nv.athinput`, `run_lhxnv_3d.sh`) | prepared | **USER said no** (11-06); do not run |
| B4 | GPU gates of athpp-doface-1105 (doppler_face) | gate/gpu_gates.sh | pending, low priority (repro dropped) |
| B5 | Full-tensor BSG-column check (cstan: R_on-R_off +1e-6, Q small) at production 3-D | none | candidate only |

B1/B2 need the rst 30 file (`doface_1105/prod_lhxo/rst/bsg3d.00030.rst`, 6.25 GB, md5 not computed) if the USER re-arms them.

---

## 3. Raven sizing summary

A100 40 GB, 4 per node (your NOTE-2026-11-05c/11-07). One block per GPU wherever the Caltech per-block memory is above 20 GB.

| run | blocks | Raven nodes | ranks x blocks | fits? | queue note |
|---|---|---|---|---|---|
| RSG E | 24 | 6 | 24 x 1 | ~19 GB/block (nvidia-smi incl. pools) | ~42 h estimate (11-07) |
| HG | 24 | 6 | 24 x 1 | ~25 GB/block | same |
| R136a1 B | 4 | 1 | 4 x 1 | unmeasured (est. 18-30 GB/block) | ~41 h for 1 node |
| r1m30 | 6 | 1 | 2 x 3 (as Caltech) | ~18 GB/GPU at 3 blocks | after WB gates + USER |
| BSG tests | 16 | 1 | 4 x 4 | small | per approved arm |
| HG/RSG small-state gates | 24 small | 1 GPU (gpudev 15 min possible) | 1 x 24 | hgb ran on 1 H200 | gpudev |

Fair share: per your 11-07 note, two 6-node chains for ~4 days add ~65 % to MPA's usage and lower priority for all of MPA. **The USER should choose how much of the list runs.** My suggested order, if not all of it: (1) small-state gates (1 GPU/1 node), (2) R136a1 (1 node), (3) RSG and HG (6 nodes each), (4) r1m30 after its gates, (5) BSG arms only if re-approved.

---

## 4. TEST LIST for Viper (priorities and gates)

USER rule 11-07: **"once those pass we should make them default"**: the r7 candidate keys and the full Eddington tensor (`vet_gd_tangential` on cs). **Caltech recommendation (not a USER decision):** make all r7 keys default except `time2_k1_redist`. Keep that one opt-in until the stage-1 negativity has a real fix: it moves ~15 % of L gross (2.0e37 erg/s, ~1.6 cells), donors are biased to the deeper side (88 %, rad_m1_time2.cpp:825), and it FATALs on x1-split decks.

Gate rules for every test (USER standing rules):
- GPU gates are required, not only CPU.
- Comparisons need `bin_double=true` (float32 dumps fake differences).
- Layout invariance is a required gate.
- Energy budget gate: R = dEtot - int(L_in - L_top) dt relative to int L_in dt, slope over the full window (your rule: switch if slope < -1e-4 L).
- keys-absent bitwise vs the production binary on GPU.
- Restart gate: with `rad_m1/rst_state_ext=true` (or time_scheme=be). hesdirk2 restarts without it DIFF even with no keys, so a DIFF there says nothing about the key.

### T1 (priority 1) r7 candidate on HG and RSG
Keys: `rad_m1/cs_cross_mono=true rad_m1/implicit_opac_newton=false mesh/cs_seam_rho_guard=4 rad_m1/implicit_predictor_tcap=1.1 rad_m1/time2_k1_redist=0.1 rad_m1/implicit_fb_local=true rad_m1/implicit_fb_local_rounds=4` + `rad_m1/implicit_opac_lag_passes=3` (= "r8"). All on prod-1107.
HG result so far (6 H200, rst 139, to 5.018e7): BE redos 1 vs 153, NOT POSITIVE 1 vs 184, dt 107 vs 81, wall -32 %, budget closed. Changes beyond noise: E(i280-320) +2.45 %, deep KE -0.41 % (both still growing at 1 d), L_top +0.9 %. Most of it comes from keeping D_r,lat on. Which answer is right is not shown.
Gates still open (seamE AUDIT verdict "safe to keep testing; not for production yet"):
- (a) Restart gate with rst_state_ext on the full key set incl. fb_local (GPU, hgb 1 GPU, 4 vs 2+2 cycles).
- (b) Layout gate: he_hg_k1rd b2x2 showed an excess (dens 3.9e-5 vs 3.5e-6 baseline) traced to opac lag + Newton off; check with passes=3 and with lag none; add a fb_local layout case (none exists). Cases in `tst/layout_invariance/cases.py` (layoutfix-1106). bins.cfg needs an `he_old` binary for ref_old cases, or li.py reports ERROR.
- (c) keys-absent GPU gate on the prod-1107 binary.
- (d) A longer r7 vs r7-noise (rst_pert_erad 1e-10) vs control run on HG (>= 1.5 d sim), to see whether the near-photosphere E and deep KE drifts saturate. Needs 6 nodes.
- (e) RSG: fb_local gives floors x2.6 in both arms (cstan_1106). Explain or accept before default. Caltech RSG r7 arm 11-07: see the placeholder in 2.1.
- (f) Fix the false code comments and add the x1-split FATAL in redist (AUDIT "Problems").

### T2 (priority 1) full Eddington tensor on cs (`rad_m1/vet_gd_tangential=true`)
Done on Caltech (cstan_1106/REPORT.md): stable on HG and RSG real states; no energy source (R_on - R_off <= 1e-5 of L_in); +6-8 % cost; ~2-3 % median local change in E, F near tau ~1. NOT shown more accurate: E is 2-4 % further from the SC J.
Gates still open:
- (a) **Rerun your own He 128^2 sp wedge energy test** (NOTE-2026-11-07-viper-caltech-vet-gd-tangential: XQ2nt R = +1.34, L_top +127 %) with the current compact face form on prod-1107. This is the reason it was off.
- (b) Flux-reference / beam (thin-sphere) accuracy test with an analytic reference; SC H output for a flux comparison.
- (c) Layout cases he_hg_tan, rsg_tan, he_col_gd_tan (written, not run).
- (d) GPU keys-absent bitwise on RSG, accretor, dhj MHD, R136a1 (only HG done on GPU).

### T3 (priority 2) RSG `rad_m1/implicit_predictor_ecap=2`
GPU A/B ecap2 vs base with the hst on (energy budget, Etot drift, Picard, s/cycle), 2 members, >= 1-2 h; then ecap 1.5; HG check. Caltech result: fallback 0.78 -> 0.41/cycle, bad cells 8x fewer; energy NOT measured (hst was off by mistake). Branch rsgdiag-1107 246e5e400 (in prod-1107). Full RSG state = 6 nodes; rsb small state for the first look.

### T4 (priority 2) multi-D WB (`hydro/wb_dirs=all` + `problem/wb_lobe_cut`), then the r1m30 resume
Done: keys-absent bitwise; poison test bitwise; face-potential fix c26c8e7f4; hydro 10-rot A/B; MHD 0.4 rot; static C32/C64.
Open: sp dhj gate, accretor GPU keys-absent and runtime gate, full layout gate (only 6- and 24-block cs decks tested). Then the USER go for the r1m30 resume from rst 118 with `hydro/wb_dirs=all problem/wb_lobe_cut=true` named (plus the r1m30 keys in r1m30_cmd.txt; `roche_wb_tan=false` stays). USER chose the blend (wb_weight). Speed ideas (AM probe fused with AmTransport, Taylor walk) are proposals, not done.

### T5 (priority 3) R136a1
- (a) `problem/he_sponge_mode=all` at a restart (same dt, hot v_t cells gone). USER decision; a short A/B is done (arm M).
- (b) Gravity-source fix: honour `hydro/wb_rmax` in HeStarGravity, properly, shared across he_star runs. Start from the patch. Gates: keys-absent bitwise, 50-70 Rsun residual, HG/RSG top check, GPU.
- (c) Diagnose the E_tot drops (counters_rollback; posfix after dt collapse on floor cells).

### T6 (priority 3) BSG
Only the arms the USER re-approves (section 2.5 table). No reproduction work.

### T7 (priority 1, cheap) `rad_m1/rst_state_ext` default
Default true in prod-1107. Remaining gates: the full 6-GPU HG state, RSG, an sp deck (R136a1/BSG), 2+2 vs 4 on GPU. Cost +9 % rst size (7.26e9 vs 6.67e9 bytes on hgb). The reader accepts old rst files unchanged. Also confirm the other 3 flips on GPU: cs_seam_pos_global (31/36 cs cases bitwise on CPU), ck_impl_global_ctl (dhj 7/7 CPU), wb_lobe_cut.

---

## 5. Standing USER rules that matter on Raven
- GPU gates required (CPU alone is not enough); timing comparisons on GPU, same binary, interleaved.
- Comparisons only with `bin_double=true` dumps.
- Layout invariance is a required gate before production.
- No accuracy sacrifice for speed; speed-ups must be accuracy-neutral (opt-in otherwise).
- Kill exact PIDs only, after reading the full command line; never by pattern.
- Share code across cart/sp/cs and ideal/Ohmic/ambipolar; port cs-only features to sp.
- No WB of any kind (radial included) outside the L1 equipotential in any Roche run (dhj, accretor); default after gates, name at restarts.
- R136a1: name `rad_m1/time2_adm_poslim=true` at every restart. Not for HG.
- r1m30: name `mhd/wall_wb_onesided=true`; no resume on the analytic tide without the USER.
- Link stop rule: summary NON-CONVERGED=[1-9] only.
- No production changes (keys, binaries) without the USER. Name every non-default key on the restart line.
- Very long tlim (1e9) for productions; longest walltime per link.

---

## 6. Data transfer Caltech -> Raven (and back)

### 6.1 Route
There is no automated Caltech <-> MPCDF route (NOTE-2026-11-05: past transfers went through the user). /raven/ptmp is mounted on viper, so one copy into `/raven/ptmp/jinma` serves both machines. **The copy needs the USER:** an interactive ssh to MPCDF with password + OTP. Neither session can do it. Run it from a Caltech login node (or a Caltech data-transfer node if one exists) in tmux/screen. One authenticated master connection then carries all rsyncs.

`~/.ssh/config` on Caltech (MPCDF reaches Raven from outside through the gateway; drop ProxyJump if a direct login works for you):
```
Host mpcdf-gate
  HostName gate1.mpcdf.mpg.de
  User jinma
  ControlMaster auto
  ControlPath ~/.ssh/cm-%r@%h:%p
  ControlPersist 12h
Host raven
  HostName raven.mpcdf.mpg.de
  User jinma
  ProxyJump mpcdf-gate
  ControlMaster auto
  ControlPath ~/.ssh/cm-%r@%h:%p
  ControlPersist 12h
```
Then `ssh raven true` once (password + OTP for gate and raven). After that, inside tmux:
```
J=/resnick/groups/carnegie_poc/jingze; R=/raven/ptmp/jinma
ssh raven "mkdir -p $R/caltech_1107/{deps,rst/r136a1,rst/w121,rst/bsg,rst/small} $R/raven_gap_1105/{hg,rsgE}/rst"
cd $J
# 1) deps bundle 1.17 GB (decks, tables, hst, scripts) + md5s
rsync -av --partial raven_copy_1107/raven_deps_1107.tar raven_copy_1107/raven_deps_1107.tar.md5 raven_copy_1107/hg.md5 raven_copy_1107/rsg.md5 raven:$R/caltech_1107/deps/
# 2) small: R136a1 12.1 GB, r1m30 0.34 GB, p1 0.29 GB
rsync -av --partial r136a1_1104/runB/prod/rst/r136.00008.rst raven:$R/caltech_1107/rst/r136a1/
rsync -av --partial w121r1_1101/runs/r1m30/rst/dhj.00118.rst raven:$R/caltech_1107/rst/w121/r1m30_dhj.00118.rst
rsync -av --partial w121r1_1101/runs/p1/rst/dhj.00302.rst raven:$R/caltech_1107/rst/w121/p1_dhj.00302.rst
# 3) RSG 45.4 GB and HG 97.5 GB, into your prepared run dirs
rsync -av --partial rsg_E_1020/prodH/rst/rsg.01295.rst raven:$R/raven_gap_1105/rsgE/rst/
rsync -av --partial hgglob_1021/prodG/rst/hegiant.00139.rst raven:$R/raven_gap_1105/hg/rst/
# 4) test states: hgb 6.7 GB, rsb 4.0 GB, R136a1 rst 5 12.7 GB, BSG rst 58 6.2 GB
rsync -av --partial prod_1103b/he/runs/hgb/rst/hegiant.00001.rst raven:$R/caltech_1107/rst/small/hgb_hegiant.00001.rst
rsync -av --partial prod_1103b/he/runs/rsb/rst/rsg.00001.rst raven:$R/caltech_1107/rst/small/rsb_rsg.00001.rst
rsync -av --partial r136a1_1104/runB/prod/rst/r136.00005.rst raven:$R/caltech_1107/rst/r136a1/
rsync -av --partial doface_1105/prod_lhxo/rst/bsg3d.00058.rst raven:$R/caltech_1107/rst/bsg/
```
Total: steps 1-3 ~157 GB, step 4 ~30 GB; ~187 GB in all, not multi-TB. Over ssh at an assumed 20-50 MB/s (not measured), that is ~1-2.5 h. Order by priority if time is short: deps, R136a1, RSG, HG, then the rest.

Viper/Raven side after arrival: unpack `raven_deps_1107.tar` under `$R/caltech_1107/` and run `md5sum -c` against the values in section 0 in a compute job, not on the login node (md5 of 97 GB is heavy). The hst files for `compare_hst` are in the deps tar (`hgglob_1021/prodG/*.hst`, `rsg_E_1020/prodH/*.hst`, ...); copy them to `raven_gap_1105/{hg,rsgE}/caltech_hst/` as your README expects.

### 6.2 Code
Code goes through GitHub (`origin = github.com/jing-ze-ma/athenak`), not the copy. Wait for prod-1107 and the code branches to be pushed (section 1.2).

### 6.3 Back to Caltech (when Caltech GPUs return, est. ~11-17)
On Raven: `touch FINAL` (your gaplink.sh) so the running link ends at its wall-limit rst. Then the USER pulls the newest rst from a Caltech login node: `rsync -av --partial raven:/raven/ptmp/jinma/raven_gap_1105/rsgE/rst/rsg.NNNNN.rst $J/rsg_E_1020/prodH/rst/` (similarly `hg/rst -> hgglob_1021/prodG/rst`, R136a1 -> `r136a1_1104/runB/prod/rst`, r1m30 -> `w121r1_1101/runs/r1m30/rst`). Also send back the Raven hst so the Caltech pages stay continuous. Send the md5 with the file.

### 6.4 Avoid forked histories (USER decision on Caltech)
Caltech jobs pending on priority: `hgGH3b` 6133592 (HG), `rsgEH2` 6133591 (RSG), `r1H1` 5840649 / `r1H2` 5840650, `testH1` 6133593, `r136B2` 6264757 (R136a1).
- HG and RSG coordinators exit on the existing `prodG/STOP` / `prodH/STOP`. **Keep both while Raven owns those runs.**
- **r136B2 (`r136a1_1104/runB/cont/cont.sbatch`) would continue R136a1 from the newest rst automatically.** It exits if `runB/prod/STOP` exists, and there is none now. If R136a1 moves to Raven, create that STOP file (or `scontrol hold 6264757`). Also: cont.sbatch does **not** name `rad_m1/time2_adm_poslim=true` on its command line; it relies on the value carried in the rst. That breaks the USER's naming rule; fix it before it is ever used.
- r1H1/r1H2 (`r1holders_1101/hold7d_r1.sh`) and testH1 are plain holders: they start nothing by themselves.
