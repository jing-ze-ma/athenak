# TASK 11-07 (USER go): BSG 3-D run WITHOUT vstar on Raven — start this FIRST

User (Jingze), 11-07 ~10:30 PST: "tell viper to start a bsg run without vstar first". Priority: start before the other items in HANDOVER-2026-11-07-caltech-to-raven.md.

## What
Fresh start t = 0 -> 57.4 d of the BSG Linhao-exact 3-D model, identical to the finished Caltech production `prod_lhxo` (ended rc 0 at 57.4 d 11-07 06:14 PST) except that the radiation Doppler/work velocity is the fluid velocity instead of the stellar-frame vstar:
- deck: `caltech-to-raven-1107/bsgLI_lhx_exact_novstar.athinput` (this branch) = `doface_1105/prep3d/bsgLI_lhx_exact_o.athinput` with TWO changes (USER 11-07: "also try to use full tensor for it"): (1) line removed: `<rad_m1> doppler_vel = vstar_nbe`; (2) `<rad_m1> vet_gd_tangential = true` (was false): the full Eddington tensor on sp (tangential anisotropy D_tt-D_pp + D_tp from the gd rays above the cut). Everything else identical, including `<problem> he_bc_linhao_fis = iter` (kept on purpose: the only differences vs prod_lhxo are vstar and the full tensor). Note: the older Caltech-prepared `bsgLI_lhx_exact_nv.athinput` also dropped `he_bc_linhao_fis`; do not use it.
- code: branch `vstarfix-1106` 3b1d3dc19 (now on origin), the same commit as the production binary (Caltech md5 20efaed08d5ab3789073112597af3ed9, built `athena_gpu_he_star_m1_v3_fmadoff`: he_star_m1 pgen, FMA off). Build it on Raven with the same cmake options as your he_star_m1 builds, -DKokkos_ARCH for A100, FMA off. (prod-1107 also contains newer code; do NOT use it for this run: the comparison must be same code, vstar on vs off.)
- tables / IC (absolute Caltech paths in the deck: repoint to your copies): `bsgli_1014/files/linhao_profile_full.txt`, `ic_ma2026_f0.txt`, `rosseland_ma2026_x0.7_z0.008_padded.txt`, `planck_ma2026_x0.7_z0.008_padded.txt`; all four are in `raven_deps_1107.tar` (arriving now via the user's transfer to `/raven/ptmp/jinma/caltech_1107/deps/`; first file copied). Deck lines 253, 272-274.
- command: `athena -i bsgLI_lhx_exact_novstar.athinput time/tlim=4.95936e6` (plus repointed table paths), env as in `caltech-to-raven-1107/bsg_cmd.txt`.

## Size
256 x 256 x 256 in 16 MeshBlocks 256x64x64 (deck comment: 1 per A100 on 4 nodes x 4, or 4 per GPU). Caltech: 4 ranks x 4 blocks on 4 H200, ~25 GB/GPU, 61 H200 GPU-h for 57.4 d (~15 h wall). On Raven: 1 node, 4 ranks x 4 blocks (25 GB fits 40 GB A100) ≈ 1.9x slower ≈ ~30 h wall; or 4 nodes x 4 ranks x 1 block if queue allows (faster, more nodes).

## STEP 0 (USER 11-07: "make sure full tensor works on it first"): validate the full tensor on this BSG deck BEFORE the long run
Short arms from t = 0 with this deck (same binary, same node layout), ~1-2 h each, all compared at equal sim time with bin_double dumps:
- A = deck as is (vstar off, tensor ON); B = same with vet_gd_tangential = false; Bn = B with E x (1+1e-10 u) (noise member; Caltech key rst_pert_erad or your equivalent).
- Required to pass, A vs B (vs the B-Bn noise):
  1. runs clean: no FATAL, no NON-CONVERGED=[1-9], no nan; NOT POSITIVE / NOT ADMISSIBLE / BE-redo / fallback rates not higher than B beyond noise; dt not lower;
  2. the term is actually ON: the tangential term is applied in most solves (not dropped by the mesh-wide fallback); report the active fraction (implicit_fb_diag / FBDIAG_SUM or your equivalent) and the size of D_tt-D_pp, D_tp where tau_top < 30;
  3. energy gate: R = dEtot - int(L_in - L_top) dt; slope(A) - slope(B) within 1e-4 L_in and not growing (your He-wedge source was +1.34 / +3.2e-2 with wall_interp, which is ON in this deck);
  4. no imprint at the poles / wedge walls / block edges in E, F_r, F_t, T_eff maps (sp-specific: check the axis and the wall_interp walls);
  5. restart 2 + 2 vs 4 cycles bitwise (tensor on), and a quick 2-layout check if cheap.
- If all pass: start the long run (tensor on). If any fails: report the numbers, start the long run with vet_gd_tangential = false (vstar off only) so the user's main question still gets answered, and say so in the NOTE.

## Full-tensor energy gate during the long run (REQUIRED, your NOTE-2026-11-07 vet-gd-tangential)
Your He wedge showed an energy source with vet_gd_tangential on sp (T7wi sink, XQ2nt R=+1.34; with wall_interp +3.2e-2); this deck has vet_gd_wall_interp = true. Caltech saw no source on cs (HG/RSG) or in short sp runs (R136a1, BSG column), but never on your He state. So: run the energy gate R = dEtot - int(L_in - L_top) dt relative to int L_in dt from the start; STOP and report if |slope| > 1e-4 L_in over the first ~2 d (or if L_top runs away). Optional if nodes allow: a second arm with vet_gd_tangential = false (vstar off only) to separate the two changes; the user asked for one run with both, so the second arm is lower priority.

## Stop rules / checks
- Stop on FATAL, `NON-CONVERGED=[1-9]` in the summary line, or nan (link stop rule; per-pass "NON-CONVERGED after" with DT_RETRY is not a stop).
- Compare against Caltech prod_lhxo hst (`doface_1105/prod_lhxo/bsg3d.{hydro,user}.hst` in the deps tar): KE, L_top light curve (amplitude: prod_lhxo ~0.98-1.02, Ma+2026 0.91-1.09), R_ph, dt. The question is whether vstar suppresses the light-curve amplitude (ours ~1/4 of Ma+2026). Report at 10, 30 and 57.4 d.
- Output cadence as in the deck (bin_double for any comparison dumps).
- Kill exact PIDs only; GPU gates rule n/a (this is a production-style arm with gated code).

Co-ordinated from Caltech (session of 11-07); reply with a NOTE on this branch.
