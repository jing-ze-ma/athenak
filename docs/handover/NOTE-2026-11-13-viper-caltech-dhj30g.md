# NOTE 2026-11-13 viper -> Caltech (USER asked): please run 30 G MHD for the dhj planets that disagree with the Fe I data

## Why
TOI-1518b 30 G MHD (from the hydro restart at rot 191, 12 rotations, Raven) brings the Fe I blueshift from 8.07 (hydro) to
3.00 km/s [16-84: 2.94-3.06] (Kp-marginalised, full, pc6) against the observed 3.82 +- 0.50 (Seidel 2026). 1 G stays at 7.74
and 10 G was not compatible. So a 30 G dipole fixes most of the hydro over-blueshift. Hydro round-1 status (Fe I pc6 blueshift,
km/s, model vs Seidel 2026; viper/Raven runs, FE_HYDRO_1112):

| planet | hydro model (rot) | observed | | 30 G |
|---|---|---|---|---|
| TOI-1518b | 8.07 (191) | 3.82 +- 0.50 | ~8 sigma | done on Raven (3.00) |
| WASP-178b | 8.80 (92) | 4.79 +- 1.19 | ~3.4 sigma | **Raven** (being set up now) |
| WASP-189b | 9.02 (108) | 1.05 +- 1.59 | ~5 sigma | **Caltech, please** |
| WASP-76b | 7.62 (305) | 6.95 +- 0.36 | ~1.8 sigma | no (roughly consistent) |
| TOI-2109b | 5.00 (21) | none | | no |

## Request
1. Run 30 G MHD for **WASP-189b** from your latest hydro round-1 restart of WASP-189b (state its rotation), and for any other
   planet in your dhj set whose hydro Fe I is > ~2 sigma from the observation and is not in the Raven column above
   (tell us which before starting, so we do not duplicate). WASP-121b is your production: your call.
2. Use the **latest binary: fork/prod-1118 = 1d1972805** (contains the dhj line, the dhjfast speed-ups c3ee3a77c and the
   MHD keys below), built as your usual Caltech recipe with PROBLEM = deep_hot_jupiter_rt (CUDA Release, MPI on,
   `-ffp-contract=off --fmad=false`, HOPPER90). Smoke 10 cycles first (rc 0, FATAL 0).
3. Recipe = our TOI-1518b 30 G run (Raven /raven/ptmp/jinma/toi1518b_mhd30_hyd_1112, deck run.athinput, keys run.cfg):
   - fresh MHD start from the hydro restart (remap of hydro rst), polar dipole set by the pgen at t = 0: B_r = bbot cos(theta)
     (r0/r)^3 with `problem/bbot_gauss` chosen so that the field is **30 G at R_p** (TOI-1518b: bbot_gauss = 36.566373 at r0);
     ohmic resistivity from the eos table, RKG2 STS, same W121-MHD deck lineage + the planet's keys.
   - MHD keys: `mhd/wb_dirs=all problem/wb_lobe_cut=true problem/ck_impl_global_ctl=true mesh/cs_seam_pos_global=true
     mesh/cs_seam_rho_guard=4 problem/roche_wb_tan=false mhd/wall_wb_onesided=true mesh/cs_seam_fc_order=3
     mesh/cs_vertex_3cell=true mhd/wb_plm_range=true mhd/top_ghost_prim_copy=true mhd/wb_plm_smooth=true
     mhd/lhlld_alfmach=true mhd/valfven_max_cgs=5e6 mhd/top_bc_bt=sourcesurface problem/ck_impl_every=14
     problem/ck_impl_xstep=28 problem/ck_vef_every=520 time/cfl_number=0.5`.
   - run to t0 + 12 rotations, then Fe I with fe3d_gpu (fehires_fast, pcuts 1e-6, pc6) + ana.py on the last 4 dumps;
     report the "Kp-marg full (a) median [16, 84]" and limb mor/eve, as in our ana_summary.
4. Holders: r1H1 (1 x 4 H200) is the natural place when it starts. One planet fits one node.

Reply with a NOTE (planet, restart/rotation, binary md5, smoke result, ETA), and the Fe I result when done. Times in PST.
