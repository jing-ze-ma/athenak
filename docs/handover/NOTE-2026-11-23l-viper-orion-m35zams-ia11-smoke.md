# NOTE 2026-11-23l viper -> orion: USER — (1) smoke the new IA11 M35ZAMS deck on Raven; (2) remove the GL phi-period constraint

Times PST. Thanks for 23k (successors queued: viper pr136f 13610529, pool 13610530; Raven pT9L 32910619, pHG262 32910620).

## 1. New fresh run: M35ZAMS (Will's T42L5.0, Schultz+2022) on IA11, for Caltech
USER (11-23 ~09:30-10:40): prepare Will's other M35 star "in the same style as the runs currently running on Caltech", on IA11,
with the lateral wedge DOUBLED at the paper's cell size (S22 has 128 x 128 over 0.064 x 0.064 rad; now 256 x 256 over
0.128 x 0.128 rad), and rays "close to 120": vet_gd_ls = 10 (S10, 120 directions).
Files in docs/handover/m35zams_1123/ (MD5SUMS):
- `IA11_m35zams.athinput` (md5 53f15c72): Caltech paths, 4 blocks 384x128x128 (1 per H200). Built by `mkdeck_zams.py` from
  Caltech's AS-RUN IA9_m35mms.athinput (md5 52da1397) + your 4 IA11 keys (fwb_stage/fwb_t2 false, w_fse_tau, gtime_mom;
  `template_IA11_m35mms.athinput`, checked key by key against your ia11_m35mms.athinput). Only star + grid keys differ:
  mesh 384 log r 4.7294e11-5.7031e11 (6.798-8.200 Rsun, S22) x theta [1.506708, 1.634885] x phi [0, 0.128]; mu 0.6220 +
  temp_unit_kelvin / arad / tfloor (1000 K) from it; implicit_flux_x1min = L/(4 pi r_in^2) = 2.35463e14; he_ic_file =
  fresh_prod_1120/m35zams/ic_m35zams_ideal.txt (ideal IC of 11-20, Caltech has it, gzipped); rays vet_gd_gl_nmu 0 + vet_gd_ls 10.
- Two values with no recorded FRESHPREP rule; please confirm or replace with your builder's rule:
  he_ic_balance_rmax 5.353221e11 (IC radius at T = 2.137e5 K, the temperature at M35MMS's 8.776579e11) and
  wb_rmax 6.691526e11 = 1.25 x that (the M35MMS and T19L6.4 FRESHPREP ratio; > x1max as in M35MMS).
- `raven_smoke_IA11_m35zams.athinput` (md5 ab323bb3): the same with Raven paths (IC already copied to
  /raven/ptmp/jinma/stag2/m35zams_1123/, md5 b480f870; your ia9b/cal tables), meshblock nx3 64 (8 blocks for 40 GB A100s;
  64 < the 127-cell near-tangent reach, so expect clamped-read warnings: vet_gd_allow_clamp is on), hst only.
  Both decks are also in /raven/ptmp/jinma/stag2/m35zams_1123/.
- USER: please smoke it on Raven with your IA11 A100 build (calqc_cu bd9553537, 35b66b73), like your ia11go smokes
  (200 cycles, fresh), when your tests leave nodes free. Please post: FATAL/CW/fallbacks/dt_retry, dt, s/cycle,
  L_top/L_in, L_bot/L_in, peak GPU memory, and the clamp count (and, if a 4-node 384x128x128-block layout fits in 40 GB,
  a clamp-free number too). Viper's own fresh IA9 smoke (55 cycles, 2 MI300A nodes, GL 7x18) stopped at cycle 0 on the
  constraint below; the n2 seed band came out fine (r 5.370-5.444e11, l 16-64).
- Expected cost at Caltech: dt ~3.5 s (the 11-15 ZAMS runs), 10 d ~ 245k cycles, ~0.9 s/cycle by M35MMS block scaling
  => ~60 h on 1 node (4 H200). Not launched; USER decides after the smoke.

## 2. USER: "figure out how to make that constraint not there anymore"
- The GL product set FATALs unless 2 pi/vet_gd_gl_nphi divides the wedge's phi period (rad_m1_vetgd.cpp:978-986 on
  bd9553537: m = per nphi/2 pi must be an integer). With GL 7x18 the phi extent must be a multiple of 20 deg, so a
  wedge narrower than 20 deg at 256 cells (here 0.128 rad) cannot use the Caltech ray set; we fell back to S10.
  Please make the GL set usable on any wedge (e.g. the inexact interpolated wall map that the other sets already use,
  vet_gd_wall_interp, or an exact map for the nearest rotation), with a gate showing the walls are as accurate as on a
  20 deg wedge.
- Related (for information, lower priority): the vet_gd near-tangent lateral reach (here 127 cells, set by dr/r vs dtheta)
  forces MeshBlocks >= the reach for an exact sweep, so the paper's 128 x 128 wedge can only run exactly as ONE block on
  ONE GPU. If the raw band can cover a reach of several blocks exactly, say so.

-- viper
