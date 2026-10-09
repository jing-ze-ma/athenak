# NOTE 2026-10-09 Caltech: BSG true reproduction (answers TASK-2026-10-09-caltech-bsg-truerepro)

## PRODUCTION STARTED 2026-10-09T09:57:39 (PDT) on hpc-sm-02-11 — viper/Raven copies can be cancelled (first start wins)
- Chain 4286643 -> 4286644 -> 4286645 (12 h links, afterany), 1 node x 2 H200, 2 ranks x 8 MeshBlocks, CUDA_VISIBLE_DEVICES=SLURM_LOCALID
  wrapper, --cpu-bind=cores; fresh start (t = 0) from bsg3d_truerepro2_hr_lm.athinput (SETUP.sh, md5 2e213845 with Caltech paths);
  no command-line physics keys. Run dir /resnick/groups/carnegie_poc/jingze/bsg_1009/prod, link script bsg_1009/link_bsg.sh
  (viper rules: rst_info last_time keys on restarts, STOP on rc/FATAL/NaN/NON-CONVERGED, DONE at rst t >= 4.96e6).
- Binary: mem-1009 6c5d8fb2, PROBLEM=he_star_m1, CUDA 12.9 / gcc 13.2 / hpcx, Kokkos HOPPER90, host -ffp-contract=off:
  athena_gpu_he_star_m1_6c5d8fb2_nofma, md5 e9873137c319b5f29d85680b535a450c.
- Smoke 4286551 (nlim 10, same binary/input): rc 0, 0 FATAL/NaN/NON-CONVERGED, Picard mean 3.7 max 5, dt 87.37 s,
  GPU peak 66537 MiB (65.0 GiB, = DeltaAI 3348593), ~1.6 s/cycle after the 5.9 s setup (cpu time 21.7 s / 10 cycles)
  -> ~25 h for 57k cycles.
- He giant N897 stop record (TASK section 1, user 10-09): cancelled 09:38 PDT, jobs 4242927 (running link 2) + 4242928/4242929
  (pending); CANCEL file in n897. Resume point /resnick/groups/carnegie_poc/jingze/hegiant_1007/n897/rst/hegiant.00048.rst,
  t 2030400.0 s (23.50 d), cycle 112735 (live run was 65 cycles further); binary 7adb12d3 + hydro/sp_x2_periodic_image=true.

## Energy budget, unseeded run (prod_noseed, t 0-6.45 d)
Answers TASK-2026-10-09-caltech-bsg-energy-budget items 1-4. Script /resnick/groups/carnegie_poc/jingze/bsg_1009/budget/budget_bsg.py (reuses pagework/pipe/bsg_figs.py:
read_hst, dump_list, radial_faces, Kappa, r_photo; same tau / own-column / L(r) definitions as cmp57_dumps.py), wrapper
budget.sub (expansion, CPU). Job 4295465. Next link NOTE: `sbatch /resnick/groups/carnegie_poc/jingze/bsg_1009/budget/budget.sub /resnick/groups/carnegie_poc/jingze/bsg_1009/prod_seed`
(add `--since-days T0` for the link window; live-safe: partial hst rows and dumps < 3 min old are skipped).

Summary (start-up transient, 6.45 of 57.4 d; qualitative only):
- The budget CLOSES: L_in - <L_top> - dEtot/dt = +0.15 % L over 0-6.46 d and +0.01 % L over the last 2 d. The L_top deficit
  (0.93 L mean, 0.86 L last 0.5 d) is energy STORED in the envelope (dEtot/dt = +0.066 L whole run, +0.069 L last 2 d), not lost.
- The thin atmosphere passes the flux unchanged: in dump 12 (6.0 d) L(r) = 0.909 L at own tau_R = 1 (50.7 Rsun, nearly
  spherical in the unseeded run) and 0.912-0.913 L from 55 to 79 Rsun (no drop; dump outer-face F0 0.9124 = hst L_top 0.9124).
  So L_top < L is set below the photosphere (envelope relaxation), not by the 53-80 Rsun atmosphere / half-range top faces.
  (Page's 1.021 L at own tau = 1 was dump 11, 5.5 d; it has swung to 0.909 at 6.0 d: transient.)
- Below 40 Rsun the single-dump L(r) is dominated by transient convective fluxes (gas+radadv -0.8 L at 30 Rsun), and the cell
  lab flux m1_f1 disagrees with the transported face flux F0 by O(L) at 25-35 Rsun (shell means of large cancelling E v_r terms;
  m1_f1 is a derived cell state). diff_f / tot_f (face) are the transport fluxes; tot_c is the page Fig.-5c definition.
- hst caveats (HeStarHist @ 6c5d8fb2): L_* columns are WEDGE sums (x12 for whole star); L_top/L_bot/L_mid are the M1
  transported face flux F0 (comoving form; lab differs by O(v/c)); L_bot is the first INTERIOR face (is+1), the inner face
  itself is L_in = 1.00000 L (also from m1_face). Etot = e + KE + rho Phi_c + E_rad with Phi_c = GM(1/r_in - 1/r) (etotgrav).
  P_esrc = P_sponge = 0 here; wall mass energy terms ~1e-8 L.

Whole run (0-6.46 d; item 1 includes the last 0.5 d and 2 d means):

**Energy budget** `/resnick/groups/carnegie_poc/jingze/bsg_1009/prod_noseed`, hst t = 0.000-6.459 d (559 rows, hst dt 1024 s); L = 5.5815e+38 erg/s; wedge Omega = 1.0472 sr, all L whole-star (x 4pi/Omega = 12.00), cgs.

1. **L_top/L** (hst `L_top` = A F0 on the outer face 80.0 Rsun, M1 transported face flux, comoving form): last 0.5 d 0.8552, last 2 d 0.9308, window 0.00-6.46 d 0.9323; last row 0.8191; min/max over the last 2 d 0.816/1.022.
2. **L_in/L** (hst `L_in` = F_in A(r_in), r_in 20.0 Rsun): 1.00000 (last row; mean over window 1.00000); newest m1_face dump inner face 1.00000, first interior face (hst L_bot) 1.1841 (last row).
3. **L(r)/L**, dump 00012 (t = 6.001 d), solid-angle (dOmega) shell means x 4 pi r^2. Columns: `diff_c` = lab-frame cell F_r - (4/3) E v_r (`m1_f1`); `diff_f` = transported face flux F0 (`m1_face`, comoving, mean of the two faces of the cell); `radadv` = (4/3) E v_r; `gas` = 5/2 P_g v_r + 1/2 rho v^2 v_r; `tot_c` = page Fig.-5c total = lab F_r + gas (= diff_c + radadv + gas); `tot_f` = diff_f + radadv + gas. Mean-tau photosphere r_ph = 50.69 Rsun (page A9-A11); the adopted total is `tot_c` below r_ph and `tot_f` above it.

   | r [Rsun] | diff_c | diff_f | radadv | gas | tot_c | tot_f |
   |---|---|---|---|---|---|---|
   | 25 | -1.1749 | 0.9062 | 0.1403 | 0.1800 | -0.8546 | 1.2266 |
   | 30 | 0.4726 | 0.8996 | -0.5032 | -0.3226 | -0.3532 | 0.0738 |
   | 35 | 0.7924 | 0.9280 | -0.3928 | -0.1209 | 0.2786 | 0.4143 |
   | 40 | 0.8843 | 0.9332 | -0.0033 | -0.0008 | 0.8802 | 0.9291 |
   | 45 | 0.8934 | 0.8965 | 0.1012 | 0.0569 | 1.0515 | 1.0547 |
   | 50 | 0.9111 | 0.9111 | -0.0019 | -0.0010 | 0.9082 | 0.9082 |
   | tau_R=1 own col (r 50.7-50.7, median 50.74) | 0.9124 | 0.9124 | -0.0010 | -0.0021 | 0.9092 | 0.9092 |
   | 55 | 0.9129 | 0.9130 | -0.0008 | -0.0002 | 0.9119 | 0.9120 |
   | 60 | 0.9129 | 0.9129 | -0.0005 | -0.0000 | 0.9124 | 0.9124 |
   | 70 | 0.9127 | 0.9127 | -0.0001 | -0.0000 | 0.9126 | 0.9126 |
   | 79 | 0.9125 | 0.9124 | 0.0002 | 0.0000 | 0.9127 | 0.9126 |

   Flatness of the adopted total: 25 Rsun-r_ph min/max -0.9275/1.0956; r_ph-top min/max 0.9092/0.9126; outer face (dump) diff_f 0.9124 vs hst L_top at the nearest row 0.9124; drop from L(r_ph) to the top -0.33 %.
4. **Budget** over t = 0.000-6.459 d (6.459 d), wedge sums -> fraction of L: <L_in> 1.0000, <L_top> 0.9323, dEtot/dt 0.0661 (Etot = internal + kinetic + gravitational (rho Phi_c) + E_rad, hst `Etot`), P_esrc 0.00e+00, P_sponge 0.00e+00, top mass outflow 9.825e+19 g (wedge; energy 5.48e-09 L), removed top inflow -5.948e+16 g, inner-face mass 4.253e+17 g (energy 5.95e-12 L). **Residual L_in - <L_top> - dEtot/dt = +0.0015 L; with P_esrc - P_sponge and the wall mass terms +0.0015 L.**

Pass criteria (steady state t > 30 d): |L_top/L - 1| <= 0.03 over 2 d; |residual| <= 0.01; L(r) flat 25 Rsun-photosphere-top, no drop > 2 %. Flag rule: t > ~10 d with L_top < 0.97 L while L(tau=1) ~ 1. This run: t_end = 6.46 d (start-up transient: judge qualitatively).

Last 2 d (4.47-6.46 d), item 4:

4. **Budget** over t = 4.468-6.459 d (1.991 d), wedge sums -> fraction of L: <L_in> 1.0000, <L_top> 0.9308, dEtot/dt 0.0691 (Etot = internal + kinetic + gravitational (rho Phi_c) + E_rad, hst `Etot`), P_esrc 0.00e+00, P_sponge 0.00e+00, top mass outflow 9.825e+19 g (wedge; energy 1.78e-08 L), removed top inflow -1.272e+16 g, inner-face mass 3.710e+14 g (energy 1.68e-14 L). **Residual L_in - <L_top> - dEtot/dt = +0.0001 L; with P_esrc - P_sponge and the wall mass terms +0.0001 L.**
