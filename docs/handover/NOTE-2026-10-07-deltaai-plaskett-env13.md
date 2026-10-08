# NOTE 2026-10-07 (DeltaAI): Plaskett setup audit -> env13 (branch accretor-1007) -- q13_s1 QUEUED

Follows NOTE-2026-10-07-deltaai-plaskett (q12_s1 / q12f_s1) and the viper handover (accretor-handover/).
The accretor work now runs on DeltaAI only.

## Why: the q12 setup had bugs

- **Spherical envelope top (main bug).** `ry_per_accretor` truncated the IC envelope at a spherical r_top
  found along phi 90 (9.18 Rsun). But the envelope follows the Roche equipotentials, which bulge along the binary axis.
  The equilibrium photosphere is 9.485 Rsun at phi 0 and 9.304 at phi 180, against 9.00 at phi 90.
  - At t = 0 the gas at r_top had rho = 94 (phi 0) and 14 (phi 180), against about 3e-4 at phi 90, next to the hot ambient.
  - This gas blew out at v_r +73 km/s, reaching r 11.4 by 0.025 orbit, and settled to a surface at 9.4-9.6.
  - The mass above R_acc in the wedge was 42 rho_s Rsun^3, against a stream inflow of 35.8 per orbit. So the q12
    fluxes through R_acc (M 8.29 / 8.65 of 17.85, j 0.49 j_K) are dominated by this bulge breathing, not by accretion.
  - R_acc = 9.0 lies inside the star on the binary axis.
- **Stream width.** The default c_s/Omega gave 1.529 Rsun, against 0.735 for a stream with the Ryu+2025 width at L1,
  pressure-supported out to r_out.
- **env_rho_ph = 0.3 was 5.5x too high.** It was a carry-over from RY Per: A = (c/Omega)^2 with no 2 pi, 650 km/s.
  It is now 0.055.
- **wb_rmax** was set for an 800 km/s impact. The Plaskett normal impact speed is 339 km/s.
- **Input comments.** The tlim comment said 10 orbits but the value was 18.6, and the bin dt was P/10.75. Several
  header comments still described RY Per.
- **Spin.** The headers listed arms of 1.0 and 7.2. 7.2 is RY Per's observed spin (Barai+04), which is above
  Plaskett's critical 4.72. All runs used 1.0. **User decision 10-07: keep spin 1, no fast-spin arm.**
- **Temperatures.** These now follow the Wade+2026 MESA track at the onset of Case A.
  - Gainer: 28.3 kK, log L 4.73, giving c_ph 19.41 km/s.
  - Donor: 26.2 kK, log L 4.84, giving c_s,don 18.68 km/s (mu 0.62).
  - The old values were 30 kK and 33 kK.

## Code: accretor-1007 8cecb89a (on top of accretor-1006 b56e6e3c)

All new keys are off by default, so runs without them are bitwise unchanged.
- **`problem/env_top_mode = sphere | equipotential`.** Phi_wb is capped at Phi_s + env_top_hp c_ph^2 at every phi, and the
  ambient is anchored on the photospheric equipotential.
- **`problem/env_wb_depth`.** x1 WB reconstruction and WB gravity are switched off outside the equipotential at
  R_acc - depth. This uses the new `Hydro::wb_phimax`, checked next to wb_rmax in the PLM-WB kernels.
- **`problem/r_meas`.** Sets the face for MR/JR and Menv/Jenv. The default is R_acc exactly.
- **Build.** DeltaAI target ryp_gpu, `athena_ryp_gpu_8cecb89ac4eb`, md5 dbbef84bb366998914264b111a8a22d3.

## Inputs: plaskett-1007/env13/

- **Inputs.** `plaskett_env13.athinput` (md5 313d7f0b1156519a8eaa7063162d07bc) and `env13_nostream.athinput`.
- **Derivation.** `RESULTS.md` and `plaskett_setup_1007.py` cover sound speeds, L1 curvature, the ballistic fan,
  the transverse width, rho scaling, the photosphere table, the penetration depth and the grid.
- **Key values:**
  - env_cs_ph 19.4132, env_cs_stream 18.6791, t_don 26200
  - stream_width 0.735, env_rho_ph 0.055
  - env_top_mode equipotential, env_wb_depth 0.41 (1.5 d_pen); there is no wb_rmax
  - env_amb_rho 1e-7, dfloor = rho_amb 1e-8
  - r_meas 9.615675, which is above max r_ph + 10 H_p
  - nx1 640 plateau grid (640 x 4 x 2048 = 5.24e6 cells); spin 1.0
  - output bin P/20, rst P/10, and tlim 10 orbits in the file. The q13 runs override tlim.
- **Note.** RESULTS.md derives the grid at nx1 500 by default. The input uses `NX1=640`.

## Validation: job 3334844, 1 node, 4 GH200, PASS

- **T1 (defaults bitwise).** The old binary and the new binary give the same q12 data. Only the bin/rst headers differ,
  because they list the new default keys. hst is identical.
- **T2 (envelope only, old setup, 0.05 orbit).** The sphere top blows out: v_r 253 km/s, MR swing 0.40, edges reach
  11.4 / 13.5. The equipotential top holds: max v 6.7 km/s, MR 3e-7, edges steady.
- **T3 (env13 smoke, stream on).** dt0 1.919e-6, rc 0, 0 FATAL.
- **T4 (env13 envelope only, 0.05 orbit).** Max v 1.6 km/s, MR 7e-9.

## Production: q13_s1

- **Run.** q13_s1 is env13 with the stream, half an orbit (`time/tlim=0.2291`).
- **Layout.** 1 node ghx4-interactive, 4 GH200. The script is `q13_1n.sub`.
- **Jobs.** Job **3336769** is queued, with an estimated start of 2026-10-08 01:50 CDT. At about 51 cycles/s
  (2.7e8 zone-cycles/s) the run should take about 42 min.
- **Backup.** Job 3336438 (ghx4, 2 nodes, same output dir; it refuses to run if the dir exists) will be cancelled
  once q13_s1 has run.
- **Data.** /work/nvme/bivj/jma20/plaskett_1007/test1007/run/q13_s1/.
- Results and the comparison page against q12 will be added to this NOTE.
