# R136a1 first look (Caltech -> viper/Raven), 11-05

Files here (md5s in MD5SUMS):
- r136a1_128.athinput (md5 51a091f9): sp wedge, 512 (log r, fine band 0.0099 Rsun near 24-25 Rsun) x 128^2 over pi/3 x pi/3,
  4 MeshBlocks 512x64x64, he_star_m1, FRESH START. Phase 1: Rosseland/Planck only, NO line force (wind Mdot/v_inf not physical).
- ic_r136a1.txt: our ic_star IC from the Keszthelyi+2025 R136 MESA grid, 281.8 Msun v200 MLT++ midMS (log L 6.81, M 229 Msun,
  X_s 0.58 table / 0.555 profile), r 23.14-75 Rsun.
- rosseland_mesakap_a09_x0.58_z0.0064_topsbase.txt: MESA a09 Rosseland (TOPS below MESA's range); planck_tops_x0.58_z0.0064_a09.txt: TOPS Planck.
- caltech_cpu_smoke_512x32x32.log: CPU smoke (100 cycles, 512x32x32, 4 ranks): 0 NON-CONVERGED, Picard mean 2.9/max 5,
  no fallbacks/floors, dt 14.3 -> 11-12 s.

Paths in the deck: replace the /resnick/... prefixes (problem/he_ic_file, he_opac_table, he_planck_table and the comment)
with the Raven ic/ dir. Command line used on Caltech for the fresh start:
  athena -i r136a1_128.athinput -d <rundir> -t <wall> time/tlim=1.728e6
(no other command-line keys; outputs as in the deck). Smoke: add time/nlim=150 output2/dt=1e9 output3/dt=1e9 output5/dt=1e9 output6/dt=1e9.

Binary: the Caltech plan used prod-1102 **dc2dc20c** he_star (= the he_star tree of prod-1102 8cd7f0a9f; he_star is unchanged
between them). Please build **8cd7f0a9f** (or dc2dc20c) for this fresh start, NOT 23a55658: fresh-start defaults changed in
between (e.g. split_c2p_once and other prod1103 fresh defaults). The RSG/HG continuations stay on 23a55658.

Known (Caltech): no GPU smoke yet (no free GPU). Expected ~0.4-0.6 s/cycle on 4 H200 at 2.1 M cells/GPU (estimate);
GPU memory per block not measured (128^2 x 512 / 4 blocks = 2.1 M cells/block; RSG E 19 GB/block at 2.25 M cells/block).
Watch: R_ph(t) (the MLT++-supported compact envelope may inflate toward ~60 Rsun when the frozen convective flux ramps
off); the code capped Gamma in ~1390 fine atmosphere cells during the photosphere fit; x1min sits in MESA's MLT++
convective envelope (10 % frozen convective flux at the wall). Setup details: Caltech r136a1_1104/SETUP.md (sec. 7).
NOTE: Caltech also has a 128^2 runner armed (it starts ~11-06 04:00 on a free Caltech node). Tell us if Raven runs it and we will
hold the Caltech copy to avoid a duplicate.
