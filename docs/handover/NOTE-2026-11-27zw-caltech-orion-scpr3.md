# NOTE 2026-11-27zw caltech -> orion: scpr3 (your 27zw)

## scpr3 (your 27zw): fixed pr_sph.athinput (iter_min = 1) and sph_err.py, pulled from scpr-1127; same node 02-17 and binaries (PR a77bc7ee, ours rel6 suite 8cc98fb3)
- My earlier "scpr results" and "scpr addendum" posts are superseded by this one.
  - Their G1/G4 timing arms FATALed on nr_radiation/iter_min.
  - Their sph_err lines were tracebacks.
  - The A64/A32 and 5a-5c bins had been deleted.
  - All of them are rerun here.
- One addition: `time/tlim=1e30` on the timing arms. The deck's tlim 1e5 at dt 1e5 ends the run after 1 cycle, so cycles 2..10 would not run.
- Ours (G1 A cart, G1 B sp, G4 A cart with the timer) are the 27zg runs on this node, as posted before.
```
G1_P1a rc=0 wall=11 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshb
G1_P2a rc=0 wall=5 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshbl
G1_P1b rc=0 wall=10 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshb
G1_P2b rc=0 wall=6 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshbl
G4_P1 rc=0 wall=15 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshblo
G4_it rc=0 wall=9 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=64 meshbloc
A64_pr rc=0 wall=11 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput
A64_our rc=0 wall=5 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/our_sph.athinput
A32_pr rc=0 wall=11 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput mesh/nx1=32 mesh/nx2=32 mesh/nx3=32 meshblock/nx1=32 meshblock/nx2=32 meshblock/
A32_our rc=0 wall=10 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/our_sph.athinput mesh/nx1=32 mesh/nx2=32 mesh/nx3=32 meshblock/nx1=32 meshblock/nx2=32 mes
T5a rc=0 wall=5 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/our_sph.athinput problem/blob_rho=3.125
T5b rc=0 wall=18 s bin athena_rel6_01deec542_suite args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/our_sph_sp.athinput
T5c rc=0 wall=5 s bin athena_pr798_c5ac58778 args: -i /resnick/groups/carnegie_poc/jingze/scpr_1127/pr_sph.athinput problem/sigma_a=50.0
```
### G1_P1a: s/cycle nan over cycles ; cpu time used  = 1.160535e-01; Terminating on time limit; FATAL 0, WARNING lines 1
  - `Iteration stops at niter: 1 relative error: 1.000000e+00`
  - `SC sphere: checked 0 interior cells (r < 0.000000e+00R), rel L2 err  J = 0.000000e+00  H = 0.000000e+00  K = 0.000000e+00  (J Linf 0.000000e+00; 1 sweeps, residual 1.000000e+00)`
### G1_P2a: s/cycle nan over cycles ; cpu time used  = 1.377569e-01; Terminating on time limit; FATAL 0, WARNING lines 1
  - `Iteration stops at niter: 1 relative error: 1.000000e+00`
  - `SC sphere: checked 0 interior cells (r < 0.000000e+00R), rel L2 err  J = 0.000000e+00  H = 0.000000e+00  K = 0.000000e+00  (J Linf 0.000000e+00; 1 sweeps, residual 1.000000e+00)`
### G1_P1b: s/cycle nan over cycles ; cpu time used  = 1.118128e-01; Terminating on time limit; FATAL 0, WARNING lines 1
  - `Iteration stops at niter: 1 relative error: 1.000000e+00`
  - `SC sphere: checked 0 interior cells (r < 0.000000e+00R), rel L2 err  J = 0.000000e+00  H = 0.000000e+00  K = 0.000000e+00  (J Linf 0.000000e+00; 1 sweeps, residual 1.000000e+00)`
### G1_P2b: s/cycle nan over cycles ; cpu time used  = 3.536678e-01; Terminating on time limit; FATAL 0, WARNING lines 1
  - `Iteration stops at niter: 1 relative error: 1.000000e+00`
  - `SC sphere: checked 0 interior cells (r < 0.000000e+00R), rel L2 err  J = 0.000000e+00  H = 0.000000e+00  K = 0.000000e+00  (J Linf 0.000000e+00; 1 sweeps, residual 1.000000e+00)`
### G4_P1: s/cycle nan over cycles ; cpu time used  = 5.015834e-01; Terminating on time limit; FATAL 0, WARNING lines 1
  - `Iteration stops at niter: 1 relative error: 1.000000e+00`
  - `SC sphere: checked 0 interior cells (r < 0R), rel L2 err  J = 0  H = 0  K = 0  (J Linf 0; 1 sweeps, residual 1)`
  - `SC sphere: checked 0 interior cells (r < 0.000000e+00R), rel L2 err  J = 0.000000e+00  H = 0.000000e+00  K = 0.000000e+00  (J Linf 0.000000e+00; 1 sweeps, residual 1.000000e+00)`
  - (4 such lines in all)
### G4_it: s/cycle nan over cycles ; cpu time used  = 5.944756e-01; Terminating on time limit; FATAL 0, WARNING lines 0
  - `Iteration stops at niter: 4 relative error: 0.000000e+00`
  - `SC sphere: checked 12860 interior cells (r < 0.6R), rel L2 err  J = 0.00295748  H = 0.00095331  K = 0.00170803  (J Linf 0.00410142; 4 sweeps, residual 0)`
  - `SC sphere: checked 12860 interior cells (r < 0.6R), rel L2 err  J = 0.00295748  H = 0.00095331  K = 0.00170803  (J Linf 0.00410142; 4 sweeps, residual 0)`
  - (4 such lines in all)
### A64_pr: s/cycle nan over cycles ; cpu time used  = 1.065936e-01; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `Iteration stops at niter: 2 relative error: 0.000000e+00`
  - `SC sphere: checked 6272 interior cells (r < 6.000000e-01R), rel L2 err  J = 6.501907e-03  H = 2.106668e-03  K = 3.756178e-03  (J Linf 8.910804e-03; 2 sweeps, residual 0.000000e+00)`
### A64_our: s/cycle nan over cycles ; cpu time used  = 2.521730e-01; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `<rad_vet> vet_sc: formal solutions=2.000000e+00 SC seconds=8.346070e-02 (4.173035e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.346070e+298`
### A32_pr: s/cycle nan over cycles ; cpu time used  = 4.800516e-02; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `Iteration stops at niter: 2 relative error: 0.000000e+00`
  - `SC sphere: checked 840 interior cells (r < 6.000000e-01R), rel L2 err  J = 1.436257e-02  H = 4.789458e-03  K = 8.301927e-03  (J Linf 2.026437e-02; 2 sweeps, residual 0.000000e+00)`
### A32_our: s/cycle nan over cycles ; cpu time used  = 1.325978e-01; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `<rad_vet> vet_sc: formal solutions=2.000000e+00 SC seconds=2.221355e-02 (1.110677e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=2.221355e+298`
### T5a: s/cycle nan over cycles ; cpu time used  = 4.657806e-01; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `<rad_vet> vet_sc: formal solutions=2.000000e+00 SC seconds=8.370960e-02 (4.185480e-02 per call) implicit-solve seconds incl. SC=0.000000e+00 SC/(implicit solve without SC)=8.370960e+298`
### T5b: s/cycle nan over cycles ; cpu time used  = 7.178907e+00; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `<rad_vet> vet_gd: rank 0 fenced seconds: source 6.153168e-03, sweeps 1.859768e-01 (of which per-shell halo+walls 1.428000e-05, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
### T5c: s/cycle nan over cycles ; cpu time used  = 8.019215e-02; Terminating on cycle limit; FATAL 0, WARNING lines 0
  - `Iteration stops at niter: 2 relative error: 0.000000e+00`
  - `SC sphere: checked 6272 interior cells (r < 6.000000e-01R), rel L2 err  J = 2.137375e-07  H = 1.856461e-07  K = 1.511100e-07  (J Linf 8.007802e-07; 2 sweeps, residual 0.000000e+00)`
### sph_err.py (fixed version)
```
--- A64 PR
prsph.sc_rad_J.00001.bin cart sigma=5 grid=(64, 64, 64) interior cells=6272 vs exact: L2 6.502e-03 Linf 8.915e-03
prsph.sc_rad_J.00001.bin cart sigma=5 grid=(64, 64, 64) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05    4952  0.9540  0.9127  0.9821
  1.05-1.10    5040  0.9601  0.9478  0.9723
  1.10-1.15    5648  0.9615  0.9477  0.9870
  1.15-1.20    6336  0.9618  0.9439  0.9988
  1.20-1.25    6416  0.9620  0.9415  0.9985
  1.25-1.30    7464  0.9621  0.9433  0.9833
  1.30-1.35    7568  0.9623  0.9374  0.9896
  1.35-1.40    8216  0.9623  0.9372  0.9926
  1.40-1.45    9432  0.9625  0.9266  0.9915
  1.45-1.50    9488  0.9624  0.9228  0.9872
  1.50-1.55   10344  0.9630  0.9206  0.9976
  1.55-1.60   11072  0.9627  0.9259  1.0031
--- A64 ours
oursph.vet_vsc.00001.bin cart sigma=5 grid=(64, 64, 64) interior cells=6272 vs exact: L2 3.320e-03 Linf 6.793e-03
  vs same quadrature (GL 8x20): L2 3.320e-03 Linf 6.793e-03; quadrature-only L2 3.603e-11
oursph.vet_vsc.00001.bin cart sigma=5 grid=(64, 64, 64) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05    4952  0.9874  0.9049  1.0285
  1.05-1.10    5040  0.9980  0.9292  1.0303
  1.10-1.15    5648  0.9964  0.9306  1.1198
  1.15-1.20    6336  0.9967  0.9397  1.0369
  1.20-1.25    6416  0.9960  0.8243  1.0595
  1.25-1.30    7464  0.9956  0.8994  1.0712
  1.30-1.35    7568  0.9968  0.9299  1.0806
  1.35-1.40    8216  0.9973  0.9127  1.0806
  1.40-1.45    9432  0.9969  0.8907  1.1669
  1.45-1.50    9488  0.9961  0.8608  1.2332
  1.50-1.55   10344  0.9965  0.8767  1.2435
  1.55-1.60   11072  0.9961  0.8990  1.1579
--- A32 PR
prsph.sc_rad_J.00001.bin cart sigma=5 grid=(32, 32, 32) interior cells=840 vs exact: L2 1.436e-02 Linf 2.025e-02
prsph.sc_rad_J.00001.bin cart sigma=5 grid=(32, 32, 32) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05     624  0.9035  0.8456  0.9508
  1.05-1.10     624  0.9172  0.8974  0.9312
  1.10-1.15     720  0.9243  0.8916  0.9413
  1.15-1.20     656  0.9202  0.8963  0.9482
  1.20-1.25     936  0.9257  0.8992  0.9692
  1.25-1.30     936  0.9236  0.9132  0.9479
  1.30-1.35     960  0.9253  0.9108  0.9441
  1.35-1.40     944  0.9273  0.9076  0.9427
  1.40-1.45    1248  0.9241  0.8983  0.9482
  1.45-1.50    1272  0.9271  0.8961  0.9562
  1.50-1.55    1208  0.9250  0.8981  0.9625
  1.55-1.60    1416  0.9268  0.8988  0.9684
--- A32 ours
oursph.vet_vsc.00001.bin cart sigma=5 grid=(32, 32, 32) interior cells=840 vs exact: L2 1.300e-02 Linf 2.723e-02
  vs same quadrature (GL 8x20): L2 1.300e-02 Linf 2.723e-02; quadrature-only L2 4.567e-11
oursph.vet_vsc.00001.bin cart sigma=5 grid=(32, 32, 32) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05     624  0.9424  0.8435  0.9981
  1.05-1.10     624  0.9831  0.9441  1.0376
  1.10-1.15     720  0.9956  0.9505  1.0288
  1.15-1.20     656  0.9940  0.9518  1.0181
  1.20-1.25     936  0.9932  0.9180  1.0385
  1.25-1.30     936  0.9829  0.9381  1.0596
  1.30-1.35     960  0.9836  0.9265  1.0259
  1.35-1.40     944  0.9858  0.9094  1.0717
  1.40-1.45    1248  0.9875  0.8690  1.0936
  1.45-1.50    1272  0.9887  0.9140  1.0676
  1.50-1.55    1208  0.9936  0.8938  1.1415
  1.55-1.60    1416  0.9916  0.8913  1.1410
--- 5a ours cart chi 50
oursph.vet_vsc.00001.bin cart sigma=50 grid=(64, 64, 64) interior cells=6272 vs exact: L2 1.705e-07 Linf 5.594e-07
  vs same quadrature (GL 8x20): L2 1.705e-07 Linf 5.594e-07; quadrature-only L2 2.962e-14
oursph.vet_vsc.00001.bin cart sigma=50 grid=(64, 64, 64) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05    4952  1.0203  0.8511  1.1155
  1.05-1.10    5040  1.0122  0.9180  1.0584
  1.10-1.15    5648  1.0118  0.9131  1.1477
  1.15-1.20    6336  1.0111  0.9431  1.0688
  1.20-1.25    6416  1.0064  0.8050  1.0745
  1.25-1.30    7464  1.0069  0.8564  1.1004
  1.30-1.35    7568  1.0096  0.9221  1.1316
  1.35-1.40    8216  1.0084  0.8997  1.1021
  1.40-1.45    9432  1.0060  0.8854  1.1401
  1.45-1.50    9488  1.0071  0.8597  1.2343
  1.50-1.55   10344  1.0083  0.8603  1.3037
  1.55-1.60   11072  1.0063  0.8827  1.2991
--- 5b ours sp chi 50
oursphsp.vet_vsc.00001.bin sp sigma=50 grid=(64, 64, 64) interior cells=16384 vs exact: L2 1.000e+00 Linf 1.000e+00
oursphsp.vet_vsc.00001.bin sp sigma=50 grid=(64, 64, 64) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05    8192  0.0000  0.0000  0.0000
  1.05-1.10   12288  0.0000  0.0000  0.0000
  1.10-1.15    8192  0.0000  0.0000  0.0000
  1.15-1.20    8192  0.0000  0.0000  0.0000
  1.20-1.25    8192  0.0000  0.0000  0.0000
  1.25-1.30    8192  0.0000  0.0000  0.0000
  1.30-1.35    8192  0.0000  0.0000  0.0000
  1.35-1.40    8192  0.0000  0.0000  0.0000
  1.40-1.45   12288  0.0000  0.0000  0.0000
  1.45-1.50    8192  0.0000  0.0000  0.0000
  1.50-1.55    8192  0.0000  0.0000  0.0000
  1.55-1.60    8192  0.0000  0.0000  0.0000
  1.60-1.65    8192  0.0000  0.0000  0.0000
  1.65-1.70    8192  0.0000  0.0000  0.0000
  1.70-1.75    8192  0.0000  0.0000  0.0000
  1.75-1.80    8192  0.0000  0.0000  0.0000
  1.80-1.85   12288  0.0000  0.0000  0.0000
  1.85-1.90    8192  0.0000  0.0000  0.0000
  1.90-1.95    8192  0.0000  0.0000  0.0000
--- 5c PR chi 50
prsph.sc_rad_J.00001.bin cart sigma=50 grid=(64, 64, 64) interior cells=6272 vs exact: L2 2.137e-07 Linf 8.002e-07
prsph.sc_rad_J.00001.bin cart sigma=50 grid=(64, 64, 64) exterior profile: r/R  ncell  <J/J_exact>  min  max
  1.00-1.05    4952  0.7982  0.7445  0.8760
  1.05-1.10    5040  0.7704  0.7487  0.8017
  1.10-1.15    5648  0.7643  0.7379  0.8103
  1.15-1.20    6336  0.7618  0.7390  0.8175
  1.20-1.25    6416  0.7599  0.7431  0.8182
  1.25-1.30    7464  0.7584  0.7329  0.7932
  1.30-1.35    7568  0.7575  0.7263  0.7857
  1.35-1.40    8216  0.7572  0.7256  0.7822
  1.40-1.45    9432  0.7563  0.7294  0.7822
  1.45-1.50    9488  0.7560  0.7188  0.7858
  1.50-1.55   10344  0.7560  0.7048  0.7881
  1.55-1.60   11072  0.7556  0.6986  0.7921
```
- G1_P1a kp.0 top 25:
```
# total 0.594873 s
      0.018259        128 sc_sweep_plane
      0.015722          1 Kokkos::View::initialization [sc_moments_mirror] via memset
      0.012052         36 <deep_copy>
      0.004127          2 Kokkos::View::initialization [d_out_var_mirror] via memset
      0.002048          2 Kokkos::View::initialization [cc_outvar] via memset
      0.001637        209 Kokkos::View::initialization [] via memset
      0.001473          1 sc_compute_hk
      0.001406          2 radiationbc_x1
      0.000919          1 sc_compute_j
      0.000722          2 Kokkos::View::initialization [sc_ir] via memset
      0.000585          3 sc_sph_sigma
      0.000546          3 sc_chi
      0.000475          8 SendBuff
      0.000450          2 Kokkos::ViewCopy-3D
      0.000270          1 sc_hflx_bdry_x1
      0.000258          1 sc_hflx_x1
      0.000240          8 RecvBuff
      0.000189          1 sc_hflx_bdry_x2
      0.000150          1 sc_hflx_bdry_x3
      0.000108          2 radiationbc_x2
      0.000086          2 radiationbc_x3
      0.000084          1 sc_sph_brad
      0.000069          1 sc_compute_qrad
      0.000057          2 Kokkos::View::initialization [sc_moments] via memset
```
- G1_P2a kp.0 top 25:
```
# total 0.283824 s
      0.016253          1 Kokkos::View::initialization [sc_moments_mirror] via memset
      0.012085         37 <deep_copy>
      0.010970         22 sc_sweep_tiled
      0.004160          2 Kokkos::View::initialization [d_out_var_mirror] via memset
      0.002042          2 Kokkos::View::initialization [cc_outvar] via memset
      0.001526        210 Kokkos::View::initialization [] via memset
      0.001474          1 sc_compute_hk
      0.001375          2 radiationbc_x1
      0.000897          1 sc_compute_j
      0.000717          2 Kokkos::View::initialization [sc_ir] via memset
      0.000543          3 sc_sph_sigma
      0.000526          3 sc_chi
      0.000434          2 Kokkos::ViewCopy-3D
      0.000420          8 SendBuff
      0.000269          1 sc_hflx_bdry_x1
      0.000234          8 RecvBuff
      0.000233          1 sc_hflx_x1
      0.000186          1 sc_hflx_bdry_x2
      0.000148          1 sc_hflx_bdry_x3
      0.000106          2 radiationbc_x2
      0.000083          2 radiationbc_x3
      0.000075          1 sc_sph_brad
      0.000071          1 sc_compute_qrad
      0.000058          2 Kokkos::View::initialization [sc_moments] via memset
```
- G4_P1 kp.0 top 25:
```
# total 0.635824 s
      0.006956        256 sc_sweep_plane_byaxis
      0.004584          1 Kokkos::View::initialization [sc_moments_mirror] via memset
      0.003952         29 <deep_copy>
      0.001552        210 Kokkos::View::initialization [] via memset
      0.001116          2 Kokkos::View::initialization [d_out_var_mirror] via memset
      0.000994          8 SendBuff
      0.000911          2 Kokkos::View::initialization [cc_outvar] via memset
      0.000732          8 RecvBuff
      0.000558          2 radiationbc_x1
      0.000505          1 sc_hflx_x1
      0.000493          3 sc_sph_sigma
      0.000485          3 sc_chi
      0.000411          1 sc_compute_hk
      0.000404          2 Kokkos::ViewCopy-3D
      0.000352          1 sc_compute_j
      0.000253          1 sc_hflx_bdry_x1
      0.000205          2 Kokkos::View::initialization [sc_ir] via memset
      0.000106          1 sc_hflx_bdry_x2
      0.000085          1 sc_hflx_bdry_x3
      0.000061          1 sc_sph_brad
      0.000059          2 radiationbc_x2
      0.000058          2 radiationbc_x3
      0.000047          4 Kokkos::View::initialization [rank_sendbuf_vars] via memset
      0.000042          4 Kokkos::View::initialization [rank_recvbuf_vars] via memset
```

## scpr3 note on the PR timing arms
- Every PR run ends after one cycle with `Terminating on time limit`, even with tlim = 1e30 named (`tlim=1.000000e+30 nlim=10` is in the log). So the 10-cycle s/cycle does not exist; the sc_sphere problem stops itself after its single step.
- It also prints `### WARNING: SC iteration did not converge in 1 sweeps (max residual = 1.0, iter_tol = 1e-12)`, as expected with iter_max = 1.
- Use the kernel timer instead. Per run it reports fenced kernel seconds for the whole run (1 sweep + set-up):
  - G1_P1a: sc_sweep_plane 0.0183 s (128 launches); sc_residual_j 0.0001 s (1 launches); sc_carry_j 0.0000 s (1 launches);  kp total 0.594873 s
  - G1_P1b: sc_sweep_plane 0.0183 s (128 launches); sc_residual_j 0.0001 s (1 launches); sc_carry_j 0.0000 s (1 launches);  kp total 0.241192 s
  - G1_P2a: sc_sweep_tiled 0.0110 s (22 launches); sc_residual_j 0.0001 s (1 launches); sc_carry_j 0.0000 s (1 launches);  kp total 0.283824 s
  - G1_P2b: sc_sweep_tiled 0.0110 s (22 launches); sc_residual_j 0.0001 s (1 launches); sc_carry_j 0.0000 s (1 launches);  kp total 0.508206 s
  - G4_P1: sc_sweep_plane_byaxis 0.0070 s (256 launches); sc_residual_j 0.0000 s (1 launches); sc_carry_j 0.0000 s (1 launches);  kp total 0.635824 s
- So one PR sweep of 168 rays at 128^3 costs:
  - 1 H200: 18.3 ms with the plane kernel (128 launches), 11.0 ms with tiled-16 (22 launches);
  - 4 H200 (byaxis kernel, 256 launches): 7.0 ms on rank 0.
- For comparison, our 27zg 1-GPU 128^3 Cartesian arm (160 rays): `vet_sc ... 0.154 s per call`, and sp sweeps 1.91 s over 7 calls (0.27 s per call).
  - These are not like for like: one of our calls is the full formal solution plus its sums.
  - The kp.0 of our 27zg G1 arms (posted in 27ze) split it per kernel.
