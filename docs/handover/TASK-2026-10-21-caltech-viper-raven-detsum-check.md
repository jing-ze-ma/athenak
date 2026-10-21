# TASK 2026-10-21 caltech -> viper + raven: check the FINAL reproducible setup on your machines (USER)

Caltech found and fixed why RSG E / HG global runs were NOT reproducible (C2 != C from cycle ~57) and restarts were NOT bitwise.

## Root cause
Floating-point MPI_Allreduce(MPI_SUM) in the evolution path, above all the adaptive MLT closure's shell sums in
he_star_m1 HsAdaptiveUpdate (pgen/he_star_m1.cpp ~3767). With Mellanox HCOLL (default in HPC-X Open MPI) the summation
order can change from run to run and after a restart, so the last bits differ, and the convection amplifies that.
Proof on 1 x 4 H200: OMPI_MCA_coll_hcoll_enable=0 alone -> C2 == C (121 rows) and R == C (cycles 101-150).

## Code fix (please review and merge into fin51-r-fix)
origin/rst-detsum-1021 = f55be371 on top of origin/vetgd-xchg-fast-1021 63a7d28b (faster cs vet_gd exchange, bitwise = 4de8168d):
- 0e22e91a: he_star_m1 adaptive MLT closure sums its shell sums over ranks in rank order (MPI_Allgather + local sum).
- f55be371: M1GlobalSum2 / M1GlobalSumArr / M1GlobalBcg always rank-ordered (MPI_Allgather) on > 1 rank.
Round-off order changes, so results are NOT bitwise vs older binaries; physics unchanged.

## Caltech validation (deck rsg_E_det N=32, nlim 200, rst at cycle 100; C, C2 from t=0, R from C's rst)
2 nodes x 2 H200, HCOLL ON, UCX GPUDirect on + per-GPU NIC pinning: C2==C (201 hst rows + final bins), R==C (100 rows + bins) PASS.
Also PASS: 63a7d28b + HCOLL off + GPUDirect/pinning. Speed: f55be371 + HCOLL on is ~5 % faster than HCOLL off.

## Please do on Viper (MI300A) and Raven (A100), on your usual RSG E multi-node layout
1. Build f55be371 (he_star_m1, your usual flags).
2. Run the same triplet with your production env and HCOLL / collective library left at its default: C (t=0, nlim 200), C2 (same),
   R (restart from C's cycle-100 rst). Report: C2==C? R==C? (all hst rows + the last hydro_w/m1 bins), and s/cycle vs your current binary.
3. Also run one C/C2 pair with your CURRENT production binary, to show whether your site was affected (Caltech: yes, from cycle ~57).
4. NOTE the result on bsg-files-1009.

## Caltech production from now on (both started fresh)
RSG E (rsg_E_1020/prodF, 5240541-43) and HG global (hgglob_1021/prodF, 5240544-46): f55be371 (md5 50772674), HCOLL on,
UCX_IB_GPU_DIRECT_RDMA=yes, each rank's UCX_NET_DEVICES = its GPU's nearest IB port (Caltech H200 nodes: GPU0 mlx5_3, GPU1 mlx5_2,
GPU2 mlx5_5, GPU3 mlx5_4 per nvidia-smi topo; without pinning GPUDirect crashed with "Transport retry count exceeded").

## Addendum: floating-point MPI_SUM Allreduces NOT converted (judged diagnostic / counters / off in the RSG E deck)
rad_m1_opcheck.cpp (6); rad_m1_precond.cpp 785, 871, 1250, 1466 (mg_gc/gf); rad_m1_sph.cpp 715; rad_m1_vet.cpp 3639-3640;
rad_m1_vetgd_cs.cpp 2168-2169; rad_m1_vetgd.cpp 1613, 5946, 5949, 6503, 6514 (+ non-twin_det branches 6999, 8088);
rad_m1_implicit.cpp 9644, 9665, 9674, 12216, 15228, 15470, 16081, 17068, 17341, 17424; rad_m1_time2.cpp 1390;
hydro_mergelat.cpp 440; dt_retry.cpp 631; turb_driver.cpp (4); he_star_m1.cpp 3052, 3077, 4082-4083, 4308-4309.
Any that feeds the evolution in another deck (e.g. HG global: vet_col+vet_gd half blend, hesdirk2) needs the same HsDetSum treatment;
please run the triplet on the HG global deck as well. Debug tooling (state hashes per cycle/task): origin rsg-rst-dbg-1021.
