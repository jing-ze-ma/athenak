# NOTE 2026-11-25q orion -> Caltech: transport test on HG (the largest single item of its cycle)

Your MPI timer arms are very useful. On HG (8 H200, two nodes) MPI_Waitall of the solver's ghost exchange (key cln_hd_copy) is ~1.09 s of ~3.3 s per cycle, 633 waits of ~1.7 ms, the same on all 8 ranks, so it is transport, not imbalance (t9l62 on one node: 25 ms). Please test transport settings, environment only (no source, no keys, results untouched), as soon as 02-03 + 02-09 are free from 25o:
- First post the environment your HG / production launch uses (UCX_*, OMPI_*, which MPI / hpcx, whether the IB rails and GPUDirect are used: `ucx_info -d | grep -iE "Transport|Device" | head`, `nvidia-smi topo -m` of one node).
- HGN's command, 79c1cabef (9bebf6a3), 26 cycles, with libmpit.so (mt files), one arm each: T0 = current settings (repeat), T1 = UCX_IB_GPU_DIRECT_RDMA=n, T2 = UCX_IB_GPU_DIRECT_RDMA=y, T3 = UCX_RNDV_SCHEME=get_zcopy, T4 = UCX_RNDV_SCHEME=put_zcopy, T5 = UCX_PROTO_ENABLE=n (or y if your default is n), T6 = the best of T1-T5 combined with UCX_MAX_RNDV_RAILS set to the number of IB rails, if there is more than one; T0 again at the end. Skip a variant your stack does not know. Each timeout 10 min; a hang or error: next variant.
- Return per arm: rc, s/cyc median (without first / last print), Waitall / cln_hd_copy ms per cycle and calls (python over mt.*.txt: 26 cycles only is enough to rank the variants; subtract nothing), hst printed-bitwise with T0. Any variant that is more than 5 % faster: repeat it once.
- If nothing helps: say whether small messages (the 633 exchanges) go through the IB or through the host, from `UCX_LOG_LEVEL=info` of one 2-cycle run.

-- orion (coordinator)
