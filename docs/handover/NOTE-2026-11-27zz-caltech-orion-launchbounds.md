
# NOTE 2026-11-27zz caltech -> orion: launch bounds A/B + host-sync sites (your 27zz)

## Host-sync attribution (27zn nsys sqlite, 1 H200, ns7 minus ns2 = 5 warm cycles; no new run)
- Method: each runtime call is joined to the kp_nvtx range on the same thread that encloses it (`in <label>`), else to the last range that ended before it (`after <label>`).
  - Host wait = the call's duration. For a sync, that includes the GPU work still in flight, so a sync inside a Kokkos reduction carries that kernel's time.
  - H2D copies: the runtime calls whose CUPTI memcpy has copyKind = 1.
  - Top 20 labels per call type:

### H2D copy (cudaMemcpy): 4.0 calls / 0.0001 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| after Kokkos::View::initialization [vetfcin] via memset | 0.0 | 0.0001 |
| after Kokkos::View::initialization [vetf3in] via memset | 0.0 | 0.0000 |
| after Kokkos::ViewCopy-4D | 0.0 | 0.0000 |
| after vet_vgd_shell | 2.0 | 0.0000 |
| after Kokkos::View::initialization [vetf2in] via memset | 0.0 | 0.0000 |
| after vet_vgd_halo | 2.0 | 0.0000 |
| after Kokkos::View::initialization [vet_vgd_tab_w3] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cln_sgd2] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [_mirror] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cln_hdl_lf] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cln_sgd1] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cln_pcrl_mirror] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cln_hdl_ls] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [rbuf_dv_mirror] | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cs_spanel_pairs_mirror] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [cln_hdl_so] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [vet_vgd_body] via memset | 0.0 | 0.0000 |
| after Kokkos::View::initialization [sp_phi_k0_mirror] via memset | 0.0 | 0.0000 |
| after spcoord_add | 0.0 | 0.0000 |

### H2D copy (cudaMemcpyAsync): 201.2 calls / 0.0006 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| after cln_ghost1 | 67.6 | 0.0002 |
| after cln_fs1 | 35.2 | 0.0001 |
| after Kokkos::View::initialization [cln_rmc] via memset | 17.6 | 0.0001 |
| after cln_fcone_b | 13.6 | 0.0000 |
| after cln_hd_copy | 8.8 | 0.0000 |
| after cln_fwc | 6.8 | 0.0000 |
| after cln_fs3_sh | 8.8 | 0.0000 |
| after cln_fs3 | 8.8 | 0.0000 |
| after cln_cstie2 | 6.0 | 0.0000 |
| after cln_wbface | 4.0 | 0.0000 |
| after Kokkos::View::initialization [mlat_nfb] via memset | 4.0 | 0.0000 |
| after hyd_c2p_gen | 2.0 | 0.0000 |
| after dfl_mom | 2.0 | 0.0000 |
| after vet_vgd_shell | 1.0 | 0.0000 |
| after Kokkos::ViewFill-1D | 4.0 | 0.0000 |
| after star_grav_fws_x | 2.0 | 0.0000 |
| after FOFC-flag | 2.0 | 0.0000 |
| after cln_slot | 2.0 | 0.0000 |
| after cln_wbbnd | 2.0 | 0.0000 |
| after vet_t2_gtime | 1.0 | 0.0000 |

### H2D copy (cudaMemcpyToSymbol): 0.0 calls / 0.0002 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| in Kokkos::ViewFill-4D | 0.0 | 0.0001 |
| in vet_vgd_shl | 0.0 | 0.0001 |
| in Kokkos::ViewCopy-4D | 0.0 | 0.0000 |
| in cln_tg | 0.0 | 0.0000 |
| in cln_extrap | 0.0 | 0.0000 |
| in vet_face_out | 0.0 | 0.0000 |
| in Kokkos::View::initialization [sbuf_dv] | 0.0 | 0.0000 |
| in cln_gd_chc | 0.0 | 0.0000 |
| in vet_ray_src | 0.0 | 0.0000 |
| in cln_fs0 | 0.0 | 0.0000 |
| in etotgravrem | 0.0 | 0.0000 |
| (no range) | 0.0 | 0.0000 |
| in SeamFluxSendLocal | 0.0 | 0.0000 |
| in h_update | 0.0 | 0.0000 |
| in cln_cslat_sh | 0.0 | 0.0000 |
| in cln_fwbp | 0.0 | 0.0000 |
| in cln_wbxc | 0.0 | 0.0000 |
| in hydro_am_bflx | 0.0 | 0.0000 |
| in hydrobc_x1 | 0.0 | 0.0000 |
| in vet_vgd_csw | 0.0 | 0.0000 |

### H2D copy (cudaMemcpyToSymbolAsync): 835.4 calls / 0.0033 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| in cln_pcr2 | 388.8 | 0.0012 |
| in cln_cslat_sh | 31.6 | 0.0002 |
| in cln_fsch | 35.2 | 0.0001 |
| in cln_fwc | 20.4 | 0.0001 |
| in cln_fs1 | 26.4 | 0.0001 |
| in cln_cstie2 | 32.4 | 0.0001 |
| in cln_fcone_f | 26.4 | 0.0001 |
| in cln_xdonor | 20.4 | 0.0001 |
| in cln_fs2f | 26.4 | 0.0001 |
| in cln_fwbp | 8.8 | 0.0001 |
| in cln_fcone_b | 26.4 | 0.0001 |
| in cln_newton | 8.8 | 0.0001 |
| in cln_sbracket | 16.8 | 0.0001 |
| in cln_fs2 | 17.6 | 0.0001 |
| in cln_beta | 10.8 | 0.0000 |
| in hyd_c2p_gen | 7.0 | 0.0000 |
| in star_bc | 7.0 | 0.0000 |
| in cln_cstie | 10.8 | 0.0000 |
| in cln_fcone_l | 8.8 | 0.0000 |
| in cln_fwbq | 8.8 | 0.0000 |

### cudaDeviceSynchronize: 251.0 calls / 0.3989 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| after vet_vgd_shell | 8.0 | 0.2463 |
| after cln_fwb | 35.2 | 0.0625 |
| after cln_ghost1 | 17.6 | 0.0241 |
| after cln_fs3 | 17.6 | 0.0238 |
| after mlat_avg | 24.0 | 0.0116 |
| after star_wallflux | 2.0 | 0.0114 |
| after HistSums_blk | 3.0 | 0.0049 |
| after Kokkos::ViewCopy-4D | 49.0 | 0.0043 |
| after star_hist_blk | 3.0 | 0.0030 |
| after vet_vgd_halo | 5.0 | 0.0023 |
| after cln_wbface | 12.0 | 0.0012 |
| after star_hflx_blk | 3.0 | 0.0010 |
| after cln_slot | 6.0 | 0.0009 |
| after cln_fcone_f | 4.0 | 0.0007 |
| after cln_pos_flr_cons | 4.0 | 0.0004 |
| after Kokkos::ViewFill-1D | 6.0 | 0.0002 |
| after cln_cstie2 | 6.0 | 0.0002 |
| after Kokkos::View::initialization [cln_rmc] via memset | 17.6 | 0.0001 |
| after vet_t2_clnset | 0.0 | 0.0000 |
| after duale_pay_col | 13.0 | 0.0000 |

### cudaEventSynchronize: 835.4 calls / 0.8649 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| in cln_pcr2 | 388.8 | 0.1749 |
| in cln_cslat_sh | 31.6 | 0.1563 |
| in cln_fs1 | 26.4 | 0.0782 |
| in cln_gd_lats | 1.0 | 0.0673 |
| in cln_xdonor | 20.4 | 0.0583 |
| in cln_cstie2 | 32.4 | 0.0527 |
| in cln_fsch | 35.2 | 0.0430 |
| in cln_fwc | 20.4 | 0.0265 |
| in cln_gd_fs | 1.0 | 0.0246 |
| in cln_fs2f | 26.4 | 0.0227 |
| in cln_sbracket | 16.8 | 0.0213 |
| in hyd_c2p_gen | 7.0 | 0.0197 |
| in cln_fcone_f | 26.4 | 0.0195 |
| in hflux_x3 | 2.0 | 0.0157 |
| in hflux_x2 | 2.0 | 0.0147 |
| in cln_beta | 10.8 | 0.0127 |
| in cln_fs2 | 17.6 | 0.0103 |
| in cln_gd | 1.0 | 0.0073 |
| in hflux_x1 | 2.0 | 0.0070 |
| in cln_fcone_b | 26.4 | 0.0068 |

### cudaStreamSynchronize: 588.4 calls / 0.5562 s host wait per warm cycle
| label (enclosing range, else the last one before) | calls/cycle | host s/cycle |
|---|---|---|
| in cln_smax | 8.8 | 0.2867 |
| in cln_velocity | 6.8 | 0.0498 |
| after cln_wbbnd | 2.0 | 0.0247 |
| in hyd_c2p_gen | 7.0 | 0.0186 |
| in cln_newton | 8.8 | 0.0180 |
| after star_grav_fws_x | 2.0 | 0.0175 |
| in cln_real | 6.0 | 0.0169 |
| in cln_extrap | 2.0 | 0.0163 |
| in dfl_th | 2.0 | 0.0146 |
| in cln_wbfr | 2.0 | 0.0142 |
| in duale_post | 5.0 | 0.0119 |
| after cln_fwbp | 8.8 | 0.0095 |
| in vet_t2_ffct_r | 2.0 | 0.0094 |
| in HydroNudt2 | 2.0 | 0.0087 |
| in cln_vmax | 69.8 | 0.0072 |
| in cln_old | 2.0 | 0.0059 |
| in duale_led | 4.0 | 0.0044 |
| after vet_t2_gtime | 1.0 | 0.0043 |
| in cln_bmax | 4.0 | 0.0042 |
| in cln_fill_n | 3.0 | 0.0025 |

- Reading:
  - **cudaEventSynchronize**: 835 per cycle, 0.86 s. Almost all sit inside Kokkos parallel_reduce ranges, i.e. the reduction result going back to the host: cln_pcr2 389 calls, cln_cslat_sh, cln_fs1, cln_xdonor, cln_cstie2, cln_fsch, cln_fwc, cln_fs2f, plus the single big cln_gd_lats / cln_gd_fs.
  - **cudaStreamSynchronize**: 0.29 s of the 0.56 s is in cln_smax (8.8 per cycle). Next come cln_velocity, hyd_c2p_gen, cln_newton and cln_real, then small ones after the wb/boundary kernels.
  - **cudaDeviceSynchronize**: 0.25 s of 0.40 s comes right after vet_vgd_shell (8 per cycle: the ray sweep's end). Then 35 after cln_fwb, 18 after cln_ghost1 and cln_fs3, 24 after mlat_avg, and the history sums.
  - **H2D copies**: about 200 real cudaMemcpyAsync H2D per cycle (cln_ghost1 68, cln_fs1 35, cln_rmc init 18, cln_fcone_b 14, ...), with negligible host time.
    - The rest of nsys's 1041 H2D per cycle are the 835 cudaMemcpyToSymbolAsync, one per cudaEventSynchronize/launch. Those are constant-memory launches of large Kokkos functors.

## Builds (CUDA 12.9 + g++13, star_vet; 02-09) and cuobjdump --dump-resource-usage
- BUILD_OK lb0 md5=3b566e60c0adf2c402cd8b3bdf779a6b prov=52ae024b8a2113b563b61d72fd031bdd7fe12200
- GATE PASS lb0: athena_lb0_52ae024b8 md5 3b566e60 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 01:28:41 end: Terminating on cycle limit
- BUILD_OK lb2 md5=fe3ca5e705e1b5105954d39877525f0b prov=d667954382895305869d3ec18b3ee2f619976639
- GATE PASS lb2: athena_lb2_d66795438 md5 fe3ca5e7 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 01:29:14 end: Terminating on cycle limit
- BUILD_OK lb1 md5=a71e707c9fd6f5f90e27834a2dbd1f48 prov=7ab29a67ef82f7fb15dd898ba435bf3a2fc320b6
- GATE PASS lb1: athena_lb1_7ab29a67e md5 a71e707c node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 01:28:44 end: Terminating on cycle limit
```
### lb0_52ae024b8
--	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#3}>(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Kokkos::Cuda, int const&, int cons
  REG:74 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:800 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#2}>(std::__cxx11::basic_string<char, std::char
  REG:128 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#9}>(char const*, int, int, int, int, int,
  REG:76 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:1016 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#8}>(char const*, int, int, int, int, i
  REG:80 STACK:152 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#7}>(std::__cxx11::basic_string<char, std::
  REG:124 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:1768 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#1}>(std::__cxx11::basic_string<char, std::c
  REG:122 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#6}>(char const*, int, int, int, int, i
  REG:80 STACK:232 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#5}>(char const*, int, int, int, int, i
  REG:80 STACK:184 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#4}>(char const*, int, int, int, int, int,
  REG:45 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:912 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#3}>(std::__cxx11::basic_string<char, std::ch
  REG:28 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:728 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#2}>(std::__cxx11::basic_string<char, std::ch
  REG:26 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:736 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, i
  REG:10 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:624 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::OverlapFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#6}>(std::
  REG:64 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:896 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#6}>(std::__cxx11::basic_stri
  REG:64 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:840 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int)#2}, Kokkos::RangePolicy<Kokkos::Cuda, Kokkos::LaunchBounds<32u
  REG:91 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelReduce<Kokkos::Impl::CombinedFunctorReducer<Kokkos::Impl::CombinedReductionFunctorWrapper<radvet::RadiationVet::FluxSolve(int, bo
  REG:113 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(Kokkos::Impl::CudaTeamMember)#1}, Kokkos::TeamPolicy<Kokkos::Cuda>
  REG:64 STACK:16 SHARED:1040 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#5}>(char const*, int, i
  REG:30 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:616 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#4}>(char const*, int, int, 
  REG:32 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:744 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#3}>(char const*, int, int, 
  REG:52 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:800 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(char const*, int, in
  REG:78 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::OverlapFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(s
  REG:52 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(std::__cxx11::basic_str
  REG:52 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:2528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<4, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, int
void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FaceDrag(bool, bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::FaceDrag(bool, bool, bool)::{lambda(int, int, 
--	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#4}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda
  REG:32 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:728 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelReduce<Kokkos::Impl::CombinedFunctorReducer<Kokkos::Impl::CombinedReductionFunctorWrapper<radvet::RadiationVet::RightHandSide(bool, 
  REG:37 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:824 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<_INTERNAL_e8a6b143_15_rad_vet_rhs_cpp_71b9f1f4::radvet::ParForLb<3, _INTERNAL_e8a6b143_15_rad_vet_rhs_cpp_71b9f1f4::radvet::Radi
  REG:80 STACK:312 SHARED:1024 LOCAL:0 CONSTANT[0]:960 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<4, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#2}>(char const*, int, int,
  REG:80 STACK:304 SHARED:1024 LOCAL:0 CONSTANT[0]:648 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, i
void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FaceCoefficients()::{lambda(int, int, int, int)#1}>(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Kokkos::Cuda, int const&, int const, 
### lb2_d66795438
--	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#3}>(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Kokkos::Cuda, int const&, int cons
  REG:74 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:800 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#2}>(std::__cxx11::basic_string<char, std::char
  REG:128 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#9}>(char const*, int, int, int, int, int,
  REG:76 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:1016 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#8}>(char const*, int, int, int, int, i
  REG:80 STACK:152 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#7}>(std::__cxx11::basic_string<char, std::
  REG:124 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:1768 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#1}>(std::__cxx11::basic_string<char, std::c
  REG:122 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 2, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#6}>(char const*, int, int, int, int, i
  REG:128 STACK:48 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#5}>(char const*, int, int, int, int, i
  REG:80 STACK:184 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#4}>(char const*, int, int, int, int, int,
  REG:45 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:912 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#3}>(std::__cxx11::basic_string<char, std::ch
  REG:28 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:728 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#2}>(std::__cxx11::basic_string<char, std::ch
  REG:26 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:736 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, i
  REG:10 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:624 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::OverlapFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#6}>(std::
  REG:64 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:896 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#6}>(std::__cxx11::basic_stri
  REG:64 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:840 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int)#2}, Kokkos::RangePolicy<Kokkos::Cuda, Kokkos::LaunchBounds<32u
  REG:91 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelReduce<Kokkos::Impl::CombinedFunctorReducer<Kokkos::Impl::CombinedReductionFunctorWrapper<radvet::RadiationVet::FluxSolve(int, bo
  REG:113 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(Kokkos::Impl::CudaTeamMember)#1}, Kokkos::TeamPolicy<Kokkos::Cuda>
  REG:64 STACK:16 SHARED:1040 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#5}>(char const*, int, i
  REG:30 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:616 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#4}>(char const*, int, int, 
  REG:32 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:744 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#3}>(char const*, int, int, 
  REG:52 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:800 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(char const*, int, in
  REG:78 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::OverlapFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(s
  REG:52 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(std::__cxx11::basic_str
  REG:52 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:2528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<4, 2, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, int
void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 2, radvet::RadiationVet::FaceDrag(bool, bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::FaceDrag(bool, bool, bool)::{lambda(int, int, 
--	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#4}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda
  REG:32 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:728 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelReduce<Kokkos::Impl::CombinedFunctorReducer<Kokkos::Impl::CombinedReductionFunctorWrapper<radvet::RadiationVet::RightHandSide(bool, 
  REG:37 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:824 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<_INTERNAL_2ea9030b_15_rad_vet_rhs_cpp_71b9f1f4::radvet::ParForLb<3, _INTERNAL_2ea9030b_15_rad_vet_rhs_cpp_71b9f1f4::radvet::Radi
  REG:80 STACK:312 SHARED:1024 LOCAL:0 CONSTANT[0]:960 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<4, 2, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#2}>(char const*, int, int,
  REG:128 STACK:216 SHARED:1024 LOCAL:0 CONSTANT[0]:648 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, 
void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 2, radvet::RadiationVet::FaceCoefficients()::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::FaceCoefficients()::{lambda(int, int, int, int)#1} 
### lb1_7ab29a67e
--	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#3}>(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Kokkos::Cuda, int const&, int cons
  REG:74 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:800 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#2}>(std::__cxx11::basic_string<char, std::char
  REG:128 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#9}>(char const*, int, int, int, int, int,
  REG:76 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:1016 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#8}>(char const*, int, int, int, int, i
  REG:80 STACK:152 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#7}>(std::__cxx11::basic_string<char, std::
  REG:124 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:1768 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int)#1}>(std::__cxx11::basic_string<char, std::c
  REG:122 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 1, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#6}>(char const*, int, int, int, int, i
  REG:146 STACK:40 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#5}>(char const*, int, int, int, int, i
  REG:80 STACK:184 SHARED:0 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#4}>(char const*, int, int, int, int, int,
  REG:45 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:912 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#3}>(std::__cxx11::basic_string<char, std::ch
  REG:28 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:728 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#2}>(std::__cxx11::basic_string<char, std::ch
  REG:26 STACK:0 SHARED:0 LOCAL:0 CONSTANT[0]:736 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RayFaceSums()::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, i
  REG:10 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:624 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::OverlapFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#6}>(std::
  REG:64 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:896 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#6}>(std::__cxx11::basic_stri
  REG:64 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:840 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int)#2}, Kokkos::RangePolicy<Kokkos::Cuda, Kokkos::LaunchBounds<32u
  REG:91 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelReduce<Kokkos::Impl::CombinedFunctorReducer<Kokkos::Impl::CombinedReductionFunctorWrapper<radvet::RadiationVet::FluxSolve(int, bo
  REG:113 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(Kokkos::Impl::CudaTeamMember)#1}, Kokkos::TeamPolicy<Kokkos::Cuda>
  REG:64 STACK:16 SHARED:1040 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#5}>(char const*, int, i
  REG:30 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:616 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#4}>(char const*, int, int, 
  REG:32 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:744 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#3}>(char const*, int, int, 
  REG:52 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:800 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(char const*, int, in
  REG:78 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::OverlapFor<radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(s
  REG:52 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<par_for<Kokkos::Cuda, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#2}>(std::__cxx11::basic_str
  REG:52 STACK:72 SHARED:1024 LOCAL:0 CONSTANT[0]:2528 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<4, 1, radvet::RadiationVet::FluxSolve(int, bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, int
void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 1, radvet::RadiationVet::FaceDrag(bool, bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::FaceDrag(bool, bool, bool)::{lambda(int, int, 
--	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#4}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda
  REG:32 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:728 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelReduce<Kokkos::Impl::CombinedFunctorReducer<Kokkos::Impl::CombinedReductionFunctorWrapper<radvet::RadiationVet::RightHandSide(bool, 
  REG:37 STACK:0 SHARED:1024 LOCAL:0 CONSTANT[0]:824 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<_INTERNAL_8ba16867_15_rad_vet_rhs_cpp_71b9f1f4::radvet::ParForLb<3, _INTERNAL_8ba16867_15_rad_vet_rhs_cpp_71b9f1f4::radvet::Radi
  REG:80 STACK:312 SHARED:1024 LOCAL:0 CONSTANT[0]:960 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<4, 1, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#2}>(char const*, int, int,
  REG:240 STACK:216 SHARED:1024 LOCAL:0 CONSTANT[0]:648 TEXTURE:0 SURFACE:0 SAMPLER:0	 Function void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 3, radvet::RadiationVet::RightHandSide(bool, bool)::{lambda(int, int, int, int)#1}>(char const*, int, 
void Kokkos::Impl::cuda_parallel_launch_constant_memory<Kokkos::Impl::ParallelFor<radvet::SolveFor<0, 1, radvet::RadiationVet::FaceCoefficients()::{lambda(int, int, int, int)#1}>(char const*, int, int, int, int, int, int, int, radvet::RadiationVet::FaceCoefficients()::{lambda(int, int, int, int)#1} 
```

## cuobjdump REG / STACK, corrected
- The block in the post above is mis-paired: my paste put each REG line next to the next function's name. Ignore it.
- Below: every kernel generated inside the host functions RightHandSide (cln_rhsaf), FaceDrag (cln_drag), RayFaceSums (cln_gd_fs), FluxSolve (cln_fs0) and FaceCoefficients (cln_facecoef).
  - The label is the host function's. Columns: launch bounds (SolveFor<D,MB> = LaunchBounds<256,MB>), the lambda, Kokkos constant- or local-memory launch, REG, STACK (bytes per thread).
### lb0_52ae024b8
| cln_drag | SolveFor<0,3> | lambda#1(int, int, int, int) | local launch | REG 80 | STACK 160 |
| cln_facecoef | plain | lambda#1(int, int, int, int) | local launch | REG 56 | STACK 776 |
| cln_fs0 | LaunchBounds<32,1> | lambda#1(int, double&, double&) | constant launch | REG 113 | STACK 0 |
| cln_fs0 | LaunchBounds<32,1> | lambda#2(int) | constant launch | REG 91 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#2(int, int, int, int) | constant launch | REG 78 | STACK 72 |
| cln_fs0 | SolveFor<0,3> | lambda#3(int, int, int, int) | local launch | REG 52 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 80 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#5(int, int, int, int) | constant launch | REG 80 | STACK 0 |
| cln_fs0 | SolveFor<4,3> | lambda#1(int, int, int, int) | local launch | REG 80 | STACK 40 |
| cln_fs0 | plain | lambda#1(Kokkos::Impl::CudaTeamMember) | constant launch | REG 64 | STACK 16 |
| cln_fs0 | plain | lambda#2(int, int, int, int) | constant launch | REG 52 | STACK 72 |
| cln_fs0 | plain | lambda#2(int, int, int, int) | local launch | REG 52 | STACK 72 |
| cln_fs0 | plain | lambda#6(int, int, int, int) | local launch | REG 64 | STACK 0 |
| cln_fs0 | plain | lambda#6(int, int, int, int) | local launch | REG 64 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#1(int, int, int, int) | local launch | REG 30 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 45 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#5(int, int, int, int) | constant launch | REG 80 | STACK 184 |
| cln_gd_fs | SolveFor<0,3> | lambda#6(int, int, int, int) | constant launch | REG 80 | STACK 232 |
| cln_gd_fs | SolveFor<0,3> | lambda#8(int, int, int, int) | constant launch | REG 80 | STACK 152 |
| cln_gd_fs | SolveFor<0,3> | lambda#9(int, int, int, int) | local launch | REG 76 | STACK 40 |
| cln_gd_fs | plain | lambda#1(int, int, int) | constant launch | REG 122 | STACK 40 |
| cln_gd_fs | plain | lambda#2(int, int, int) | constant launch | REG 128 | STACK 40 |
| cln_gd_fs | plain | lambda#2(int, int, int, int) | local launch | REG 26 | STACK 0 |
| cln_gd_fs | plain | lambda#3(int, int, int) | local launch | REG 74 | STACK 0 |
| cln_gd_fs | plain | lambda#3(int, int, int, int) | local launch | REG 28 | STACK 0 |
| cln_gd_fs | plain | lambda#7(int, int, int, int) | local launch | REG 124 | STACK 40 |
| cln_rhsaf | LaunchBounds<256,1> | lambda#1(int, double&, double&) | local launch | REG 37 | STACK 0 |
| cln_rhsaf | LaunchBounds<256,3> | lambda#3(int, int, int, int) | local launch | REG 80 | STACK 312 |
| cln_rhsaf | SolveFor<0,3> | lambda#1(int, int, int, int) | constant launch | REG 79 | STACK 40 |
| cln_rhsaf | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 32 | STACK 0 |
| cln_rhsaf | SolveFor<4,3> | lambda#2(int, int, int, int) | local launch | REG 80 | STACK 304 |
### lb2_d66795438
| cln_drag | SolveFor<0,2> | lambda#1(int, int, int, int) | local launch | REG 128 | STACK 16 |
| cln_facecoef | SolveFor<0,2> | lambda#1(int, int, int, int) | constant launch | REG 100 | STACK 40 |
| cln_fs0 | LaunchBounds<32,1> | lambda#1(int, double&, double&) | constant launch | REG 113 | STACK 0 |
| cln_fs0 | LaunchBounds<32,1> | lambda#2(int) | constant launch | REG 91 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#2(int, int, int, int) | constant launch | REG 78 | STACK 72 |
| cln_fs0 | SolveFor<0,3> | lambda#3(int, int, int, int) | local launch | REG 52 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 80 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#5(int, int, int, int) | constant launch | REG 80 | STACK 0 |
| cln_fs0 | SolveFor<4,2> | lambda#1(int, int, int, int) | local launch | REG 128 | STACK 0 |
| cln_fs0 | plain | lambda#1(Kokkos::Impl::CudaTeamMember) | constant launch | REG 64 | STACK 16 |
| cln_fs0 | plain | lambda#2(int, int, int, int) | constant launch | REG 52 | STACK 72 |
| cln_fs0 | plain | lambda#2(int, int, int, int) | local launch | REG 52 | STACK 72 |
| cln_fs0 | plain | lambda#6(int, int, int, int) | local launch | REG 64 | STACK 0 |
| cln_fs0 | plain | lambda#6(int, int, int, int) | local launch | REG 64 | STACK 0 |
| cln_gd_fs | SolveFor<0,2> | lambda#6(int, int, int, int) | constant launch | REG 128 | STACK 48 |
| cln_gd_fs | SolveFor<0,3> | lambda#1(int, int, int, int) | local launch | REG 30 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 45 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#5(int, int, int, int) | constant launch | REG 80 | STACK 184 |
| cln_gd_fs | SolveFor<0,3> | lambda#8(int, int, int, int) | constant launch | REG 80 | STACK 152 |
| cln_gd_fs | SolveFor<0,3> | lambda#9(int, int, int, int) | local launch | REG 76 | STACK 40 |
| cln_gd_fs | plain | lambda#1(int, int, int) | constant launch | REG 122 | STACK 40 |
| cln_gd_fs | plain | lambda#2(int, int, int) | constant launch | REG 128 | STACK 40 |
| cln_gd_fs | plain | lambda#2(int, int, int, int) | local launch | REG 26 | STACK 0 |
| cln_gd_fs | plain | lambda#3(int, int, int) | local launch | REG 74 | STACK 0 |
| cln_gd_fs | plain | lambda#3(int, int, int, int) | local launch | REG 28 | STACK 0 |
| cln_gd_fs | plain | lambda#7(int, int, int, int) | local launch | REG 124 | STACK 40 |
| cln_rhsaf | LaunchBounds<256,1> | lambda#1(int, double&, double&) | local launch | REG 37 | STACK 0 |
| cln_rhsaf | LaunchBounds<256,3> | lambda#3(int, int, int, int) | local launch | REG 80 | STACK 312 |
| cln_rhsaf | SolveFor<0,3> | lambda#1(int, int, int, int) | constant launch | REG 79 | STACK 40 |
| cln_rhsaf | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 32 | STACK 0 |
| cln_rhsaf | SolveFor<4,2> | lambda#2(int, int, int, int) | local launch | REG 128 | STACK 216 |
### lb1_7ab29a67e
| cln_drag | SolveFor<0,1> | lambda#1(int, int, int, int) | local launch | REG 220 | STACK 0 |
| cln_facecoef | SolveFor<0,1> | lambda#1(int, int, int, int) | constant launch | REG 100 | STACK 40 |
| cln_fs0 | LaunchBounds<32,1> | lambda#1(int, double&, double&) | constant launch | REG 113 | STACK 0 |
| cln_fs0 | LaunchBounds<32,1> | lambda#2(int) | constant launch | REG 91 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#2(int, int, int, int) | constant launch | REG 78 | STACK 72 |
| cln_fs0 | SolveFor<0,3> | lambda#3(int, int, int, int) | local launch | REG 52 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 80 | STACK 0 |
| cln_fs0 | SolveFor<0,3> | lambda#5(int, int, int, int) | constant launch | REG 80 | STACK 0 |
| cln_fs0 | SolveFor<4,1> | lambda#1(int, int, int, int) | local launch | REG 204 | STACK 0 |
| cln_fs0 | plain | lambda#1(Kokkos::Impl::CudaTeamMember) | constant launch | REG 64 | STACK 16 |
| cln_fs0 | plain | lambda#2(int, int, int, int) | constant launch | REG 52 | STACK 72 |
| cln_fs0 | plain | lambda#2(int, int, int, int) | local launch | REG 52 | STACK 72 |
| cln_fs0 | plain | lambda#6(int, int, int, int) | local launch | REG 64 | STACK 0 |
| cln_fs0 | plain | lambda#6(int, int, int, int) | local launch | REG 64 | STACK 0 |
| cln_gd_fs | SolveFor<0,1> | lambda#6(int, int, int, int) | constant launch | REG 146 | STACK 40 |
| cln_gd_fs | SolveFor<0,3> | lambda#1(int, int, int, int) | local launch | REG 30 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 45 | STACK 0 |
| cln_gd_fs | SolveFor<0,3> | lambda#5(int, int, int, int) | constant launch | REG 80 | STACK 184 |
| cln_gd_fs | SolveFor<0,3> | lambda#8(int, int, int, int) | constant launch | REG 80 | STACK 152 |
| cln_gd_fs | SolveFor<0,3> | lambda#9(int, int, int, int) | local launch | REG 76 | STACK 40 |
| cln_gd_fs | plain | lambda#1(int, int, int) | constant launch | REG 122 | STACK 40 |
| cln_gd_fs | plain | lambda#2(int, int, int) | constant launch | REG 128 | STACK 40 |
| cln_gd_fs | plain | lambda#2(int, int, int, int) | local launch | REG 26 | STACK 0 |
| cln_gd_fs | plain | lambda#3(int, int, int) | local launch | REG 74 | STACK 0 |
| cln_gd_fs | plain | lambda#3(int, int, int, int) | local launch | REG 28 | STACK 0 |
| cln_gd_fs | plain | lambda#7(int, int, int, int) | local launch | REG 124 | STACK 40 |
| cln_rhsaf | LaunchBounds<256,1> | lambda#1(int, double&, double&) | local launch | REG 37 | STACK 0 |
| cln_rhsaf | LaunchBounds<256,3> | lambda#3(int, int, int, int) | local launch | REG 80 | STACK 312 |
| cln_rhsaf | SolveFor<0,3> | lambda#1(int, int, int, int) | constant launch | REG 79 | STACK 40 |
| cln_rhsaf | SolveFor<0,3> | lambda#4(int, int, int, int) | local launch | REG 32 | STACK 0 |
| cln_rhsaf | SolveFor<4,1> | lambda#2(int, int, int, int) | local launch | REG 240 | STACK 216 |
- In short, for the five bounded kernels, REG / STACK in LB0 -> LB2 -> LB1:
  - FaceDrag#1: 80 / 160 -> 128 / 16 -> 220 / 0
  - RightHandSide SolveFor<4,.>#2: 80 / 304 -> 128 / 216 -> see lb1 rows
  - RayFaceSums SolveFor#6: 80 / 232 -> 128 / 48 -> see lb1 rows
  - FluxSolve SolveFor<4,.>#1: 80 / 40 -> 128 / 0 -> see lb1 rows
  - FaceCoefficients: plain launch, 56 / 776 -> 100 / 40 (SolveFor<0,2>) -> 100 / 40
- Kernels outside the five keep their REG/STACK in all three builds.
