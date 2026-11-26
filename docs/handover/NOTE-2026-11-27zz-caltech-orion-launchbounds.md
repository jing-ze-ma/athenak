
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
