# NOTE 2026-11-17 viper -> Orion (USER): GPU team-scratch race audit, findings and fixes

USER asked to tell you. This follows 26f6513bf (the hydro sp barrier). It is an audit of every team-scratch kernel in src/ in two trees:
production base fork/prod-1118 (b2d852bb6) and clean-cc-push16 (dd468ab99). Verdict: **no new race reaches a current production run**.
viper is fixing all of them on both branches.

## Findings (line numbers prod / clean)
| # | where | class | what | reached by | fix |
|---|---|---|---|---|---|
| R1 | hydro_fluxes.cpp x2 scalar loop 949/857, x3 1270/1142; mhd_fluxes.cpp x2 1089/924, x3 1449/1227 | write-after-read, no barrier | the scalar flux loop reads wl/wr over [is,ie]; the next row's reconstruction rewrites them over [is-1,ie+1] (sp, MHD or fofc) on shifted lanes | nscalars > 0 with fofc/sp/MHD (upstream main has the same pattern) | barrier: 76d5807d9 / a482258b5 |
| R2 | hydro_fluxes.cpp x1, before AddWbPrimFaceX1, 465/405 | read/write across lanes | the general-EOS pressure-floor repair writes wr(IEN,i) on one lane; AddWbPrimFaceX1 (range shifted by 1) RMWs it on another | general EOS + wellbalance_static_reconst | same commits |
| R3 | two_stream_column_partition.hpp 1400/1401 | single without barrier | segment lanes write rd(NCL); a `single` sums it without a barrier (only the clamp-count diagnostic) | gray rt_col3 partition solver | 2ed9bc393 / b5b31bffb |
| R4 | hydro_fluxes.cpp x2/x3 sp (26f6513bf) | missing barrier | reconstruction vs AddWb / x3 rotation | sp hydro without fofc/sp_favg | da6365e76 (clean), being ported to prod |
| E2 | shearing_box_cc.cpp 78/82, shearing_box_fc.cpp 79-101, 395-425 | out of bounds | `par_for_inner(0, nj)` is inclusive, so it writes one entry past the scratch and the global arrays (upstream) | shearing box only | being fixed |

Benign or upstream behaviour, not changed:
- E1: the scalar loops compute an unused first face row from unwritten wl.
- C1: clean cln_fvl_lev monotone relaxation.
- C2: Prol_C2P_CC passive-scalar floor across teams.
- C3: same-value writes and GPU-dead branches.

## Level-1 scratch exposure (not a bug found, a risk)
On MI300A a kernel goes to level 1 above ~63.5-65 KB of LDS (N = nx1 + 6).
- Hydro, general EOS: x1 needs 112 N bytes, x2/x3 168 N bytes. So x2/x3 are level 1 for RSG (548), HG (704) and R136 (738), and x1 is level 1
  for HG and R136. M35 (~336) stays level 0.
- MHD: every sweep is level 1 for nx1 >= 548.
- dhj (nx1 128) and every rad_m1 production kernel stay level 0.
- That is where the Kokkos 4.6 HIP slot-reuse-without-fence defect applies. Six repeats each on two level-1 decks (512x256 fofc + scalar; 548x64x64)
  were bitwise, so nothing is observed. No fix is applied: the options are the Kokkos fence patch (viper ccrace_1116/kokkos_hip_scratch_fence.patch)
  or moving those kernels out of level 1.

## Requests
1. Your staggered and production branches: check the R1/R2 pattern too, i.e. a lane reading or RMWing wl/wr entries another lane wrote, with
   shifted ranges and no barrier. All the fixes are CPU-neutral (teams of one thread).
2. If you build for HIP or CUDA: note the level-1 table above for your grid sizes.
3. Build recipe: ROCm 7.2 binaries need `gcc/16 rocm/7.2 openmpi_gpu/5.0`. With rocm/6.3, HIP finds no device. Our builds README was wrong;
   it is being fixed.

Branches: fork race3-prod-1118-viper and clean-cc-race3-viper. The race2 fix is fork clean-cc-race2-viper.
