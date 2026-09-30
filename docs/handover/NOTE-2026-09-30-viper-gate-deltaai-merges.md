# NOTE 2026-09-30 (viper): gate of the DeltaAI merges + the FOFC inline-lambda fix

**Result: all gates passed on viper.** Merged: fix-fofc-inline-lambda (4d4b0e0b + cba4e797) into rt-integration.

## What was gated
1. **CUDA c2p_track fix** (952abbfe, merged ad992f7b): `c2p_track::Before/After` moved from the header into `src/utils/c2p_track.cpp`.
2. **Cubed-sphere seam flux send waits** (2dfe4e95, merged 793e03c3): ClearSend/ClearRecv also wait on the seam requests on a uniform cubed sphere.
3. **This fix** (branch `fix-fofc-inline-lambda`): the same pattern as the c2p_track fix, applied to:
   - `FofcTrialRemoveGrav`: moved from `src/hydro/fofc_etotgrav.hpp` into the new `src/hydro/fofc_etotgrav.cpp`. It is called from hydro_fofc.cpp and mhd_fofc.cpp. Under CUDA, an MHD run with fofc and etotgrav on would probably crash the way c2p_track did.
   - `pgen_eos::HostGamma1FromP`: moved from `src/pgen/pgen_eos_utils.hpp` into the new `src/pgen/pgen_eos_utils.cpp`. It is called from linear_wave.cpp and box_convection.cpp (2 TUs in a box_convection build).

   Nothing else in the code changed. The header `FofcGravDiv` stays `KOKKOS_INLINE_FUNCTION` because it contains no lambda.

## Inline-lambda survey (src/, headers containing KOKKOS_LAMBDA, transitive includes)
- **Fixed** (non-template inline, used in 2 or more TUs): the two functions listed in item 3.
- **Single TU per build, not changed:** these headers are included only by the user pgen (dhj, box or red_giant):
  - two_stream_rt.hpp: CkConserveFlux, CkNonconvLoc, CkNcRefCmp, TsrtCkChain, picket_fence_two_stream_RT_pass, CkCadStore, CkCadLin
  - correlated_k.hpp: ck_rt_selftest, ck_selftest, ck_build_rosseland_table
  - two_stream_column_*.hpp: RTCol3TeamLaunch, RTCol3Launch, CkWarmSeed, CkImplStep
  - utils/runaway_scan.hpp:75 Scan, included only by hydro_tasks.cpp. If the MHD tasks ever include it, it needs the same fix.
- **Templates, not changed:** the risk under nvcc is unclear.
  - athena.hpp par_for / par_for_outer; rad_m1_parfor.hpp par_for_lb; two_stream_rt.hpp par_reduce_clip3/4 and TsrtGreyRt; two_stream_column_ck.hpp CkParFor4. These are wrappers that take the caller's lambda.
  - multigrid.hpp:272 CalculateDefect, :294 CalculateFASRHS (member templates defined in the class) and :622 Multigrid::Smooth. These member templates are included in 13 TUs; their instantiations may be shared across TUs.
  - eos/primitive_solver_hyd.hpp:236, :359 (class template PrimitiveSolverHydro; GR/dyngrmhd family).

## Gates (viper, ROCm 7.2, apudev 1 node x 2 GPUs, HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1)
**Binaries** (`/viper/ptmp2/jinma/builds/bin/`):
- pre-merge: `athena_{dhj,box}_gpu72_a41c901e`. Its src/ is identical to 9a4a7d62, the commit just before the DeltaAI merges.
- DeltaAI tip: `*_6ee19692`.
- tip + this fix: `*_cba4e797`.

**Runs** are in `/viper/ptmp2/jinma/gate_deltaai_0930`: jobs 12031857/8/9 and 12032284; `compare.sh` does the comparison.

| gate | test | result |
|---|---|---|
| a | dhj hydro (WASP-121b 1x), restart w121prod_0929/w1x/rst/dhj.00600.rst (rot 300), 107 cycles, 2 ranks | tip and fix vs pre: hydro.hst, final rst, cycle/dt lines **bitwise identical**; 0 FATAL/nan |
| b | He box (bench-2026-09-29-hebox input, fofc + etotgrav on, so it calls FofcTrialRemoveGrav in hydro), 30 cycles | hydro.hst, user.hst, rst 00000/00001, cycle lines **bitwise identical** |
| c | dhj MHD (w121_mhd_0929 package, bbot 3 G = 0.846 code, max_eta 1e13, lhlld), fresh start, 100 cycles; fofc 0 and fofc 1 at nghost 3 (calls FofcTrialRemoveGrav in MHD::FOFC) | clean (0 FATAL/nan); mhd.hst, rst 00601/00602, cycle lines **bitwise identical** in both arms |
| d | tst CPU on the fix branch: hydro fofc blast/blast_mpi/cs/sod/sp/wb, restart_bitwise; mhd balsara_vortex, fofc blast/blast_mpi/cs, lhlld_lowmach; rad dhj_ck, dhj_ck_mpi, dhj_srclim | **15/15 pass**. Style: the new and changed files are cpplint-clean (the pre-existing tree-wide count is unchanged) |

**Notes on the gates:**
- dhj user.hst columns 12-14 (Efloor, Mfloor, Efloor_rt) differ at 1e-15 relative. This is run-to-run noise: a repeat of the pre-merge binary also differs in user.hst while its rst is bitwise identical. The floor-diagnostic reductions are not deterministic.
- dhj with `fofc = 1` at the production nghost 2 stops at startup with "FOFC and plm reconstruction requires at least 3 ghost zones", in every binary.
- HIP was never affected by the inline-lambda problem. This gate only shows the fixes change nothing on viper. The CUDA crash itself is fixed only by construction; a CUDA MHD run with fofc on, on DeltaAI or Caltech, would confirm it.
