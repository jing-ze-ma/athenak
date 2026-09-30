---
name: m1-perf-merged-0929
description: m1-perf-0928 (implicit_opac_newton DEFAULT ON, He box -58..-65 %) merged into rt-integration d26b7364 on 09-29 (user go)
metadata:
  type: project
---
Merged on viper 09-29 ~08:45 as d26b7364 (pushed to fork), after the viper HIP check (NOTE-2026-09-29-viper-m1-perf-check.md,
3aaaefc6): Picard 2.03 vs 6.76, He box 48.9 vs 117.5 ms/cycle at 1 GPU, hst agree (1-KE 2.7e-7). src after the merge is
identical to the tested branch (rt-integration had only docs since 401875f0).
Effect: <rad_m1>/implicit_opac_newton defaults ON wherever implicit_opac_update + implicit_gas_newton (He box, sp wedge, MHD).
Restarts lacking the key keep it off. Recheck the He presn inner-wall instability with it on ([[he-presn-m1-wedge-0929]]);
he-presn-m1 is based before the merge and must be rebased/merged onto d26b7364.

one_pass on top of Newton (viper gate c, NOTE-2026-09-29-viper-onepass-cost.md): -6/-8.5 % at cfl 0.3, 0 accepted at 0.9 (+2 % re-probes). USER 09-29: stays OPT-IN.

**09-29 ~23:45 DEFECT FOUND in the merged default:** the Newton face term can make the implicit row diagonal NEGATIVE (non-M-matrix) where face flux is large and d kappa/dT has the wrong sign (He presn wedge FeCZ plume); caused the negative-E stage failures there. Guard fix proposed; stopgap implicit_opac_newton = false. See [[he-presn-m1-wedge-0929]].

**GUARD MERGED 09-30 ~00:50 (23ec1fb8, pushed; NOTE-2026-09-29-m1-opn-guard.md):** implicit_opac_newton_guard 0.5 + guard_mode 2 (diagonal test + neighbour entry sign) default; guard 0 = old rows bitwise. He box: 0 drops, bitwise, speed kept (39.6-40.1 ms both). Wedge (rst 00016 -> 19450): 0 floor clips, 0 NC, Picard 8.1 (Newton off 10.3, unguarded 19.4), 8 drops/solve; t 19000 vs Newton-off 5e-10. Fix fast-forwarded onto he-presn-m1 (b4b7c231, binary athena_he_gpu72_0b8b6c0d); He presn input still Newton OFF (c60e9960) pending user.

**09-30 ~02:30 DeltaAI fixes on rt-integration:** (1) CUDA dhj MHD startup segfault (inline c2p_track lambdas; 952abbfe merged ad992f7b; FofcTrialRemoveGrav same pattern unfixed); (2) cubed-sphere seam flux Isend requests never waited on uniform meshes -> Cray MPICH aborts ~8100 cycles (2dfe4e95 merged 793e03c3; waits only, bitwise). viper HIP/OpenMPI tolerated both; C256 production must use >= 793e03c3. Merged by DeltaAI without a viper combined gate.

**09-30 ~03:15 viper gate of DeltaAI merges PASSED + FOFC fix merged (b6a6eff3, note d87ac849):** dhj hydro, He box, dhj MHD (fofc 0/1) bitwise on HIP vs pre-merge; 15/15 CPU tests. Inline-lambda-in-header fixed: FofcTrialRemoveGrav (fofc_etotgrav.cpp) and HostGamma1FromP (pgen_eos_utils.cpp). Left: multigrid.hpp templates (13 TUs, nvcc risk unclear). CUDA fix untested on GPU (needs an MHD fofc run on DeltaAI/Caltech).
