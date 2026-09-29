# NOTE from DeltaAI: dhj MHD segfaults at startup under CUDA (c2p_track) -- fix on branch fix-c2p-track-cuda

**DeltaAI, 2026-09-29 ~19:00 CDT.** Affects every CUDA build (DeltaAI GH200, Caltech H200) of rt-integration since
55e4d001 (09-28, ConToPrim floor bookkeeping) when a dhj run uses `<mhd>`. HIP (viper) is unaffected, which is why the
viper smoke 12027359 passed.

## Symptom
Both w121 MHD scan inputs (b3_e13_f13, b10_e14_sts_f13; binary 350eee1c) died with SIGSEGV right after
`initial B -- ...`, before any output. The crash was the same without ohmic_resistivity, with hlld, without the bin
outputs, and at 1 or 2 ranks (job 3272794). gdb:
```
#0  0x0000000000000000 in ?? ()
#1  __nv_hdl_wrapper_t<... &c2p_track::Before ...>::__nv_hdl_wrapper_t(...)
#2  mhd::MHD::C2PTrack(Driver*, int, bool)
#3  mhd::MHD::ConToPrim(Driver*, int)
#4  Driver::InitBoundaryValuesAndPrimitives
```

## Cause
`c2p_track::Before` / `After` were `inline` functions in `utils/c2p_track.hpp` that contain a `KOKKOS_LAMBDA` (an nvcc
extended host-device lambda). The header is compiled into both hydro_tasks.cpp and mhd_tasks.cpp; the linker keeps one
copy of each inline function, and the MHD call then goes through an uninitialised lambda-wrapper pointer. The hydro
path happened to get a working copy.

## Fix (branch `fix-c2p-track-cuda`, 952abbfe, on 350eee1c; not merged)
Before/After are declared in the header and defined once in the new `src/utils/c2p_track.cpp` (added to
src/CMakeLists.txt). The code is otherwise unchanged. Enable/Weight stay inline (they have no lambda).
- DeltaAI binary athena_dhj_dev_952abbfe6e29, md5 b7fcbbc40366d5cdf659f9631c0c73b5.
- Smoke 3272818 (2 GH200, 300 cycles, fresh -i start): rc 0, 0 FATAL, 0 NaN, 0 NOT-CONVERGED for both arms.
  - b3_e13_f13: dt 1.778 s (Ohmic, as predicted), 16.5 ms/cycle.
  - b10_e14_sts_f13: dt 3.09 -> 3.79 s, 9 RKG stages, 60.7 ms/cycle.
- **Merge decision is the user's.** Please cherry-pick it into rt-integration if the user agrees; Caltech needs it for
  any MHD run.

## Same pattern, not fixed
`FofcTrialRemoveGrav` (src/hydro/fofc_etotgrav.hpp) is an inline header function with a KOKKOS_LAMBDA, used by both
hydro_fofc.cpp and mhd_fofc.cpp. An MHD run with `fofc` on under CUDA will probably crash the same way. The scan has
fofc = 0. The same one-.cpp fix applies.

## Scan status on DeltaAI
- The 12 short arms + noise twin (b3_e12_f13_r2, 1 rank) run with the fixed binary.
- They are 2-GPU jobs on ghx4, with per-arm time limits of 45 min to 3 h (twin: 1 GPU). One arm at a time is moved to
  ghx4-interactive (1 job per user).
- Wall time per arm is about 8-70 min.
