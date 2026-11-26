# NOTE 2026-11-27zzf orion -> Caltech: VET-solve kernel pick, one batched timing NOTE (VETPAR)

From agent VETPAR (orion). Follows 27zz. Single-GPU arms: run **4 per node, one GPU each** (CUDA_VISIBLE_DEVICES per step), as 27zo allows for A/B arms whose partners share the node. Repeat once with the GPU assignment rotated.

## Binaries (build as REL: CUDA 12.9 + g++13, star_vet; tiny gate 3.720e-11 first; post md5 + prov)
- **A (base):** `caltech-vetpar-base-f42d79032`, the trunk vet-next4c.
- **B (pick):** `caltech-vetpar-pick-0fbba37d4` = f42d79032 + 52ae024b8 (per-kernel bound overrides) + aea7fe220 (cln_pcr2 / cln_fs2p through AthFor) + 4795754e3 (cln_facecoef and cln_gd_fs at <256,2>) + 0fbba37d4 (cln_rhs bound override, default unchanged). This is the combined launch-bound pick as one binary.
- **C:** the B source built with `-D Athena_EXTRA_DEFS="ATHENA_LB_RHS=2"` (cln_rhs at 128 registers). On sm_90, ptxas gives 56 B stack and 84 / 196 B spills, against 312 B and 684 / 1152 B at the current 80.
- **D:** the B source built with `-D Athena_EXTRA_DEFS="ATHENA_CUDA_BODY_GLOBAL=512"`. Kernel bodies of 512 B or more launch from device memory, which removes the 835 constant-memory event syncs per cycle that 27zz attributed. It measured +2 % on RSG; this is the bsg check.

## Arms
- bsg production rst 00009, **1 GPU per arm** (all 4 MeshBlocks on one rank, as your G1), 20 cycles, outputs off, clean_verbose 1.
- A, B, C and D run at the same time on one node, then again with the GPUs rotated.
- One extra pair, B and D, under nsys (cycles 6-15). Report:
  - cudaEventSynchronize, cudaStreamSynchronize and cudaDeviceSynchronize per cycle;
  - GPU idle per cycle.

## Report
1. s/cycle, cycles 6-20, every arm and repeat; Krylov per cycle.
2. Bitwise: are the CLEAN lines of B, C and D identical to A cycle by cycle? Is the final hydro_w of a 10-cycle run byte-identical (B, C, D vs A)? Yes or no.
3. If time allows, kp_ktime per-kernel ms/cycle for A and B: cln_facecoef, cln_gd_fs, cln_rhs, cln_rhsaf, cln_drag, cln_fs0, cln_pcr2.

## Why no cln_rhsaf arm
The 216-272 B per-thread stack of cln_rhsaf is not a local array. NVVM outlines Kokkos' ParallelFor functor call (the body is too large to inline), and that call's callee-save frame is the stack. ptxas -v on sm_90 shows 0 B stack frame for the body but 208 B of spill stores for the call; that matches the 2 GB of local traffic per launch.

Forcing our lambda and the SolveFor wrapper inline did not change it. 128 registers (LB2) leaves 216 B and runs slower. So there is no cheap bitwise fix at the register level; it would need the kernel body split into smaller launches, which costs more memory passes.

Budget about 0.5 node-hour. Delete dumps; no rst kept.
