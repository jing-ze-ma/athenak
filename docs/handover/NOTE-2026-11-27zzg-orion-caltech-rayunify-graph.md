# NOTE 2026-11-27zzg orion -> Caltech: addendum to 27zza, the sp sweep as one CUDA graph (one more binary)

From RAYUNIFY (orion). Run after 27zza (and 27zzb if queued). Idle nodes only.

**RG** = branch **`caltech-rayunify-bee5261d2`** = RU (0e1601792) + one commit. Each pass of the sp shell loop (shell kernels, ring copies, IPC exchange kernels) is captured once as a CUDA graph and replayed while its inputs hold. Configurations it cannot capture (MPI halo with host waits) run launch by launch, and the log says so once.
- Checked so far: HIP bitwise vs RU (1 rank; 2 ranks raw band + IPC + ring); CPU suite quick + rays bitwise; nvcc 12.1 compile-checked.
- Not tested: the CUDA-only parts: the fused vet_vgd_ipc_x with the device-side sequence number, and the second-stream interior fork inside the capture.
- On MI300A it gave no gain: the stream stayed full. H200 at 4 GPUs is where the launch gaps were (your 27ze kp_H4_B: ~23 ms/build of gaps).

Build like RU (CUDA 12.9 + g++13, star_vet + built_in_pgens; tiny gate pass1 3.720e-11). Run S4 and S8 of 27zza with RG; RU from 27zza is the reference.

Report:
1. s/cycle (61-120) and the vet_gd fenced sweep line, RG vs RU.
2. **Bitwise:** are the CLEAN lines of RG and RU identical cycle by cycle (S4 and S8)? Any FATAL / graph line in the log (grep -i graph), verbatim.
3. kp is not meaningful with graphs (the tool fences each kernel; graphs hide them), so give s/cycle only. If you have nsys at hand: GPU idle per cycle on rank 0, RU vs RG, S4.

-- RAYUNIFY (orion)
