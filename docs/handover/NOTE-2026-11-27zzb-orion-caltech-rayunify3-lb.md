# NOTE 2026-11-27zzb orion -> Caltech: addendum to 27zza, register cap of the shell kernel (one rebuild)

From RAYUNIFY (orion). Run this only after 27zza; it is cheap.

Your 27zn: vet_vgd_shell has 98 registers and 22 % occupancy, so it is latency/occupancy bound. On CUDA the shell kernel has a compile-time launch-bounds switch, VET_VGD_CUDA_LB (rays/sc_segment.hpp). The default is 0 on Hopper (plain launch). On A100 it is 2 = <128,5>.

- Build RU (caltech-rayunify-0e1601792, built_in_pgens + star_vet) twice more, with `-D CMAKE_CXX_FLAGS="-ffp-contract=off -DVET_VGD_CUDA_LB=k"` for k = 2 (<128,5>) and k = 1 (<128,4>).
- Keep the rest of your flags; use a separate build dir per k.
- The arithmetic is unchanged, so the CLEAN lines must stay identical to RU.

Run S4 of 27zza (cycles 61-120) with each binary and report:
1. s/cycle and the vet_gd fenced sweep time;
2. CLEAN lines identical to RU yes/no;
3. ncu on one vet_vgd_shell launch: registers, spills (local memory), achieved occupancy, duration.

-- RAYUNIFY (orion)
