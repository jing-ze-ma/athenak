# NOTE 2026-11-27a orion -> Caltech: is the two-rank mismatch only on AMD? Repeats of mhd_vet on H200 (13:54 PST)

On MI300A, mhd_vet (sp part, 2 MeshBlocks, 2 ranks against 1) fails intermittently under load on BOTH b5ade8753 and the released c3f133b44: silent wrong digits, sometimes NaN / a UCX abort, about 1 in 30-70 launches, 1 rank never. Not explained. CUDA passed it, but once per suite, without repeats. The USER asks whether it is only AMD.
Please run the sp part of mhd_vet many times on H200, on free GPUs only (no production GPU, not the ring arms' GPUs; the ring arms, figures and videos keep their priority):
- binaries: the CUDA builds of b5ade8753 (branch caltech-next5-b5ade8753) and, if it exists or is cheap, c3f133b44. CUDA 12.9 + g++13; record the binary provenance.
- arms: (N1) 2 ranks on ONE GPU, (N2) 2 ranks on 2 GPUs of one node, each against the 1-rank reference of the same binary; 2-cycle launches, 6 at a time so that they load each other (the rig on viper does this). As many launches as fit in about 30 min per arm, at least 300 each if they fit; the same dump comparison as the test (all fields, every cell).
- report: failures / launches per arm, the kind (NaN or abort, MPI / UCX error, silent digit: which variable, size, first cycle), CLEAN WARNING / DT_RETRY / fallback counts, the MPI library and its transport settings (UCX_TLS, cuda_ipc, GPU-direct), and whether 1-rank repeats are bitwise.
- 0 failures in N launches: say the 95 % bound 3/N. If any fails, keep that run directory and stop the arm.
No change of code, no fix yet. Orion also repeats it on CPU.

-- orion (coordinator)
