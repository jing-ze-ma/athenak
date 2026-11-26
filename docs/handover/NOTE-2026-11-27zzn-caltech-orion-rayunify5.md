# NOTE 2026-11-27zzn caltech -> orion: rayunify round 5, FD / FD0 against FS (your 27zzn)

## Builds (worktree wt_orion1121/fd at c9d5d3ee8; CUDA 12.9 + g++13, HOPPER90; built_in_pgens only; 02-17)
- BUILD fd rc=0 wall=671 s
- BUILD_OK fd (built_in_pgens, CMAKE_CXX_FLAGS '-ffp-contract=off') md5=dfcfa8081bcd3eac9c7107b23511d959 prov=c9d5d3ee8b79e5315cd8e365f80f723afd0a7794
- BUILD fd0 rc=0 wall=422 s
- BUILD_OK fd0 (built_in_pgens, CMAKE_CXX_FLAGS '-ffp-contract=off -DVET_SLAB_PHI_FOLD=0') md5=7bde3f99075fde452f6b8359f3dfc581 prov=c9d5d3ee8b79e5315cd8e365f80f723afd0a7794

## G1 / G4 (02-11 alone; suite binaries; kernel timer; 6 cycles, 128^3; order FS, FD, FD0)
| run | elapsed at cycle 6 [s] | peak GPU memory used (nvidia-smi, max over GPUs, MiB) | shell-path kernels, fenced s (launches) |
|---|---|---|---|
| G1c_FS | 2.138930e+00 |  | vet_vgd_shell 0.2791 (1792); vet_vgd_slab_phi 0.2751 (1792); vet_vgd_slab_fs 0.1763 (1792); vet_vgd_slab_geom 0.0011 (1); vet_vgd_slab_fz 0.0008 (7);  |
| G1c_FD | 8.034114e+00 |  | vet_vgd_shell 6.1332 (1792); vet_vgd_slab_fs 0.2101 (14); vet_vgd_slab_geom 0.0087 (1); vet_vgd_slab_fz 0.0009 (7);  |
| G1c_FD0 | 2.573892e+00 |  | vet_vgd_shell 0.4333 (1792); vet_vgd_slab_phi 0.3459 (1792); vet_vgd_slab_fs 0.2102 (14); vet_vgd_slab_geom 0.0063 (1); vet_vgd_slab_fz 0.0008 (7);  |
| G1s_FS | 2.456773e+00 |  | vet_vgd_shell 0.4180 (1792);  |
| G1s_FD | 2.488171e+00 |  | vet_vgd_shell 0.4179 (1792);  |
| G1s_FD0 | 2.471759e+00 |  | vet_vgd_shell 0.4195 (1792);  |
| G4c_FS | 2.374623e+00 |  | vet_vgd_slab_fs 0.1393 (1792); vet_vgd_slab_phi 0.0994 (1792); vet_vgd_shell_int 0.0932 (1778); vet_vgd_shell_edge 0.0505 (1778); vet_vgd_slab_geom 0.0011 (1); vet_vgd_shell 0.0008 (14); vet_vgd_slab_fz 0.0003 (7);  |
| G4c_FD | 3.428457e+00 |  | vet_vgd_shell_int 1.0773 (1778); vet_vgd_shell_edge 0.2815 (1778); vet_vgd_slab_fs 0.0534 (14); vet_vgd_shell 0.0025 (14); vet_vgd_slab_geom 0.0010 (1); vet_vgd_slab_fz 0.0003 (7);  |
| G4c_FD0 | 2.364858e+00 |  | vet_vgd_shell_int 0.1404 (1778); vet_vgd_slab_phi 0.1157 (1792); vet_vgd_slab_fs 0.0545 (14); vet_vgd_shell_edge 0.0519 (1778); vet_vgd_shell 0.0013 (14); vet_vgd_slab_geom 0.0011 (1); vet_vgd_slab_fz 0.0003 (7);  |
| G4s_FS | 2.174440e+00 |  | vet_vgd_shell_int 0.0963 (1778); vet_vgd_shell_edge 0.0514 (1778); vet_vgd_shell 0.0007 (14);  |
| G4s_FD | 2.216818e+00 |  | vet_vgd_shell_int 0.0965 (1778); vet_vgd_shell_edge 0.0513 (1778); vet_vgd_shell 0.0008 (14);  |
| G4s_FD0 | 2.182352e+00 |  | vet_vgd_shell_int 0.0957 (1778); vet_vgd_shell_edge 0.0506 (1778); vet_vgd_shell 0.0008 (14);  |
- Peak memory: nvidia-smi sampled every 0.2 s during the run (whole GPU, includes the CUDA context of about 0.5 GB).
- kp.0 top 10 per run:
  - G1c_FS:
```
# total 2.166268 s
      0.279058       1792 vet_vgd_shell
      0.275087       1792 vet_vgd_slab_phi
      0.176299       1792 vet_vgd_slab_fs
      0.136968        141 cln_rhsaf
      0.121835       2090 cln_pcr2
      0.050500         94 cln_fsjo
      0.044334        141 cln_fs0
      0.043889        484 vet_lay_red2
      0.034078        275 cln_apply_st
      0.031678       1941 <deep_copy>
```
  - G1c_FD:
```
# total 8.069717 s
      6.133245       1792 vet_vgd_shell
      0.210141         14 vet_vgd_slab_fs
      0.136468        141 cln_rhsaf
      0.121858       2090 cln_pcr2
      0.064989       1941 <deep_copy>
      0.052073         94 cln_fsjo
      0.049383        141 cln_fs0
      0.043808        484 vet_lay_red2
      0.041942          2 vet_vgd_shl
      0.039846         12 hflux_x1
```
  - G1c_FD0:
```
# total 2.607248 s
      0.433315       1792 vet_vgd_shell
      0.345858       1792 vet_vgd_slab_phi
      0.210231         14 vet_vgd_slab_fs
      0.136615        141 cln_rhsaf
      0.120122       2090 cln_pcr2
      0.052152         94 cln_fsjo
      0.048742        141 cln_fs0
      0.044229       1941 <deep_copy>
      0.043598        484 vet_lay_red2
      0.036721         48 cln_facecoef
```
  - G1s_FS:
```
# total 2.532728 s
      0.417982       1792 vet_vgd_shell
      0.159427        144 cln_rhsaf
      0.132419       2270 cln_pcr2
      0.068152        144 cln_fs0
      0.060168         96 cln_fsjo
      0.055297          7 cln_gd_fs
      0.053702       2159 <deep_copy>
      0.048399        534 vet_lay_red2
      0.036686        296 cln_apply_st
      0.034119         37 cln_velocity
```
  - G1s_FD:
```
# total 2.563241 s
      0.417904       1792 vet_vgd_shell
      0.160208        144 cln_rhsaf
      0.132541       2270 cln_pcr2
      0.068121        144 cln_fs0
      0.058394         96 cln_fsjo
      0.055375          7 cln_gd_fs
      0.053659       2159 <deep_copy>
      0.048402        534 vet_lay_red2
      0.036671        296 cln_apply_st
      0.034076         37 cln_velocity
```
  - G1s_FD0:
```
# total 2.545035 s
      0.419469       1792 vet_vgd_shell
      0.159740        144 cln_rhsaf
      0.132597       2270 cln_pcr2
      0.068211        144 cln_fs0
      0.058463         96 cln_fsjo
      0.055313          7 cln_gd_fs
      0.054190       2159 <deep_copy>
      0.048422        534 vet_lay_red2
      0.036670        296 cln_apply_st
      0.034118         37 cln_velocity
```
  - G4c_FS:
```
# total 2.647525 s
      0.139261       1792 vet_vgd_slab_fs
      0.099391       1792 vet_vgd_slab_phi
      0.093208       1778 vet_vgd_shell_int
      0.059508       2090 cln_pcr2
      0.050500       1778 vet_vgd_shell_edge
      0.039018        141 cln_rhsaf
      0.027959       4257 <deep_copy>
      0.026754        484 vet_lay_red2
      0.026206         94 cln_fsjo
      0.023339       1536 vet_vgd_hl_scatter
```
  - G4c_FD:
```
# total 3.711560 s
      1.077325       1778 vet_vgd_shell_int
      0.281481       1778 vet_vgd_shell_edge
      0.061651       2090 cln_pcr2
      0.053427         14 vet_vgd_slab_fs
      0.039031        141 cln_rhsaf
      0.029012       4257 <deep_copy>
      0.027739       1536 vet_vgd_hl_zero
      0.026945          2 vet_vgd_shl
      0.026852        484 vet_lay_red2
      0.026177         94 cln_fsjo
```
  - G4c_FD0:
```
# total 2.644230 s
      0.140387       1778 vet_vgd_shell_int
      0.115714       1792 vet_vgd_slab_phi
      0.059675       2090 cln_pcr2
      0.054501         14 vet_vgd_slab_fs
      0.051914       1778 vet_vgd_shell_edge
      0.038973        141 cln_rhsaf
      0.030028        484 vet_lay_red2
      0.027922       4257 <deep_copy>
      0.026187         94 cln_fsjo
      0.024841       1536 vet_vgd_hl_zero
```
  - G4s_FS:
```
# total 2.275785 s
      0.096336       1778 vet_vgd_shell_int
      0.066821       2270 cln_pcr2
      0.051400       1778 vet_vgd_shell_edge
      0.045336        144 cln_rhsaf
      0.039086       1785 vet_vgd_ipc_x
      0.036555       1820 vet_vgd_ring
      0.035742       4277 <deep_copy>
      0.032435         96 cln_fsjo
      0.029838        534 vet_lay_red2
      0.022529       1866 cln_hd_pack
```
  - G4s_FD:
```
# total 2.319724 s
      0.096485       1778 vet_vgd_shell_int
      0.066271       2270 cln_pcr2
      0.051337       1778 vet_vgd_shell_edge
      0.045385        144 cln_rhsaf
      0.040338       1785 vet_vgd_ipc_x
      0.037621       1820 vet_vgd_ring
      0.037520       4277 <deep_copy>
      0.030781         96 cln_fsjo
      0.030224        534 vet_lay_red2
      0.026722          2 vet_vgd_shl
```
  - G4s_FD0:
```
# total 2.278913 s
      0.095717       1778 vet_vgd_shell_int
      0.064778       2270 cln_pcr2
      0.050552       1778 vet_vgd_shell_edge
      0.045290        144 cln_rhsaf
      0.040280       1785 vet_vgd_ipc_x
      0.037078       1820 vet_vgd_ring
      0.036172       4277 <deep_copy>
      0.030372         96 cln_fsjo
      0.029601        534 vet_lay_red2
      0.021624        144 cln_fs0
```

## Peak GPU memory, G1 cart2 (128^3, 1 GPU)
- Measured by nvidia-smi every 0.1 s inside the run's own step; whole GPU, including the CUDA context. The sampler of the table above ran as a separate step and saw no device, so its column is empty.
  - FS 13029 MiB, FD 22953 MiB, FD0 22953 MiB.
  - So FD / FD0 hold about 9.9 GB more than FS, which matches the expected scratch of about 10.6 GB.

## Reading
- **FD (fold on) regresses badly on cart2.** vet_vgd_shell takes 6.13 s on 1 GPU (FS 0.279 + slab_phi 0.275 + slab_fs 0.176 = 0.73 s); elapsed 8.03 vs 2.14 s. On 4 GPUs shell_int 1.08 + edge 0.28 s (FS 0.38 s in total).
- FD0 (fold off) is close to FS: on 1 GPU shell 0.433 + slab_phi 0.346 + slab_fs 0.210 = 0.99 s (FS 0.73); on 4 GPUs 0.36 s (FS 0.38).
  - slab_fs now launches once per pass (14 launches instead of 1792): 0.21 s on 1 GPU, 0.05 s on 4.
- sp2: identical in all three (0.418 s on 1 GPU; 0.146-0.148 s on 4 GPUs).
