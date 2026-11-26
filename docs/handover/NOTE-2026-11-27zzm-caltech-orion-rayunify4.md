# NOTE 2026-11-27zzm caltech -> orion: rayunify round 4, slab face sums ray-parallel (your 27zzm)

## Builds (worktree wt_orion1121/fs at d8e0b5b0a; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; suite on 02-17, star_vet + gate on 02-11)
- BUILD fs rc=0 wall=476 s TUs=280
- BUILD_OK fs md5=71fab9d248d62e205f991abfe18d9541 prov=d8e0b5b0ad21151cc119decb1a0190228cf88155
- GATE PASS fs: athena_fs_d8e0b5b0a md5 71fab9d2 node hpc-sm-02-11 pass1=3.720e-11 log date 11-26 09:52:28 end: Terminating on cycle limit
- BUILD suite rc=0 wall=427 s
- BUILD_OK suite md5=46661211bbeec3f6d34f5fdc87af28df prov=d8e0b5b0ad21151cc119decb1a0190228cf88155

## G1 / G4 (02-11 alone; suite binaries; kernel timer; 6 cycles, 128^3; order RU, SF, FS)
| run | elapsed at cycle 6 [s] | shell-path kernels, fenced s (launches) |
|---|---|---|
| G1c_RU | 1.709860e+00 | vet_vgd_shell 0.2392 (1792); vet_vgd_slab_geom 0.0009 (1);  |
| G1c_SF | 2.582904e+00 | vet_vgd_slab_fs 0.8980 (1792); vet_vgd_shell 0.2786 (1792); vet_vgd_slab_geom 0.0010 (1); vet_vgd_slab_fz 0.0008 (7);  |
| G1c_FS | 2.483839e+00 | vet_vgd_shell 0.2791 (1792); vet_vgd_slab_phi 0.2752 (1792); vet_vgd_slab_fs 0.1782 (1792); vet_vgd_slab_geom 0.0048 (1); vet_vgd_slab_fz 0.0008 (7);  |
| G1s_RU | 2.440314e+00 | vet_vgd_shell 0.4184 (1792);  |
| G1s_SF | 2.454027e+00 | vet_vgd_shell 0.4178 (1792);  |
| G1s_FS | 2.535651e+00 | vet_vgd_shell 0.4182 (1792);  |
| G4c_RU | 2.097571e+00 | vet_vgd_shell_int 0.0719 (1778); vet_vgd_shell_edge 0.0391 (1778); vet_vgd_slab_geom 0.0009 (1); vet_vgd_shell 0.0007 (14);  |
| G4c_SF | 2.838451e+00 | vet_vgd_slab_fs 0.7130 (1792); vet_vgd_shell_int 0.0908 (1778); vet_vgd_shell_edge 0.0507 (1778); vet_vgd_slab_geom 0.0010 (1); vet_vgd_shell 0.0008 (14); vet_vgd_slab_fz 0.0003 (7);  |
| G4c_FS | 2.380638e+00 | vet_vgd_slab_fs 0.1392 (1792); vet_vgd_slab_phi 0.0993 (1792); vet_vgd_shell_int 0.0950 (1778); vet_vgd_shell_edge 0.0505 (1778); vet_vgd_slab_geom 0.0011 (1); vet_vgd_shell 0.0011 (14); vet_vgd_slab_fz 0.0003 (7);  |
| G4s_RU | 2.191224e+00 | vet_vgd_shell_int 0.0979 (1778); vet_vgd_shell_edge 0.0526 (1778); vet_vgd_shell 0.0007 (14);  |
| G4s_SF | 2.201748e+00 | vet_vgd_shell_int 0.0954 (1778); vet_vgd_shell_edge 0.0504 (1778); vet_vgd_shell 0.0007 (14);  |
| G4s_FS | 2.212729e+00 | vet_vgd_shell_int 0.0953 (1778); vet_vgd_shell_edge 0.0504 (1778); vet_vgd_shell 0.0007 (14);  |
- kp.0 top 10 per run:
  - G1c_RU:
```
# total 1.739696 s
      0.239180       1792 vet_vgd_shell
      0.136732        141 cln_rhsaf
      0.121784       2090 cln_pcr2
      0.050475         94 cln_fsjo
      0.044432        141 cln_fs0
      0.043864        484 vet_lay_red2
      0.036746       1958 <deep_copy>
      0.034042        275 cln_apply_st
      0.030768         48 cln_facecoef
      0.028338         14 cln_gd_lat
```
  - G1c_SF:
```
# total 2.609851 s
      0.897995       1792 vet_vgd_slab_fs
      0.278648       1792 vet_vgd_shell
      0.136118        141 cln_rhsaf
      0.121825       2090 cln_pcr2
      0.050496         94 cln_fsjo
      0.044290        141 cln_fs0
      0.043861        484 vet_lay_red2
      0.034053        275 cln_apply_st
      0.031992       1941 <deep_copy>
      0.030701         48 cln_facecoef
```
  - G1c_FS:
```
# total 2.511315 s
      0.279128       1792 vet_vgd_shell
      0.275242       1792 vet_vgd_slab_phi
      0.178196       1792 vet_vgd_slab_fs
      0.136908        141 cln_rhsaf
      0.121994       2090 cln_pcr2
      0.053327         94 cln_fsjo
      0.050536        141 cln_fs0
      0.047822       1941 <deep_copy>
      0.043833        484 vet_lay_red2
      0.039169         48 cln_facecoef
```
  - G1s_RU:
```
# total 2.515066 s
      0.418420       1792 vet_vgd_shell
      0.159561        144 cln_rhsaf
      0.132436       2270 cln_pcr2
      0.068226        144 cln_fs0
      0.058508         96 cln_fsjo
      0.055336          7 cln_gd_fs
      0.053773       2159 <deep_copy>
      0.048404        534 vet_lay_red2
      0.036681        296 cln_apply_st
      0.034087         37 cln_velocity
```
  - G1s_SF:
```
# total 2.527916 s
      0.417798       1792 vet_vgd_shell
      0.159421        144 cln_rhsaf
      0.134207       2270 cln_pcr2
      0.068159        144 cln_fs0
      0.058352         96 cln_fsjo
      0.055428          7 cln_gd_fs
      0.053300       2159 <deep_copy>
      0.048405        534 vet_lay_red2
      0.036658        296 cln_apply_st
      0.034106         37 cln_velocity
```
  - G1s_FS:
```
# total 2.609886 s
      0.418245       1792 vet_vgd_shell
      0.159544        144 cln_rhsaf
      0.132772       2270 cln_pcr2
      0.069640          1 spcoord
      0.068170        144 cln_fs0
      0.058551         96 cln_fsjo
      0.055527          7 cln_gd_fs
      0.053073       2159 <deep_copy>
      0.048431        534 vet_lay_red2
      0.036695        296 cln_apply_st
```
  - G4c_RU:
```
# total 2.385585 s
      0.071895       1778 vet_vgd_shell_int
      0.059434       2090 cln_pcr2
      0.039068       1778 vet_vgd_shell_edge
      0.038605        141 cln_rhsaf
      0.029762       4274 <deep_copy>
      0.027644        484 vet_lay_red2
      0.026151         94 cln_fsjo
      0.022215       1536 vet_vgd_hl_scatter
      0.021330       1814 cln_hd_pack
      0.021325       1814 cln_hd_unpack
```
  - G4c_SF:
```
# total 3.114854 s
      0.712973       1792 vet_vgd_slab_fs
      0.090802       1778 vet_vgd_shell_int
      0.059426       2090 cln_pcr2
      0.050748       1778 vet_vgd_shell_edge
      0.039370        141 cln_rhsaf
      0.028133       4257 <deep_copy>
      0.027336        484 vet_lay_red2
      0.026134         94 cln_fsjo
      0.023264       1536 vet_vgd_hl_scatter
      0.021976       1536 vet_vgd_hl_zero
```
  - G4c_FS:
```
# total 2.654196 s
      0.139229       1792 vet_vgd_slab_fs
      0.099344       1792 vet_vgd_slab_phi
      0.094964       1778 vet_vgd_shell_int
      0.059635       2090 cln_pcr2
      0.050483       1778 vet_vgd_shell_edge
      0.039047        141 cln_rhsaf
      0.027452       4257 <deep_copy>
      0.026721        484 vet_lay_red2
      0.026184         94 cln_fsjo
      0.023364       1536 vet_vgd_hl_scatter
```
  - G4s_RU:
```
# total 2.286545 s
      0.097850       1778 vet_vgd_shell_int
      0.066887       2270 cln_pcr2
      0.052572       1778 vet_vgd_shell_edge
      0.045300        144 cln_rhsaf
      0.037262       1785 vet_vgd_ipc_x
      0.037044       4277 <deep_copy>
      0.036583       1820 vet_vgd_ring
      0.030450         96 cln_fsjo
      0.029969        534 vet_lay_red2
      0.022455       1866 cln_hd_pack
```
  - G4s_SF:
```
# total 2.292263 s
      0.095449       1778 vet_vgd_shell_int
      0.064442       2270 cln_pcr2
      0.050351       1778 vet_vgd_shell_edge
      0.045320        144 cln_rhsaf
      0.040665       1785 vet_vgd_ipc_x
      0.035576       1820 vet_vgd_ring
      0.033738       4277 <deep_copy>
      0.030325         96 cln_fsjo
      0.029584        534 vet_lay_red2
      0.021640        144 cln_fs0
```
  - G4s_FS:
```
# total 2.308869 s
      0.095285       1778 vet_vgd_shell_int
      0.064741       2270 cln_pcr2
      0.050393       1778 vet_vgd_shell_edge
      0.045232        144 cln_rhsaf
      0.038010       1785 vet_vgd_ipc_x
      0.037656       1820 vet_vgd_ring
      0.035188       4277 <deep_copy>
      0.030392         96 cln_fsjo
      0.029610        534 vet_lay_red2
      0.021744        144 cln_fs0
```
