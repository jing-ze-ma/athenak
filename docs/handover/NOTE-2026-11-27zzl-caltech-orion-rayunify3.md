# NOTE 2026-11-27zzl caltech -> orion: rayunify round 3, Cartesian slab on the shell kernel (your 27zzl)

## Builds (worktree wt_orion1121/sf at 05ba2d71b; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-17)
- BUILD sf rc=0 wall=701 s TUs=280
- BUILD_OK sf md5=9ac8add536fb9046bb0981a59b051736 prov=05ba2d71b682694fb439bcdcfa0c268f4598fb17
- GATE PASS sf: athena_sf_05ba2d71b md5 9ac8add5 node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 09:15:41 end: Terminating on cycle limit
- BUILD suite rc=0 wall=698 s
- BUILD_OK suite md5=7f128ff4a66fe74b664180ff42e2d16c prov=05ba2d71b682694fb439bcdcfa0c268f4598fb17

- Timing note: this NOTE was committed at 18:02 CET (09:02 PST) and the builds started at 09:03; nothing held it up. The ETA "07:15" I gave the coordinator was my own clock error.

## 1-2. G1 / G4 (02-11 alone; suite binaries; kernel timer; 6 cycles, 128^3; RU then SF)
| run | elapsed at cycle 6 [s] | shell-path kernels, fenced s (launches) |
|---|---|---|
| G1c_RU | 2.120917e+00 | vet_vgd_shell 0.2395 (1792);  |
| G1c_SF | 2.757485e+00 | vet_vgd_slab_fs 0.8982 (1792); vet_vgd_shell 0.2787 (1792);  |
| G1s_RU | 2.492068e+00 | vet_vgd_shell 0.4197 (1792);  |
| G1s_SF | 2.477155e+00 | vet_vgd_shell 0.4178 (1792);  |
| G4c_RU | 2.261397e+00 | vet_vgd_shell_int 0.0724 (1778); vet_vgd_shell_edge 0.0407 (1778); vet_vgd_shell 0.0007 (14);  |
| G4c_SF | 2.969683e+00 | vet_vgd_slab_fs 0.7134 (1792); vet_vgd_shell_int 0.0895 (1778); vet_vgd_shell_edge 0.0495 (1778); vet_vgd_shell 0.0008 (14);  |
| G4s_RU | 2.224421e+00 | vet_vgd_shell_int 0.0966 (1778); vet_vgd_shell_edge 0.0492 (1778); vet_vgd_shell 0.0007 (14);  |
| G4s_SF | 2.211160e+00 | vet_vgd_shell_int 0.0960 (1778); vet_vgd_shell_edge 0.0484 (1778); vet_vgd_shell 0.0007 (14);  |
- kp.0 top 10 per run:
  - G1c_RU:
```
# total 2.151493 s
      0.239477       1792 vet_vgd_shell
      0.136723        141 cln_rhsaf
      0.122199       2090 cln_pcr2
      0.079764       1958 <deep_copy>
      0.051941         94 cln_fsjo
      0.049060        141 cln_fs0
      0.046524         12 hflux_x1
      0.043945        484 vet_lay_red2
      0.036612         42 Kokkos::ViewFill-4D
      0.035631         48 cln_facecoef
```
  - G1c_SF:
```
# total 2.784805 s
      0.898211       1792 vet_vgd_slab_fs
      0.278679       1792 vet_vgd_shell
      0.135953        141 cln_rhsaf
      0.121907       2090 cln_pcr2
      0.052208         94 cln_fsjo
      0.049015        141 cln_fs0
      0.044358       1941 <deep_copy>
      0.043906        484 vet_lay_red2
      0.036265         48 cln_facecoef
      0.034041        275 cln_apply_st
```
  - G1s_RU:
```
# total 2.566194 s
      0.419702       1792 vet_vgd_shell
      0.159457        144 cln_rhsaf
      0.130603       2270 cln_pcr2
      0.068265        144 cln_fs0
      0.058477         96 cln_fsjo
      0.055392          7 cln_gd_fs
      0.053806       2159 <deep_copy>
      0.048216        534 vet_lay_red2
      0.036557        296 cln_apply_st
      0.034128         37 cln_velocity
```
  - G1s_SF:
```
# total 2.551311 s
      0.417804       1792 vet_vgd_shell
      0.159377        144 cln_rhsaf
      0.132578       2270 cln_pcr2
      0.068200        144 cln_fs0
      0.058418         96 cln_fsjo
      0.055346          7 cln_gd_fs
      0.054132       2159 <deep_copy>
      0.048377        534 vet_lay_red2
      0.036652        296 cln_apply_st
      0.034123         37 cln_velocity
```
  - G4c_RU:
```
# total 2.535296 s
      0.072362       1778 vet_vgd_shell_int
      0.061701       2090 cln_pcr2
      0.040681       1778 vet_vgd_shell_edge
      0.038954        141 cln_rhsaf
      0.031800       4274 <deep_copy>
      0.027498        484 vet_lay_red2
      0.026234         94 cln_fsjo
      0.023294       1536 vet_vgd_hl_scatter
      0.022728       1814 cln_hd_unpack
      0.022354       1814 cln_hd_pack
```
  - G4c_SF:
```
# total 3.256542 s
      0.713428       1792 vet_vgd_slab_fs
      0.089532       1778 vet_vgd_shell_int
      0.059629       2090 cln_pcr2
      0.049451       1778 vet_vgd_shell_edge
      0.039318        141 cln_rhsaf
      0.028957       4257 <deep_copy>
      0.026801        484 vet_lay_red2
      0.026200         94 cln_fsjo
      0.023700       1536 vet_vgd_hl_scatter
      0.022691        768 vet_vgd_hl_scan2
```
  - G4s_RU:
```
# total 2.319797 s
      0.096616       1778 vet_vgd_shell_int
      0.064929       2270 cln_pcr2
      0.049201       1778 vet_vgd_shell_edge
      0.045682        144 cln_rhsaf
      0.041334       1785 vet_vgd_ipc_x
      0.040104       1820 vet_vgd_ring
      0.034831       4277 <deep_copy>
      0.030272         96 cln_fsjo
      0.029801        534 vet_lay_red2
      0.022552        382 vet_vgd_hxb_mark
```
  - G4s_SF:
```
# total 2.310512 s
      0.095963       1778 vet_vgd_shell_int
      0.066880       2270 cln_pcr2
      0.048397       1778 vet_vgd_shell_edge
      0.045674        144 cln_rhsaf
      0.042063       1785 vet_vgd_ipc_x
      0.037701       1820 vet_vgd_ring
      0.037577       4277 <deep_copy>
      0.030353         96 cln_fsjo
      0.029900        534 vet_lay_red2
      0.022819       1866 cln_hd_unpack
```
- **Reading:**
  - cart2, 1 GPU: SF's slab shell is 0.279 s (RU's slab 0.240 s), but the new per-shell face-sum kernel vet_vgd_slab_fs adds **0.898 s**. The Cartesian path is about 4.9x the shell time of RU (1.18 vs 0.24 s), and the run's elapsed is 2.76 vs 2.12 s.
  - cart2, 4 GPUs: slab_fs 0.713 s, against 0.113 s for shell_int + edge.
  - sp2, 1 and 4 GPUs: identical to RU (0.418 vs 0.420 s; 0.145 vs 0.146 s). As intended, the sp kernel has no slab code.
  - So the regression has moved into vet_vgd_slab_fs (1792 launches, about 0.5 ms each at 128^3 on 1 GPU).

## 3. S4 / S8 (27zy bsg decks, cycles 61-120; star_vet binaries)
| arm | GPUs | s/cycle (61-120) | zone-cycles/s | rays: source / sweeps (halo+walls, MPI wait) s/cycle | builds/cycle | cut shell | VET solve: clean stage steps s/cycle | Krylov it/cycle | stage - rays | total - stage (hydro + rest) |
|---|---|---|---|---|---|---|---|---|---|---|
| S4_RU | 4 | 1.2402 | 1.353e+07 | 0.0192 / 0.0663 (0.0010, 0.0167) | 1.1 | 135 | 0.7330 | 41.7 | 0.6476 | 0.5071 |
| S4_SF | 4 | 1.2718 | 1.319e+07 | 0.0208 / 0.0659 (0.0011, 0.0170) | 1.1 | 135 | 0.7340 | 41.7 | 0.6473 | 0.5378 |
| S8_RU | 8 | 0.8770 | 1.913e+07 | 0.0261 / 0.1131 (0.1119, 0.1090) | 1.1 | 135 | 0.5755 | 41.7 | 0.4363 | 0.3016 |
| S8_SF | 8 | 0.8597 | 1.951e+07 | 0.0250 / 0.1126 (0.1113, 0.1081) | 1.1 | 135 | 0.5719 | 41.7 | 0.4343 | 0.2878 |
- S4: SF/RU = 1.025; S8: 0.980. Both inside the 2-3 % run-to-run spread; the sweep line is equal (0.066 / 0.113 s per cycle).
- **Bitwise: yes.** The CLEAN lines of SF and RU are identical at S4 and S8 (242/242).
- Files: docs/handover/caltech_results/RAYUNIFY3/<run>/.
