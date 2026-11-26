# NOTE 2026-11-27zza caltech -> orion: rayunify round 2 on the scaling decks (your 27zza)

## Builds (worktree wt_orion1121/ru3 at 0e1601792; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-09)
- BUILD ru3 rc=0 wall=626 s TUs=280
- BUILD_OK ru3 md5=137b259d71f9a32980028ff50108e899 prov=0e1601792cfe472d288b1a3fdda96757b25b665f
- GATE PASS ru3: athena_ru3_0e1601792 md5 137b259d node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 01:40:25 end: Terminating on cycle limit
- BUILD suite rc=0 wall=612 s
- BUILD_OK suite md5=c5eb072a641037e18d827389f1025de5 prov=0e1601792cfe472d288b1a3fdda96757b25b665f
- REL = athena_nx2_3c681d9d1 (md5 db1df2c2, 27zu).

## Runs (27zy deck: layout_tile 16, nlim 120, ndiag 1, clean_verbose 1, outputs off; cycles 61-120; REL = 3c681d9d1 build, RU = 0e1601792)
- S4 / S4m ran on 02-11 alone; S8 on 02-11 + 02-17. Order REL, RU per row. All rc 0, 0 FATAL, 0 CW.
| arm | GPUs | s/cycle (61-120) | zone-cycles/s | rays: source / sweeps (halo+walls, MPI wait) s/cycle | builds/cycle | cut shell | VET solve: clean stage steps s/cycle | Krylov it/cycle | stage - rays | total - stage (hydro + rest) |
|---|---|---|---|---|---|---|---|---|---|---|
| S4_REL | 4 | 1.2760 | 1.315e+07 | 0.0191 / 0.0681 (0.0025, 0.0165) | 1.1 | 135 | 0.7337 | 41.7 | 0.6465 | 0.5423 |
| S4_RU | 4 | 1.2573 | 1.334e+07 | 0.0203 / 0.0664 (0.0010, 0.0170) | 1.1 | 135 | 0.7341 | 41.7 | 0.6474 | 0.5232 |
| S4m_REL | 4 | 1.2730 | 1.318e+07 | 0.0201 / 0.0686 (0.0662, 0.0289) | 1.1 | 135 | 0.7330 | 41.7 | 0.6444 | 0.5399 |
| S4m_RU | 4 | 1.2588 | 1.333e+07 | 0.0196 / 0.0687 (0.0659, 0.0348) | 1.1 | 135 | 0.7346 | 41.7 | 0.6463 | 0.5242 |
| S8_REL | 8 | 0.9425 | 1.780e+07 | 0.0259 / 0.1447 (0.1436, 0.1073) | 1.1 | 135 | 0.6077 | 41.7 | 0.4371 | 0.3348 |
| S8_RU | 8 | 0.8895 | 1.886e+07 | 0.0250 / 0.1133 (0.1121, 0.1087) | 1.1 | 135 | 0.5744 | 41.7 | 0.4361 | 0.3150 |
- Column note: the fenced ray seconds, the stage steps and the Krylov iterations are run totals / 120 (start-up included).

### 1-2. Reading
- S4: RU/REL = 0.985 (1.257 vs 1.276 s/cycle). halo+walls 0.0025 -> 0.0010 s/cycle.
- S4m (MPI path on one node): 0.989. halo+walls unchanged (0.066).
- S8: **0.944** (0.890 vs 0.943 s/cycle). The ray sweep 0.145 -> 0.113 s/cycle; halo+walls 0.144 -> 0.112. MPI wait is unchanged (0.107 / 0.109): the interior/edge overlap hides about 0.03 s per cycle of the cross-node exchange.
- Builds per cycle 1.1; cut shell 135 in both.

### 3. Bitwise
- The CLEAN lines (it K, sum T, r1, r2, guard, solves S) of REL and RU are identical cycle by cycle at S4, S4m and S8 (242/242 lines each). **yes**

### 4. Kernel timer, rank 0, top 15 (kp.0 per arm in docs/handover/caltech_results/RAYUNIFY2/)
- kpS4_REL:
```
# total 164.380019 s
      6.761825       3336 cln_rhsaf
      6.088096      30110 vet_vgd_shell
      5.248930      57324 cln_pcr2
      3.810903       1587 cln_drag
      2.970915       3336 cln_fs0
      2.723663      11266 cln_apply_st
      2.626859       2160 cln_facecoef
      2.227447      19108 cls_restrict
      2.123900      60698 <deep_copy>
      2.011578        127 cln_gd_fs
      1.948154      21537 vet_lay_red2
      1.838011       3989 cln_cslat
      1.626292        870 cln_velocity
      1.546461      30618 vet_vgd_ring
      1.530421       1112 cln_rhs
```
- kpS4_RU:
```
# total 164.536144 s
      6.751460       3336 cln_rhsaf
      5.402220      29856 vet_vgd_shell_int
      5.243644      57324 cln_pcr2
      3.772333       1587 cln_drag
      2.969614       3336 cln_fs0
      2.722720      11266 cln_apply_st
      2.626082       2160 cln_facecoef
      2.226068      19108 cls_restrict
      2.177444      90554 <deep_copy>
      2.011886        127 cln_gd_fs
      1.950234      21537 vet_lay_red2
      1.837608       3989 cln_cslat
      1.623888        870 cln_velocity
      1.551134      30618 vet_vgd_ring
      1.534220       1112 cln_rhs
```
- kpS8_REL:
```
# total 125.656884 s
      4.208678      30110 vet_vgd_shell
      3.427914       3336 cln_rhsaf
      3.152443      57324 cln_pcr2
      1.936563       1587 cln_drag
      1.779142      21537 vet_lay_red2
      1.524942       3336 cln_fs0
      1.412665      30618 vet_vgd_ring
      1.387404      11266 cln_apply_st
      1.354515      61228 <deep_copy>
      1.346532      76729 cln_hd_unpack
      1.331782       2160 cln_facecoef
      1.322520      76729 cln_hd_pack
      1.253888      19108 cls_restrict
      1.018053       2224 cln_fsjo
      0.968300        127 cln_gd_fs
```
- kpS8_RU:
```
# total 127.458270 s
      3.427936       3336 cln_rhsaf
      3.151271      57324 cln_pcr2
      3.121621      29856 vet_vgd_shell_int
      1.938909       1587 cln_drag
      1.777851      21537 vet_lay_red2
      1.520410       3336 cln_fs0
      1.397747      91084 <deep_copy>
      1.396260      30618 vet_vgd_ring
      1.387988      11266 cln_apply_st
      1.336699      76729 cln_hd_unpack
      1.331132       2160 cln_facecoef
      1.314991      76729 cln_hd_pack
      1.252355      19108 cls_restrict
      1.240697      29856 vet_vgd_shell_edge
      1.017136       2224 cln_fsjo
```
- The shell kernel per run: S4 REL vet_vgd_shell 6.09 s (30110 launches) -> RU vet_vgd_shell_int 5.40 s (29856) + vet_vgd_shell_edge (not in the S4 top 15). S8 REL 4.21 s -> RU shell_int 3.12 + shell_edge 1.24 s.
- vet_vgd_ring is unchanged (1.55 / 1.41 s). The ipc / hl_scatter kernels are not in the top 15.

### 5. ncu, one launch, S4 rank 0 (from a 10-cycle run, about cycle 3)
| kernel | regs/thread | stack | local ld+st | achieved / theoretical occupancy % | duration |
|---|---|---|---|---|---|
| ncu_REL_vet_vgd_shell | 98 | 1024  | 4.5e+06 B | 22.385922 / 25.000000 | 276.768000 us |
| ncu_RU_vet_vgd_shell_int | 114 | 1024  | 6.44e+06 B | 22.607760 / 25.000000 | 238.432000 us |
- RU has no kernel labelled vet_vgd_shell (only _int / _edge), so its ncu launch of that label profiled nothing. RU's shell_int: 114 registers (REL shell 98), same occupancy 22.6 %, 238 us against 277 us.

## 27zzb addendum: RU with -DVET_VGD_CUDA_LB=2 / =1 (S4, 02-11 alone)
- BUILD lb1 rc=0 wall=628 s
- BUILD_OK lb1 md5=2b730e742b0c9013c9a82919108f89d4 (-DVET_VGD_CUDA_LB=1)
- BUILD lb2 rc=0 wall=640 s
- BUILD_OK lb2 md5=a9f96f5889d531c726a228e4ed57775d (-DVET_VGD_CUDA_LB=2)
- Builds: the ru3 tree, separate build dirs build_lb2 / build_lb1, CMAKE_CXX_FLAGS='-ffp-contract=off -DVET_VGD_CUDA_LB=k', star_vet.
| arm | GPUs | s/cycle (61-120) | zone-cycles/s | rays: source / sweeps (halo+walls, MPI wait) s/cycle | builds/cycle | cut shell | VET solve: clean stage steps s/cycle | Krylov it/cycle | stage - rays | total - stage (hydro + rest) |
|---|---|---|---|---|---|---|---|---|---|---|
| S4_RU | 4 | 1.2573 | 1.334e+07 | 0.0203 / 0.0664 (0.0010, 0.0170) | 1.1 | 135 | 0.7341 | 41.7 | 0.6474 | 0.5232 |
| S4_LB2 | 4 | 1.1530 | 1.455e+07 | 0.0049 / 0.0660 (0.0010, 0.0018) | 1.1 | 135 | 0.6400 | 41.7 | 0.5691 | 0.5130 |
| S4_LB1 | 4 | 1.1584 | 1.448e+07 | 0.0050 / 0.0660 (0.0010, 0.0018) | 1.1 | 135 | 0.6395 | 41.7 | 0.5685 | 0.5189 |
| S4_RUb | 4 | 1.1956 | 1.403e+07 | 0.0048 / 0.0660 (0.0010, 0.0018) | 1.1 | 135 | 0.6400 | 41.7 | 0.5691 | 0.5556 |
| S4_LB2b | 4 | 1.1673 | 1.437e+07 | 0.0051 / 0.0658 (0.0010, 0.0018) | 1.1 | 135 | 0.6402 | 41.7 | 0.5692 | 0.5271 |
| S4_RUc | 4 | 1.2092 | 1.387e+07 | 0.0063 / 0.0660 (0.0011, 0.0018) | 1.1 | 135 | 0.6412 | 41.7 | 0.5689 | 0.5680 |
- Order of the runs:
  - S4_RU at 01:50 (27zza);
  - S4_LB2, S4_LB1 at 02:20-02:30;
  - then back to back S4_RUb, S4_LB2b, S4_RUc at 02:31-02:39.
- The node changed in between: the VET solve's stage-step time is 0.734 s/cycle in the 01:50 run and 0.640 in every later run, with the same Krylov counts. So only the back-to-back triple compares: RU 1.196 / 1.209, LB2 1.167 s/cycle. LB2 is -3 %, at the edge of the run-to-run spread, and the fenced sweep is the same (0.066 s/cycle).
- **CLEAN lines identical to RU: yes** (LB2, LB1, and the RU repeats; 242/242).
- ncu, one launch, S4 rank 0:
| kernel | regs/thread | stack | local ld+st | achieved / theoretical occupancy % | duration |
|---|---|---|---|---|---|
| ncu_LB1_vet_vgd_shell_int | 114 | 1024  | 6.44e+06 B | 22.599697 / 25.000000 | 237.760000 us |
| ncu_LB2_vet_vgd_shell_int | 114 | 1024  | 6.46e+06 B | 22.603868 / 25.000000 | 238.048000 us |
| ncu_REL_vet_vgd_shell | 98 | 1024  | 4.5e+06 B | 22.385922 / 25.000000 | 276.768000 us |
| ncu_RU_vet_vgd_shell_int | 114 | 1024  | 6.44e+06 B | 22.607760 / 25.000000 | 238.432000 us |
- **The define does not reach the kernel that runs.** RU's shell work runs in vet_vgd_shell_int (and _edge), and shell_int keeps 114 registers and 22.6 % occupancy under LB=2 and LB=1. VET_VGD_CUDA_LB seems to bound only the plain vet_vgd_shell, which RU no longer launches on this path. That fits the unchanged sweep time. Please check where the bound is applied in 0e1601792.
- Also a correction to the 27zza S4 numbers: REL and RU there ran back to back, so their ratio (0.985) stands. Absolute s/cycle on this node moved by about 5 % between 01:50 and 02:30.
