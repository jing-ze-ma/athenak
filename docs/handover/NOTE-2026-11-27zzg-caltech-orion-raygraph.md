# NOTE 2026-11-27zzg caltech -> orion: sp shell loop as one CUDA graph, S4 / S8 against RU (your 27zzg)

## Builds (worktree wt_orion1121/rg at bee5261d2; CUDA 12.9 + g++13, HOPPER90, -ffp-contract=off; 02-17)
- BUILD rg rc=0 wall=573 s TUs=280
- BUILD_OK rg md5=a74acbd8c5ea175659a558cb3886a7ec prov=bee5261d2f3f8a1135d7f221eea11dbaae7ae3a9
- GATE PASS rg: athena_rg_bee5261d2 md5 a74acbd8 node hpc-sm-02-17 pass1=3.720e-11 log date 11-26 02:31:05 end: Terminating on cycle limit
- BUILD suite rc=0 wall=557 s
- BUILD_OK suite md5=e6a374d78c33a97089a6898d4a9c793c prov=bee5261d2f3f8a1135d7f221eea11dbaae7ae3a9

## Arms (27zza deck: tile 16, nlim 120, ndiag 1, cycles 61-120; RU = 0e1601792, RG = bee5261d2; S4 on 02-11 alone, S8 on 02-11 + 02-17; RU then RG)
- S4_RU rc=0 wall=168 s end: Terminating on cycle limit FATAL 0 CW 0
- S4_RG rc=0 wall=165 s end: Terminating on cycle limit FATAL 0 CW 0
- S8_RU rc=0 wall=127 s end: Terminating on cycle limit FATAL 0 CW 0
- S8_RG rc=0 wall=128 s end: Terminating on cycle limit FATAL 0 CW 0
| arm | GPUs | s/cycle (61-120) | zone-cycles/s | rays: source / sweeps (halo+walls, MPI wait) s/cycle | builds/cycle | cut shell | VET solve: clean stage steps s/cycle | Krylov it/cycle | stage - rays | total - stage (hydro + rest) |
|---|---|---|---|---|---|---|---|---|---|---|
| S4_RU | 4 | 1.3066 | 1.284e+07 | 0.0196 / 0.0661 (0.0010, 0.0168) | 1.1 | 135 | 0.7348 | 41.7 | 0.6490 | 0.5718 |
| S4_RG | 4 | 1.2630 | 1.328e+07 | 0.0198 / 0.0679 (0.0001, 0.0172) | 1.1 | 135 | 0.7383 | 41.7 | 0.6506 | 0.5248 |
| S8_RU | 8 | 0.8985 | 1.867e+07 | 0.0259 / 0.1158 (0.1145, 0.1121) | 1.1 | 135 | 0.5786 | 41.7 | 0.4370 | 0.3198 |
| S8_RG | 8 | 0.9172 | 1.829e+07 | 0.0255 / 0.1146 (0.1133, 0.1099) | 1.1 | 135 | 0.5743 | 41.7 | 0.4342 | 0.3430 |
- Fenced ray seconds, stage steps and Krylov counts are run totals / 120.

### Reading
- S4: RG/RU = 0.967 (1.263 vs 1.307 s/cycle). The fenced sweep is the same (0.068 vs 0.066 s/cycle); halo+walls 0.0010 -> 0.0001.
  - The 0.044 s/cycle gain sits outside the stage steps. On this node and this hour the run-to-run spread was about 3 %, so that is at the edge of noise.
- S8: RG/RU = 1.021 (0.917 vs 0.899), with equal sweeps. The graph is **not used** at S8: the log says the MPI halo path launches the loop launch by launch.
- **Bitwise: yes.** The CLEAN lines of RG and RU are identical at S4 and S8 (242/242).
- Graph log lines, verbatim:
  - S4 RG:
    - `<rad_vet> vet_gd graph (rank 0): shell loops captured 1.400000e+01, replayed 2.170000e+02`
    - `<rad_vet> vet_gd graph: shell loop captured as one GPU graph per pass (replayed while its launch inputs are unchanged)`
  - S8 RG:
    - `<rad_vet> vet_gd graph (rank 0): shell loops captured 0.000000e+00, replayed 0.000000e+00`
    - `<rad_vet> vet_gd graph: shell loop launch by launch: MPI halo (host wait per shell)`
  - No FATAL or warning lines about the graph.
- nsys on S4: not done. My wrapper's quoting broke (rc 2 before start). The S4 GPU idle can be added on request.
