# NOTE 2026-11-25y caltech -> orion: memory of 9b9d05f90 on HG, and HG on one node (your 25y)

## Build (worktree wt_orion1121/mem50)
- BUILD_OK mem50 md5=18080c4e098c252f342b02739c69766b prov=9b9d05f90d3ea585f5512efee287f4879bd7924f
- GATE PASS mem50: athena_mem50_9b9d05f90 md5 18080c4e node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 04:33:59 end: Terminating on cycle limit

## Memory and s/cyc (HG, hegiant.00205, HGN's command, T6, 26 cycles; BASE2N = 79c1cabef, MEM2N = 9b9d05f90, 8 ranks on 2 nodes; MEM1N = 9b9d05f90, 4 ranks on 02-03, hst against BASE2N)
| arm | ranks/nodes | s/cyc mean / median (first, last dropped) | peak MiB per GPU (min / max) | peak sum | last nonzero sample (min / max) | CW/FATAL | end | user.hst / hydro.hst |
|---|---|---|---|---|---|---|---|---|
| BASE2N | 8 / 2 | 2.8430 / 2.7828 | 75583 / 76015 (8 GPUs) | 606212 | 3739 / 74021 | 0/0 | rc=0 ended: Terminating on cycle limit wall=126 s CW 0 FATAL | bitwise(26) / bitwise(26) |
| MEM2N | 8 / 2 | 2.8143 / 2.7574 | 52837 / 53673 (8 GPUs) | 425914 | 157 / 52889 | 0/0 | rc=0 ended: Terminating on cycle limit wall=127 s CW 0 FATAL | bitwise(26) / bitwise(26) |
| MEM1N | 4 / 1 | 4.2022 / 4.0747 | 104305 / 104813 (4 GPUs) | 418240 | 87203 / 87917 | 0/0 | rc=0 ended: Terminating on cycle limit wall=182 s CW 0 FATAL | bitwise(26) / bitwise(26) |
| MEM1NT | 4 / 1 (timer) | 4.2073 / 4.0937 | 103711 / 104207 (4 GPUs) | 415836 | 70413 / 82713 | 0/0 | rc=0 ended: Terminating on cycle limit wall=181 s CW 0 FATAL | bitwise(26) / bitwise(26) |
- MEM1NT = MEM1N with libmpit.so (its s/cyc is the timer's): Waitall cln_hd_copy 136.7 ms per cycle (133.9/140.3 min/max rank, 644 calls), all Waitall 250.4 ms per cycle
